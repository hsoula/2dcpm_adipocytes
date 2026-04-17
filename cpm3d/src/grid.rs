use std::fs;
use rand::prelude::*;
use rand_distr::{Distribution, Normal};
use rand_distr::num_traits::Pow;
use serde::{Deserialize, Serialize};

use crate::cellstate::CellState;
use crate::energy::{j, delta_volume_loss, delta_volume_gain, delta_surface, delta_sphericity};
use crate::init::{place_cells_spheres, place_cells_spheres_individual};
use crate::params::Params;

// ── Neighbour offsets ─────────────────────────────────────────────────────────

/// 26-connected Moore neighbourhood (copy-attempt candidates).
pub const MOORE_26: [(i32, i32, i32); 26] = [
    (-1,-1,-1),(-1,-1, 0),(-1,-1, 1),
    (-1, 0,-1),(-1, 0, 0),(-1, 0, 1),
    (-1, 1,-1),(-1, 1, 0),(-1, 1, 1),
    ( 0,-1,-1),( 0,-1, 0),( 0,-1, 1),
    ( 0, 0,-1),            ( 0, 0, 1),
    ( 0, 1,-1),( 0, 1, 0),( 0, 1, 1),
    ( 1,-1,-1),( 1,-1, 0),( 1,-1, 1),
    ( 1, 0,-1),( 1, 0, 0),( 1, 0, 1),
    ( 1, 1,-1),( 1, 1, 0),( 1, 1, 1),
];

/// 6-connected Von Neumann neighbourhood (surface counting).
pub const VN6: [(i32, i32, i32); 6] = [
    (-1, 0, 0),(1, 0, 0),
    ( 0,-1, 0),(0, 1, 0),
    ( 0, 0,-1),(0, 0, 1),
];

// ── Serialisable snapshot ─────────────────────────────────────────────────────

#[derive(Serialize, Deserialize)]
pub struct SaveState {
    pub mcs:    usize,
    pub params: Params,
    /// Flat grid, z-major: index = z*W*H + y*W + x
    pub grid:   Vec<u32>,
    pub cells:  Vec<CellState>,
}

// ── Simulation ────────────────────────────────────────────────────────────────

pub struct Cpm3d {
    pub p:        Params,
    pub mcs:      usize,
    pub grid:     Vec<u32>,
    pub cells:    Vec<CellState>,
    pub mcs_size: usize,
    pub rng:      StdRng,
}

impl Cpm3d {
    // ── Construction ──────────────────────────────────────────────────────────

    pub fn new(p: Params) -> Self {
        let n = p.n_cells;
        let mcs_size = p.mcs_per_step.unwrap_or(p.grid_w * p.grid_h * p.grid_d);

        // cell 0 = medium (alive=false keeps it out of demography loops)
        // If volume_sigma > 0, sample each cell's target_volume independently
        // from N(target_volume, volume_sigma²), clamped to [1, ∞).
        let mut rng_init = StdRng::from_entropy();
        let normal_dist = if p.volume_sigma > 0.0 {
            Normal::new(p.target_volume as f64, p.volume_sigma).ok()
        } else {
            None
        };
        let cells: Vec<CellState> = (0..=n)
            .map(|k| {
                let tv = if k == 0 {
                    p.target_volume
                } else if let Some(ref dist) = normal_dist {
                    let v = dist.sample(&mut rng_init).round() as i64;
                    v.max(1)
                } else {
                    p.target_volume
                };
                // Recompute target_surface to match sampled volume (sphere approximation).
                let ts = ((36.0 * std::f64::consts::PI).powf(1.0 / 3.0)
                    * (tv as f64).powf(2.0 / 3.0)).round() as i64;
                let mut c = CellState::new(k as u32, tv, ts);
                if k == 0 { c.alive = false; }
                c
            })
            .collect();

        let mut sim = Self {
            p: p.clone(),
            mcs: 0,
            grid: vec![0u32; p.grid_w * p.grid_h * p.grid_d],
            cells,
            mcs_size,
            rng: StdRng::from_entropy(),
        };

        fs::create_dir_all(&sim.p.out_dir).expect("cannot create out_dir");

        // Place spherical seeds.  When volume_sigma > 0, each cell gets its own
        // radius derived from its individually-sampled target_volume.
        if p.volume_sigma > 0.0 {
            place_cells_spheres_individual(
                &mut sim.grid,
                p.grid_w, p.grid_h, p.grid_d,
                &sim.cells,
            );
        } else {
            let radius = ((3.0 / (4.0 * std::f64::consts::PI))
                * p.target_volume as f64)
                .powf(1.0 / 3.0);
            place_cells_spheres(
                &mut sim.grid,
                p.grid_w, p.grid_h, p.grid_d,
                p.n_cells, radius,
            );
        }

        sim.recompute_stats();
        sim
    }

    pub fn from_save(s: SaveState) -> Self {
        let mcs_size = s.params.mcs_per_step
            .unwrap_or(s.params.grid_w * s.params.grid_h * s.params.grid_d);
        fs::create_dir_all(&s.params.out_dir).expect("cannot create out_dir");
        Self {
            p: s.params,
            mcs: s.mcs,
            grid: s.grid,
            cells: s.cells,
            mcs_size,
            rng: StdRng::from_entropy(),
        }
    }

    // ── Statistics ────────────────────────────────────────────────────────────

    pub fn recompute_stats(&mut self) {
        let (w, h, d) = (self.p.grid_w, self.p.grid_h, self.p.grid_d);
        for c in self.cells.iter_mut() { c.volume = 0; c.surface = 0; }

        for z in 0..d {
            for y in 0..h {
                for x in 0..w {
                    let s = self.grid[z * w * h + y * w + x] as usize;
                    if s == 0 { continue; }
                    self.cells[s].volume += 1;
                    for (dz, dy, dx) in VN6 {
                        let nx = x as i32 + dx;
                        let ny = y as i32 + dy;
                        let nz = z as i32 + dz;
                        let nb = if nx < 0 || nx >= w as i32
                                || ny < 0 || ny >= h as i32
                                || nz < 0 || nz >= d as i32 {
                            0usize
                        } else {
                            self.grid[nz as usize * w * h + ny as usize * w + nx as usize] as usize
                        };
                        if nb != s { self.cells[s].surface += 1; }
                    }
                }
            }
        }
    }

    // ── Hamiltonian helpers ───────────────────────────────────────────────────

    #[inline]
    fn delta_h_adhesion(&self, x: usize, y: usize, z: usize, s_old: u32, s_new: u32) -> f64 {
        let (w, h, d) = (self.p.grid_w as i32, self.p.grid_h as i32, self.p.grid_d as i32);
        let mut dh = 0.0f64;
        for (dz, dy, dx) in MOORE_26 {
            let nx = x as i32 + dx;
            let ny = y as i32 + dy;
            let nz = z as i32 + dz;
            let nb = if nx < 0 || nx >= w || ny < 0 || ny >= h || nz < 0 || nz >= d {
                0u32
            } else {
                self.grid[nz as usize * self.p.grid_w * self.p.grid_h
                         + ny as usize * self.p.grid_w + nx as usize]
            };
            dh += j(s_new, nb, self.p.j_cell_medium, self.p.j_cell_cell)
                - j(s_old, nb, self.p.j_cell_medium, self.p.j_cell_cell);
        }
        dh
    }

    #[inline]
    fn delta_h_volume(&self, s_old: u32, s_new: u32) -> f64 {
        let lam = self.p.lambda_vol;
        let penalty = self.p.small_volume_penalty;
        let n = self.p.small_volume_n;
        let mut dh = 0.0;
        if s_old > 0 {
            let c = &self.cells[s_old as usize];
            dh += delta_volume_loss(c.volume, c.target_volume, lam);
            dh +=  penalty / (c.volume as f64).pow(n);
        }
        if s_new > 0 {
            let c = &self.cells[s_new as usize];
            dh += delta_volume_gain(c.volume, c.target_volume, lam);
        }
        dh
    }

    /// Computes both the energy delta and returns (ds_old, ds_new) for local update.
    fn delta_h_surface_with_ds(
        &self, x: usize, y: usize, z: usize,
        s_old: u32, s_new: u32,
    ) -> (f64, i64, i64) {
        let (w, h, d) = (self.p.grid_w as i32, self.p.grid_h as i32, self.p.grid_d as i32);
        let lam = self.p.lambda_surf;
        let mut ds_old = 0i64;
        let mut ds_new = 0i64;

        for (dz, dy, dx) in VN6 {
            let nx = x as i32 + dx;
            let ny = y as i32 + dy;
            let nz = z as i32 + dz;
            let nb = if nx < 0 || nx >= w || ny < 0 || ny >= h || nz < 0 || nz >= d {
                0u32
            } else {
                self.grid[nz as usize * self.p.grid_w * self.p.grid_h
                         + ny as usize * self.p.grid_w + nx as usize]
            };
            if s_old > 0 {
                if nb != s_old { ds_old -= 1; } else { ds_old += 1; }
            }
            if s_new > 0 {
                if nb != s_new { ds_new += 1; } else { ds_new -= 1; }
            }
        }

        let mut dh = 0.0;
        if s_old > 0 {
            let c = &self.cells[s_old as usize];
            dh += delta_surface(c.surface, c.target_surface, ds_old, lam);
        }
        if s_new > 0 {
            let c = &self.cells[s_new as usize];
            dh += delta_surface(c.surface, c.target_surface, ds_new, lam);
        }
        (dh, ds_old, ds_new)
    }

    fn delta_h_sphericity(&self, s_old: u32, s_new: u32, ds_old: i64, ds_new: i64) -> f64 {
        let lam = self.p.lambda_spher;
        let mut dh = 0.0;
        if s_old > 0 {
            let c = &self.cells[s_old as usize];
            dh += delta_sphericity(c.volume, c.surface,
                                   c.volume - 1, c.surface + ds_old, lam);
        }
        if s_new > 0 {
            let c = &self.cells[s_new as usize];
            dh += delta_sphericity(c.volume, c.surface,
                                   c.volume + 1, c.surface + ds_new, lam);
        }
        dh
    }

    // ── Monte Carlo step ──────────────────────────────────────────────────────

    fn attempt(&mut self) {
        let (w, h, d) = (self.p.grid_w, self.p.grid_h, self.p.grid_d);

        let x = self.rng.gen_range(0..w);
        let y = self.rng.gen_range(0..h);
        let z = self.rng.gen_range(0..d);
        let s_old = self.grid[z * w * h + y * w + x];

        // Collect Moore-26 neighbours with a different sigma
        let mut candidates: Vec<u32> = Vec::with_capacity(26);
        for (dz, dy, dx) in MOORE_26 {
            let nx = x as i32 + dx;
            let ny = y as i32 + dy;
            let nz = z as i32 + dz;
            if nx < 0 || nx >= w as i32 || ny < 0 || ny >= h as i32 || nz < 0 || nz >= d as i32 {
                continue;
            }
            let nb = self.grid[nz as usize * w * h + ny as usize * w + nx as usize];
            if nb != s_old { candidates.push(nb); }
        }
        if candidates.is_empty() { return; }

        let s_new = candidates[self.rng.gen_range(0..candidates.len())];

        let dh_adh  = self.delta_h_adhesion(x, y, z, s_old, s_new);
        let dh_vol  = self.delta_h_volume(s_old, s_new);
        let (dh_surf, ds_old, ds_new) = self.delta_h_surface_with_ds(x, y, z, s_old, s_new);
        let dh_sph  = self.delta_h_sphericity(s_old, s_new, ds_old, ds_new);
        
        let dh = dh_adh + dh_vol + dh_surf + dh_sph;
        let accept = dh <= 0.0
            || self.rng.gen_range(0.0f64..1.0) < (-dh / self.p.temperature).exp();

        if accept {
            self.grid[z * w * h + y * w + x] = s_new;
            if s_old > 0 {
                self.cells[s_old as usize].volume  -= 1;
                self.cells[s_old as usize].surface += ds_old;
            }
            if s_new > 0 {
                self.cells[s_new as usize].volume  += 1;
                self.cells[s_new as usize].surface += ds_new;
            }
        }
    }

    pub fn run_mcs(&mut self) {
        for _ in 0..self.mcs_size { self.attempt(); }
        self.mcs += 1;
    }

    // ── Console summary ───────────────────────────────────────────────────────

    pub fn print_stats(&self) {
        let living: Vec<_> = self.cells.iter().filter(|c| c.id > 0 && c.alive).collect();
        let mean_vol = if living.is_empty() { 0.0 } else {
            living.iter().map(|c| c.volume as f64).sum::<f64>() / living.len() as f64
        };
        println!("MCS {:5}  cells={:3}  mean_vol={:.1}", self.mcs, living.len(), mean_vol);
    }

    // ── Slice PNG export ──────────────────────────────────────────────────────

    pub fn save_slice_png(&self, axis: u8, slice_idx: usize, path: &str) {
        let (w, h, d) = (self.p.grid_w, self.p.grid_h, self.p.grid_d);
        let n_cells = self.cells.len() as u32;

        let (pw, ph, pixels): (usize, usize, Vec<u8>) = match axis {
            0 => { // XY slice at z = slice_idx
                let z = slice_idx.min(d - 1);
                let px = (0..h).flat_map(|y| (0..w).flat_map(move |x| {
                    sigma_rgb(self.grid[z * w * h + y * w + x], n_cells)
                })).collect();
                (w, h, px)
            }
            1 => { // XZ slice at y = slice_idx
                let y = slice_idx.min(h - 1);
                let px = (0..d).flat_map(|z| (0..w).flat_map(move |x| {
                    sigma_rgb(self.grid[z * w * h + y * w + x], n_cells)
                })).collect();
                (w, d, px)
            }
            _ => { // YZ slice at x = slice_idx
                let x = slice_idx.min(w - 1);
                let px = (0..d).flat_map(|z| (0..h).flat_map(move |y| {
                    sigma_rgb(self.grid[z * w * h + y * w + x], n_cells)
                })).collect();
                (h, d, px)
            }
        };

        write_png(path, pw, ph, &pixels);
    }

    // ── Serialisation ─────────────────────────────────────────────────────────

    pub fn save_state(&self, path: Option<&str>) {
        let default = format!("{}/state_mcs{:06}.json", self.p.out_dir.trim_end_matches('/'), self.mcs);
        let path = path.unwrap_or(&default);
        let s = SaveState {
            mcs: self.mcs,
            params: self.p.clone(),
            grid: self.grid.clone(),
            cells: self.cells.clone(),
        };
        let json = serde_json::to_string_pretty(&s).expect("serialisation failed");
        fs::write(path, json).unwrap_or_else(|e| eprintln!("save error: {e}"));
        println!("  [json] → {path}");
    }

    pub fn load_state(path: &str) -> Self {
        let json = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
        let s: SaveState = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("cannot parse {path}: {e}"));
        println!("  [json] ← {path}  (MCS={})", s.mcs);
        Self::from_save(s)
    }
}

// ── Centroid utilities ────────────────────────────────────────────────────────

/// World-space centroid for each cell sigma (index = sigma).
///
/// Returns a `Vec` of length `n_cells`.  Index 0 (medium) is always `[0,0,0]`.
/// Entries for empty or dead cells are also `[0,0,0]`; callers should gate on
/// `cells[sigma].alive && cells[sigma].volume > 0`.
pub fn cell_centroids(
    grid:    &[u32],
    w: usize, h: usize, d: usize,
    n_cells: usize,
) -> Vec<[f32; 3]> {
    let mut sx  = vec![0f64; n_cells];
    let mut sy  = vec![0f64; n_cells];
    let mut sz  = vec![0f64; n_cells];
    let mut cnt = vec![0u64; n_cells];

    for z in 0..d {
        for y in 0..h {
            for x in 0..w {
                let s = grid[z * w * h + y * w + x] as usize;
                if s == 0 || s >= n_cells { continue; }
                sx[s]  += x as f64;
                sy[s]  += y as f64;
                sz[s]  += z as f64;
                cnt[s] += 1;
            }
        }
    }

    (0..n_cells)
        .map(|s| {
            if cnt[s] == 0 {
                [0.0f32; 3]
            } else {
                [(sx[s] / cnt[s] as f64) as f32,
                 (sy[s] / cnt[s] as f64) as f32,
                 (sz[s] / cnt[s] as f64) as f32]
            }
        })
        .collect()
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn sigma_rgb(sigma: u32, n: u32) -> [u8; 3] {
    if sigma == 0 { return [20, 20, 30]; }
    let h = (sigma as f64 / n as f64 * 360.0) as u32 % 360;
    let (r, g, b) = hsv_to_rgb_u8(h as f64, 0.75, 0.9);
    [r, g, b]
}

fn hsv_to_rgb_u8(h: f64, s: f64, v: f64) -> (u8, u8, u8) {
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;
    let (r, g, b) = match h as u32 / 60 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    (((r + m) * 255.0) as u8, ((g + m) * 255.0) as u8, ((b + m) * 255.0) as u8)
}

fn write_png(path: &str, w: usize, h: usize, rgb: &[u8]) {
    let f = std::io::BufWriter::new(fs::File::create(path).expect("png create failed"));
    let mut enc = png::Encoder::new(f, w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header().expect("png header failed")
       .write_image_data(rgb).expect("png write failed");
}
