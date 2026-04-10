use std::fs;
use rand::prelude::StdRng;
use serde::{Deserialize, Serialize};
use crate::boundary::{build_boundary, pixel_is_boundary};
use crate::energy::delta_isoperimetric;
use crate::params::Params;
use crate::cellstate::CellState;
use rand::prelude::*;
use crate::render::arrayvec_neighbours;
pub const MOORE: [(i32, i32); 8] = [
    (-1, -1), (-1, 0), (-1, 1),
    ( 0, -1),           ( 0, 1),
    ( 1, -1),  ( 1, 0), ( 1, 1),
];

/// Von Neumann (4-connected), used for perimeter counting
pub const VON_NEUMANN: [(i32, i32); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];

// ─── Serialisable state ───────────────────────────────────────────────────────

#[derive(Serialize, Deserialize)]
pub struct SaveState {
    pub mcs: usize,
    pub params: Params,
    /// Flat grid, row-major: index = r * grid_w + c
    pub grid: Vec<u32>,
    pub cells:  Vec<CellState>,
}

// ─── Simulation ───────────────────────────────────────────────────────────────

pub struct Cpm2d {
    pub p: Params,
    pub mcs: usize,
    /// sigma values: 0 = medium, 1..=n_cells = cells
    pub grid: Vec<u32>,
    /// area[0] unused; area[k] = pixel count of cell k
    //pub area: Vec<i64>,
    /// perim[k] = 4-connected perimeter of cell k
    //pub perimeter: Vec<i64>,
    // new version : target area/perimeter is cell based
    pub cells:   Vec<CellState>,
    /// boundary[w*h] with id's k of boundary cells
    pub boundary: Vec<u32>,
    mcs_size: usize,
    rng: StdRng,
}

impl Cpm2d {
    // ── Construction ─────────────────────────────────────────────────────────

    pub fn new(p: Params) -> Self {
        let n = p.n_cells;
        let mcs_size = p.mcs_per_step.unwrap_or(p.grid_w * p.grid_h);
        let cells = (0..n)
            .map(|k| CellState::new((k+1) as u32,p.target_area, p.target_perim))
            .collect();
        let mut sim = Self {
            p: p.clone(),
            mcs: 0,
            grid: vec![0u32; p.grid_w * p.grid_h],
            cells : cells,
            boundary: vec![0u32; n + 1],
            mcs_size: mcs_size,
            rng: StdRng::from_entropy()
        };
        fs::create_dir_all(&sim.p.png_dir).expect("cannot create png_dir");
        sim.place_initial_cells();
        sim.recompute_stats();
        sim.boundary = build_boundary(&sim.grid, p.grid_w, p.grid_h);
        sim
    }

    pub fn from_save(s: SaveState) -> Self {
        let n = s.params.n_cells;
        let mcs_size = s.params.mcs_per_step.unwrap_or(s.params.grid_w * s.params.grid_h);
        fs::create_dir_all(&s.params.png_dir).expect("cannot create png_dir");
        let mut boundary = build_boundary(&s.grid, s.params.grid_w, s.params.grid_h);
        Self {
            p: s.params.clone(),
            mcs: s.mcs,
            grid: s.grid,
            cells : s.cells,
            boundary: boundary,
            mcs_size: mcs_size,
            rng : StdRng::from_entropy()
        }
    }

    // ── Initial conditions ────────────────────────────────────────────────────

    fn place_initial_cells(&mut self) {
        let n_total = self.p.n_cells;
        let n_cols = (n_total as f64).sqrt().ceil() as usize;
        let n_rows = (n_total + n_cols - 1) / n_cols;
        let W = self.p.grid_w;
        let H = self.p.grid_h;

        let tile_w = W / n_cols;
        let tile_h = H / n_rows;
        let ta = (self.p.target_area / 2) as f64;
        let blob_side = (ta.sqrt()
            .round() as usize)
            .max(2)
            .min(tile_w.saturating_sub(2))
            .min(tile_h.saturating_sub(2));

        for k in 0..n_total {
            let sigma = (k + 1) as u32;
            let col = k % n_cols;
            let row = k / n_cols;

            let tx = col * tile_w;
            let ty = row * tile_h;

            let bx = tx + (tile_w - blob_side) / 2;
            let by = ty + (tile_h - blob_side) / 2;

            for dy in 0..blob_side {
                for dx in 0..blob_side {
                    let r = by + dy;
                    let c = bx + dx;
                    if r < H && c < W {
                        self.grid[r * W + c] = sigma;
                    }
                }
            }
        }
    }

    // ── Statistics ────────────────────────────────────────────────────────────

    fn recompute_stats(&mut self) {
        let W = self.p.grid_w;
        let H = self.p.grid_h;
        self.cells.iter_mut().for_each(|c| {c.area=0; c.perimeter=0;});

        for r in 0..H {
            for c in 0..W {
                let s = self.grid[r * W + c] as usize;
                if s > 0 {
                    self.cells[s].area += 1;
                    for (dr, dc) in VON_NEUMANN {
                        let nr = r as i32 + dr;
                        let nc = c as i32 + dc;
                        let nb = if nr < 0 || nr >= H as i32 || nc < 0 || nc >= W as i32 {
                            0usize
                        } else {
                            self.grid[nr as usize * W + nc as usize] as usize
                        };
                        if nb != s {
                            self.cells[s].perimeter += 1;
                        }
                    }
                }
            }
        }
    }
    /// Recompute boundary status for (r,c) and all its Moore neighbours.
    /// Call AFTER self.grid[r*W+c] has been updated to s_new.
    fn update_boundary_local(&mut self, r: usize, c: usize) {
        let W = self.p.grid_w;
        let H = self.p.grid_h;

        // The flipped pixel itself + its 8 Moore neighbours = at most 9 pixels
        let mut to_recheck: [Option<(usize, usize)>; 9] = [None; 9];
        to_recheck[0] = Some((r, c));
        for (i, (dr, dc)) in MOORE.iter().enumerate() {
            let nr = r as i32 + dr;
            let nc = c as i32 + dc;
            if nr >= 0 && nr < H as i32 && nc >= 0 && nc < W as i32 {
                to_recheck[i + 1] = Some((nr as usize, nc as usize));
            }
        }

        for cell in to_recheck.into_iter().flatten() {
            let (rr, cc) = cell;
            let s = self.grid[rr * W + cc];
            let idx = rr * W + cc;
            self.boundary[idx] = if s > 0
                && pixel_is_boundary(&self.grid, W, H, rr, cc, s)
            {
                s
            } else {
                0
            };
        }
    }
    // ── Hamiltonian helpers ───────────────────────────────────────────────────

    #[inline]
    fn j(&self, s1: u32, s2: u32) -> f64 {
        if s1 == s2 { return 0.0; }
        if s1 == 0 || s2 == 0 { self.p.j_cell_medium } else { self.p.j_cell_cell }
    }

    fn delta_h_contact(&self, r: usize, c: usize, s_old: u32, s_new: u32) -> f64 {
        let W = self.p.grid_w as i32;
        let H = self.p.grid_h as i32;
        let mut dh = 0.0f64;
        for (dr, dc) in MOORE {
            let nr = r as i32 + dr;
            let nc = c as i32 + dc;
            let nb = if nr < 0 || nr >= H || nc < 0 || nc >= W {
                0u32
            } else {
                self.grid[nr as usize * self.p.grid_w + nc as usize]
            };
            dh += self.j(s_new, nb) - self.j(s_old, nb);
        }
        dh
    }

    fn delta_h_area(&self, s_old: u32, s_new: u32) -> f64 {
        let lam = self.p.lambda_area;
        let at = self.p.target_area;
        let mut dh = 0.0f64;
        if s_old > 0 {
            let a = self.cells[s_old as usize].area;
            dh += lam * ((a - 1 - at).pow(2) - (a - at).pow(2)) as f64;
        }
        if s_new > 0 {
            let a = self.cells[s_new as usize].area;
            dh += lam * ((a + 1 - at).pow(2) - (a - at).pow(2)) as f64;
        }
        dh
    }

    fn delta_h_perim(&self, r: usize, c: usize, s_old: u32, s_new: u32) -> f64 {
        let W = self.p.grid_w as i32;
        let H = self.p.grid_h as i32;
        let lam = self.p.lambda_perim;
        let pt = self.p.target_perim;

        // dp[s] = change in 4-connected perimeter of cell s if (r,c) flips s_old→s_new
        let mut dp_old = 0i64;
        let mut dp_new = 0i64;

        for (dr, dc) in VON_NEUMANN {
            let nr = r as i32 + dr;
            let nc = c as i32 + dc;
            let nb = if nr < 0 || nr >= H || nc < 0 || nc >= W {
                0u32
            } else {
                self.grid[nr as usize * self.p.grid_w + nc as usize]
            };

            // Effect on s_old: it loses pixel (r,c)
            if s_old > 0 {
                // edge (r,c)-(nr,nc) was a perim edge for s_old if nb != s_old
                if nb != s_old { dp_old -= 1; }
                // after flip, neighbour (nr,nc) now borders s_new; if nb==s_old it gains a perim edge
                if nb == s_old { dp_old += 1; }
            }

            // Effect on s_new: it gains pixel (r,c)
            if s_new > 0 {
                if nb != s_new { dp_new += 1; }
                if nb == s_new { dp_new -= 1; }
            }
        }

        let mut dh = 0.0f64;
        if s_old > 0 && dp_old != 0 {
            let p = self.cells[s_old as usize].perimeter;
            dh += lam * ((p + dp_old - pt).pow(2) - (p - pt).pow(2)) as f64;
        }
        if s_new > 0 && dp_new != 0 {
            let p = self.cells[s_new as usize].perimeter;
            dh += lam * ((p + dp_new - pt).pow(2) - (p - pt).pow(2)) as f64;
        }
        dh
    }

    // ── Monte Carlo step ──────────────────────────────────────────────────────

    fn attempt(&mut self) {
        let W = self.p.grid_w;
        let H = self.p.grid_h;
        let T = self.p.temperature;

        let r = self.rng.gen_range(0..H);
        let c = self.rng.gen_range(0..W);
        let s_old = self.grid[r * W + c];
        let old_area  = self.cells[s_old as usize].area as f64;
        let old_perim  = self.cells[s_old as usize].perimeter as f64;

        // Collect in-bounds neighbours with a different sigma
        let mut candidates = arrayvec_neighbours(r, c, W, H, &self.grid, s_old);
        if candidates.is_empty() { return; }

        // Pick one at random
        let idx = self.rng.gen_range(0..candidates.len());
        let s_new = candidates[idx];

        let new_area  = self.cells[s_new as usize].area as f64;
        let new_perim  = self.cells[s_new as usize].perimeter as f64;

        let lambda_iso = self.p.lambda_iso;
        // ΔH
        let dh = self.delta_h_contact(r, c, s_old, s_new)
            + self.delta_h_area(s_old, s_new)
            + self.delta_h_perim(r, c, s_old, s_new)
            + delta_isoperimetric(old_area, old_perim,new_area, new_perim,lambda_iso);

        // Metropolis
        let accept = dh <= 0.0 || self.rng.gen_range(0f64..1f64) < (-dh / T).exp();
        if accept {
            self.grid[r * W + c] = s_new;
            if s_old > 0 { self.cells[s_old as usize].area -= 1; }
            if s_new > 0 { self.cells[s_new as usize].area += 1; }
            self.update_perimeter_local(r, c, s_old, s_new);
            self.update_boundary_local(r, c);
        }
    }

    pub fn run_mcs(&mut self) {
        for _ in 0..self.mcs_size {
            self.attempt();
        }
        self.mcs += 1;
    }

    // ── Local perimeter update ────────────────────────────────────────────────

    fn update_perimeter_local(&mut self, r: usize, c: usize, s_old: u32, s_new: u32) {
        let W = self.p.grid_w as i32;
        let H = self.p.grid_h as i32;

        // Pixels whose perimeter count might have changed: (r,c) + 4-neighbours
        let mut affected: Vec<(usize, usize)> = vec![(r, c)];
        for (dr, dc) in VON_NEUMANN {
            let nr = r as i32 + dr;
            let nc = c as i32 + dc;
            if nr >= 0 && nr < H && nc >= 0 && nc < W {
                affected.push((nr as usize, nc as usize));
            }
        }

        for (rr, cc) in affected {
            let s = self.grid[rr * self.p.grid_w + cc] as usize;
            if s == 0 { continue; }

            // The "old" sigma at (r,c) for the purposes of counting edges
            // For the flipped pixel itself: old was s_old; for neighbours: unchanged.
            let old_at_flip = if (rr, cc) == (r, c) { s_old } else { s as u32 };

            let mut old_p = 0i64;
            let mut new_p = 0i64;
            for (dr, dc) in VON_NEUMANN {
                let nr2 = rr as i32 + dr;
                let nc2 = cc as i32 + dc;
                let nb_old = if nr2 < 0 || nr2 >= H || nc2 < 0 || nc2 >= W {
                    0u32
                } else if (nr2 as usize, nc2 as usize) == (r, c) {
                    s_old   // what was there before the flip
                } else {
                    self.grid[nr2 as usize * self.p.grid_w + nc2 as usize]
                };
                let nb_new = if nr2 < 0 || nr2 >= H || nc2 < 0 || nc2 >= W {
                    0u32
                } else {
                    self.grid[nr2 as usize * self.p.grid_w + nc2 as usize]
                };
                if nb_old != old_at_flip { old_p += 1; }
                if nb_new != s as u32    { new_p += 1; }
            }

            let diff = new_p - old_p;
            if diff != 0 {
                self.cells[s].perimeter += diff;
            }
        }
    }

    // ── Serialisation ─────────────────────────────────────────────────────────

    pub fn save_state(&self, path: Option<&str>) {
        let default = format!("{}/state_mcs{:06}.json", self.p.frames_dir, self.mcs);
        let path = path.unwrap_or(&default);
        let s = SaveState {
            mcs: self.mcs,
            params: self.p.clone(),
            grid: self.grid.clone(),
            cells: self.cells.clone(),
        };
        let json = serde_json::to_string_pretty(&s).expect("serialisation failed");
        fs::write(path, json).unwrap_or_else(|e| eprintln!("JSON save error: {e}"));
        println!("  [json] state saved → {path}");
    }

    pub fn load_state(path: &str) -> Self {
        let json = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("Cannot read {path}: {e}"));
        let s: SaveState = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("Cannot parse {path}: {e}"));
        println!("  [json] state loaded ← {path}  (MCS={})", s.mcs);
        Self::from_save(s)
    }
}

