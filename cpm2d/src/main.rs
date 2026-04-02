//! 2-D Cellular Potts Model (CPM)
//!
//! Features
//! --------
//! - Hard (non-periodic) boundary conditions
//! - Square-block initial conditions, cells grow toward target area
//! - PNG frame output via the `image` crate
//! - ANSI coloured console display
//! - JSON serialisation / deserialisation (serde_json)
//!
//! Usage
//! -----
//!   cargo run --release                        # fresh run
//!   cargo run --release -- --load state.json  # resume from save

use std::f64;
use std::fs;

use image::{ImageBuffer, Rgb};
use rand::prelude::*;
use serde::{Deserialize, Serialize};

// ─── Parameters ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Params {
    // Grid dimensions
    pub grid_w: usize,
    pub grid_h: usize,

    // Cells
    pub n_cells: usize,
    pub target_area: i64,
    pub target_perim: i64,

    // Hamiltonian weights
    pub lambda_area: f64,
    pub lambda_perim: f64,

    // Surface tensions
    pub j_cell_medium: f64,
    pub j_cell_cell: f64,

    // Monte Carlo
    pub temperature: f64,
    /// None → grid_w * grid_h attempts per MCS
    pub mcs_per_step: Option<usize>,
    pub total_steps: usize,

    // Output
    pub png_dir: String,
    pub png_every: usize,
    pub console_every: usize,
    pub save_every: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            grid_w: 60,
            grid_h: 60,
            n_cells: 9,
            target_area: 300,
            target_perim: 300f64.sqrt() as i64,
            lambda_area: 1.0,
            lambda_perim: 0.1,
            j_cell_medium: 10.0,
            j_cell_cell: 2.0,
            temperature: 20.0,
            mcs_per_step: None,
            total_steps: 1000,
            png_dir: "frames".to_string(),
            png_every: 10,
            console_every: 10,
            save_every: 50,
        }
    }
}

// ─── Moore neighbourhood offsets (8-connected) ───────────────────────────────

const MOORE: [(i32, i32); 8] = [
    (-1, -1), (-1, 0), (-1, 1),
    ( 0, -1),           ( 0, 1),
    ( 1, -1),  ( 1, 0), ( 1, 1),
];

/// Von Neumann (4-connected), used for perimeter counting
const VON_NEUMANN: [(i32, i32); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];

// ─── Serialisable state ───────────────────────────────────────────────────────

#[derive(Serialize, Deserialize)]
pub struct SaveState {
    pub mcs: usize,
    pub params: Params,
    /// Flat grid, row-major: index = r * grid_w + c
    pub grid: Vec<u32>,
    pub area: Vec<i64>,
    pub perim: Vec<i64>,
}

// ─── Simulation ───────────────────────────────────────────────────────────────

pub struct Cpm2d {
    pub p: Params,
    pub mcs: usize,
    /// sigma values: 0 = medium, 1..=n_cells = cells
    pub grid: Vec<u32>,
    /// area[0] unused; area[k] = pixel count of cell k
    pub area: Vec<i64>,
    /// perim[k] = 4-connected perimeter of cell k
    pub perim: Vec<i64>,
    mcs_size: usize,
    rng: StdRng,
}

impl Cpm2d {
    // ── Construction ─────────────────────────────────────────────────────────

    pub fn new(p: Params) -> Self {
        let n = p.n_cells;
        let mcs_size = p.mcs_per_step.unwrap_or(p.grid_w * p.grid_h);
        let mut sim = Self {
            p: p.clone(),
            mcs: 0,
            grid: vec![0u32; p.grid_w * p.grid_h],
            area: vec![0i64; n + 1],
            perim: vec![0i64; n + 1],
            mcs_size: mcs_size,
            rng: StdRng::from_entropy()
        };
        fs::create_dir_all(&sim.p.png_dir).expect("cannot create png_dir");
        sim.place_initial_cells();
        sim.recompute_stats();
        sim
    }

    pub fn from_save(s: SaveState) -> Self {
        let n = s.params.n_cells;
        let mcs_size = s.params.mcs_per_step.unwrap_or(s.params.grid_w * s.params.grid_h);
        fs::create_dir_all(&s.params.png_dir).expect("cannot create png_dir");
        Self {
            p: s.params.clone(),
            mcs: s.mcs,
            grid: s.grid,
            area: s.area,
            perim: s.perim,
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
        self.area.iter_mut().for_each(|x| *x = 0);
        self.perim.iter_mut().for_each(|x| *x = 0);

        for r in 0..H {
            for c in 0..W {
                let s = self.grid[r * W + c] as usize;
                if s > 0 {
                    self.area[s] += 1;
                    for (dr, dc) in VON_NEUMANN {
                        let nr = r as i32 + dr;
                        let nc = c as i32 + dc;
                        let nb = if nr < 0 || nr >= H as i32 || nc < 0 || nc >= W as i32 {
                            0usize
                        } else {
                            self.grid[nr as usize * W + nc as usize] as usize
                        };
                        if nb != s {
                            self.perim[s] += 1;
                        }
                    }
                }
            }
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
            let a = self.area[s_old as usize];
            dh += lam * ((a - 1 - at).pow(2) - (a - at).pow(2)) as f64;
        }
        if s_new > 0 {
            let a = self.area[s_new as usize];
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
            let p = self.perim[s_old as usize];
            dh += lam * ((p + dp_old - pt).pow(2) - (p - pt).pow(2)) as f64;
        }
        if s_new > 0 && dp_new != 0 {
            let p = self.perim[s_new as usize];
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

        // Collect in-bounds neighbours with a different sigma
        let mut candidates = arrayvec_neighbours(r, c, W, H, &self.grid, s_old);
        if candidates.is_empty() { return; }

        // Pick one at random
        let idx = self.rng.gen_range(0..candidates.len());
        let s_new = candidates[idx];

        // ΔH
        let dh = self.delta_h_contact(r, c, s_old, s_new)
            + self.delta_h_area(s_old, s_new)
            + self.delta_h_perim(r, c, s_old, s_new);

        // Metropolis
        let accept = dh <= 0.0 || self.rng.gen_range(0f64..1f64) < (-dh / T).exp();
        if accept {
            self.grid[r * W + c] = s_new;
            if s_old > 0 { self.area[s_old as usize] -= 1; }
            if s_new > 0 { self.area[s_new as usize] += 1; }
            self.update_perim_local(r, c, s_old, s_new);
        }
    }

    pub fn run_mcs(&mut self) {
        for _ in 0..self.mcs_size {
            self.attempt();
        }
        self.mcs += 1;
    }

    // ── Local perimeter update ────────────────────────────────────────────────

    fn update_perim_local(&mut self, r: usize, c: usize, s_old: u32, s_new: u32) {
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
                self.perim[s] += diff;
            }
        }
    }

    // ── Console output ────────────────────────────────────────────────────────

    pub fn console_print(&self) {
        // ANSI background colours: 0=reset, 1=red, 2=green, 3=yellow, 4=blue, …
        const BG: &[&str] = &[
            "\x1b[0m  ",   // 0 medium
            "\x1b[41m  ",  // 1 red
            "\x1b[42m  ",  // 2 green
            "\x1b[43m  ",  // 3 yellow
            "\x1b[44m  ",  // 4 blue
            "\x1b[45m  ",  // 5 magenta
            "\x1b[46m  ",  // 6 cyan
            "\x1b[47m  ",  // 7 white
            "\x1b[101m  ", // 8 bright red
        ];
        let reset = "\x1b[0m";
        let W = self.p.grid_w;
        let H = self.p.grid_h;
        let n = self.p.n_cells;

        println!("\n─── MCS {} ───", self.mcs);
        for r in 0..H {
            for c in 0..W {
                let s = self.grid[r * W + c] as usize;
                let [red, grn, blu] = cell_colour(s, n);
                print!("\x1b[48;2;{red};{grn};{blu}m  ");
               // let idx = s.min(BG.len() - 1);
                //print!("{}", BG[idx]);
            }
            println!("{}", reset);
        }
        for k in 1..=self.p.n_cells {
            println!(
                "  cell {}: area={:4} (tgt {:3})  perim={:4} (tgt {:3})",
                k, self.area[k], self.p.target_area,
                self.perim[k], self.p.target_perim
            );
        }
    }

    // ── PNG output ────────────────────────────────────────────────────────────

    pub fn save_png(&self, path: Option<&str>) {
        let default = format!("{}/frame_{:06}.png", self.p.png_dir, self.mcs);
        let path = path.unwrap_or(&default);

        const COLOURS: &[[u8; 3]] = &[
            [230, 230, 230], // 0 medium
            [220,  60,  60], // 1 red
            [ 60, 180,  60], // 2 green
            [240, 200,  40], // 3 yellow
            [ 60, 100, 220], // 4 blue
            [180,  60, 200], // 5 magenta
            [ 60, 200, 210], // 6 cyan
            [200, 140,  60], // 7 orange
        ];

        let scale = 8u32;
        let W = self.p.grid_w as u32;
        let H = self.p.grid_h as u32;
        let img_w = W * scale;
        let img_h = H * scale + 20;

        let mut img = ImageBuffer::<Rgb<u8>, _>::new(img_w, img_h);

        // Fill background (label bar)
        for px in img.pixels_mut() {
            *px = Rgb([30u8, 30, 30]);
        }
        let n = self.p.n_cells;
        for r in 0..H {
            for c in 0..W {
                let s = self.grid[(r * W + c) as usize] as usize;
               // let col = COLOURS[s.min(COLOURS.len() - 1)];
                //let rgb = Rgb(col);
                let rgb = Rgb(cell_colour(s, n));
                for dy in 0..scale {
                    for dx in 0..scale {
                        img.put_pixel(c * scale + dx, r * scale + dy, rgb);
                    }
                }
            }
        }

        img.save(path).unwrap_or_else(|e| eprintln!("PNG save error: {e}"));
        println!("  [png] saved → {path}");
    }

    // ── Serialisation ─────────────────────────────────────────────────────────

    pub fn save_state(&self, path: Option<&str>) {
        let default = format!("state_mcs{:06}.json", self.mcs);
        let path = path.unwrap_or(&default);
        let s = SaveState {
            mcs: self.mcs,
            params: self.p.clone(),
            grid: self.grid.clone(),
            area: self.area.clone(),
            perim: self.perim.clone(),
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

// ─── Helper: collect candidate copy-source sigmas ────────────────────────────

fn arrayvec_neighbours(
    r: usize, c: usize,
    W: usize, H: usize,
    grid: &[u32],
    s_old: u32,
) -> Vec<u32> {
    let mut out = Vec::with_capacity(8);
    for (dr, dc) in MOORE {
        let nr = r as i32 + dr;
        let nc = c as i32 + dc;
        if nr >= 0 && nr < H as i32 && nc >= 0 && nc < W as i32 {
            let nb = grid[nr as usize * W + nc as usize];
            if nb != s_old {
                out.push(nb);
            }
        }
    }
    out
}
// ─── Procedural colormap ─────────────────────────────────────────────────────

/// Convert HSV (h in [0,360), s and v in [0,1]) to RGB u8 triple.
fn hsv_to_rgb(h: f32, s: f32, v: f32) -> [u8; 3] {
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;
    let (r1, g1, b1) = match h as u32 {
        0..=59   => (c, x, 0.0),
        60..=119 => (x, c, 0.0),
        120..=179 => (0.0, c, x),
        180..=239 => (0.0, x, c),
        240..=299 => (x, 0.0, c),
        _        => (c, 0.0, x),
    };
    [
        ((r1 + m) * 255.0).round() as u8,
        ((g1 + m) * 255.0).round() as u8,
        ((b1 + m) * 255.0).round() as u8,
    ]
}

/// Return the RGB colour for sigma `s` out of `n_cells` total cells.
/// sigma 0 (medium) is a neutral dark grey.
/// Cells are spread evenly around the hue wheel, with alternating
/// saturation (1.0 / 0.65) and value (0.9 / 1.0) so adjacent indices
/// are visually distinct even when n is large.
fn cell_colour(s: usize, n_cells: usize) -> [u8; 3] {
    if s == 0 {
        return [45, 45, 45]; // medium: dark grey
    }
    let k = s - 1; // 0-based cell index
    let n = n_cells.max(1);
    let hue = (k as f32 / n as f32) * 360.0;
    // Alternate saturation and brightness so neighbouring hues stay distinct
    let sat = if k % 2 == 0 { 0.85 } else { 0.55 };
    let val = if k % 2 == 0 { 1.00 } else { 0.90 };
    hsv_to_rgb(hue, sat, val)
}

// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut sim = if let Some(pos) = args.iter().position(|a| a == "--load") {
        let path = args.get(pos + 1).expect("--load requires a path argument");
        Cpm2d::load_state(path)
    } else {
        Cpm2d::new(Params::default())
    };

    let p = sim.p.clone();
    let end_mcs = sim.mcs + p.total_steps;

    println!(
        "Starting CPM2D  grid={}×{}  cells={}  T={}  MCS={}",
        p.grid_w, p.grid_h, p.n_cells, p.temperature, p.total_steps
    );

    sim.console_print();
    sim.save_png(None);

    while sim.mcs < end_mcs {
        sim.run_mcs();

        if sim.mcs % p.console_every == 0 {
            sim.console_print();
        }
        if sim.mcs % p.png_every == 0 {
            sim.save_png(None);
        }
        if sim.mcs % p.save_every == 0 {
            sim.save_state(None);
        }
    }

    println!("\nDone. Final MCS={}", sim.mcs);
    sim.save_state(Some("state_final.json"));
    sim.save_png(Some(&format!("{}/final.png", p.png_dir)));
}