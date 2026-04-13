//! Cell demography: growth, death, and birth.
//!
//! Call [`Cpm2d::step_demography`] once per MCS (after `run_mcs`) to evolve
//! the cell population:
//!   - Living cells grow linearly (target_area += growth_rate each MCS).
//!   - Each living cell dies with probability `death_prob` per MCS; death is
//!     implemented by setting target_area = 0 so the CPM gradually erodes the
//!     cell.  Once its actual area drops below `death_area_threshold` the slot
//!     is freed.
//!   - With probability `birth_prob` a new 2×2 cell is placed in a free
//!     medium region; dead slots are reused before the cells vec is extended.

use std::f64::consts::PI;
use rand::Rng;
use crate::boundary::build_boundary;
use crate::cellstate::CellState;
use crate::grid::Cpm2d;

impl Cpm2d {
    // ── Public entry point ────────────────────────────────────────────────────

    /// Run one demography step.  Call after `run_mcs`.
    pub fn step_demography(&mut self) {
        self.grow_cells();
        self.maybe_kill_cells();
        self.sweep_dead_cells();
        self.maybe_birth_cell();
    }

    // ── Growth ────────────────────────────────────────────────────────────────

    /// Increase target_area (and matching target_perimeter) for every live,
    /// non-dying cell by `growth_rate` pixels.
    fn grow_cells(&mut self) {
        let rate = self.p.growth_rate as i64;
        if rate == 0 { return; }
        for cell in self.cells.iter_mut() {
            if cell.id == 0 || !cell.alive || cell.dying { continue; }
            cell.target_area += rate;
            cell.target_perimeter = (cell.target_area as f64 * 4.0 * PI).sqrt() as i64;
        }
    }

    // ── Death ─────────────────────────────────────────────────────────────────

    /// Mark each live cell as dying with probability `death_prob`.
    /// A dying cell has target_area = 0 so CPM pressure shrinks it.
    fn maybe_kill_cells(&mut self) {
        let prob = self.p.death_prob;
        if prob <= 0.0 { return; }
        // collect indices to kill so we don't borrow rng while iterating cells
        let to_kill: Vec<usize> = (1..self.cells.len())
            .filter(|&i| {
                let c = &self.cells[i];
                c.alive && !c.dying && self.rng.gen_bool(prob.min(1.0))
            })
            .collect();
        for i in to_kill {
            self.cells[i].dying = true;
            self.cells[i].target_area = 0;
            self.cells[i].target_perimeter = 0;
        }
    }

    /// Remove dying cells whose actual area has dropped below
    /// `death_area_threshold`: zero all their grid pixels and free the slot.
    fn sweep_dead_cells(&mut self) {
        let threshold = self.p.death_area_threshold;
        let to_clear: Vec<u32> = self.cells.iter()
            .filter(|c| c.id > 0 && c.dying && c.area < threshold)
            .map(|c| c.id)
            .collect();

        if to_clear.is_empty() { return; }

        for sigma in &to_clear {
            for px in self.grid.iter_mut() {
                if *px == *sigma { *px = 0; }
            }
            let c = &mut self.cells[*sigma as usize];
            c.alive = false;
            c.dying = false;
            c.area = 0;
            c.perimeter = 0;
            c.target_area = 0;
            c.target_perimeter = 0;
        }

        // Rebuild stats and boundary after bulk grid edits
        self.recompute_stats();
        self.boundary = build_boundary(&self.grid, self.p.grid_w, self.p.grid_h);
    }

    // ── Birth ─────────────────────────────────────────────────────────────────

    /// With probability `birth_prob`, place a new 2×2 cell in a free medium
    /// spot and assign it a sigma from a dead slot (or a new slot).
    fn maybe_birth_cell(&mut self) {
        let prob = self.p.birth_prob;
        if prob <= 0.0 { return; }
        if !self.rng.gen_bool(prob.min(1.0)) { return; }

        let W = self.p.grid_w;
        let H = self.p.grid_h;

        // Collect top-left corners of free 2×2 blocks (not wall, all medium)
        let mut free: Vec<(usize, usize)> = Vec::new();
        for r in 0..H.saturating_sub(1) {
            for c in 0..W.saturating_sub(1) {
                if self.wall.contains(r,   c)
                || self.wall.contains(r+1, c)
                || self.wall.contains(r,   c+1)
                || self.wall.contains(r+1, c+1)
                { continue; }
                if self.grid[r*W+c]       == 0
                && self.grid[(r+1)*W+c]   == 0
                && self.grid[r*W+c+1]     == 0
                && self.grid[(r+1)*W+c+1] == 0
                {
                    free.push((r, c));
                }
            }
        }
        if free.is_empty() { return; }

        let idx = self.rng.gen_range(0..free.len());
        let (br, bc) = free[idx];

        // Find a dead slot (id > 0, !alive) or push a new one
        let sigma: u32 = self.cells.iter()
            .position(|c| c.id > 0 && !c.alive)
            .map(|i| i as u32)
            .unwrap_or_else(|| {
                let new_id = self.cells.len() as u32;
                self.cells.push(CellState::new(new_id, 0, 0, 0, 0));
                new_id
            });

        let ta = self.p.target_area;
        let tp = self.p.target_perim;
        {
            let cell = &mut self.cells[sigma as usize];
            cell.alive = true;
            cell.dying = false;
            cell.target_area = ta;
            cell.target_perimeter = tp;
        }

        // Place 2×2 pixels
        self.grid[br*W+bc]       = sigma;
        self.grid[(br+1)*W+bc]   = sigma;
        self.grid[br*W+bc+1]     = sigma;
        self.grid[(br+1)*W+bc+1] = sigma;

        // Update stats and boundary for the new cell and its neighbourhood
        self.recompute_stats();
        self.boundary = build_boundary(&self.grid, self.p.grid_w, self.p.grid_h);
    }

    // ── Convenience ───────────────────────────────────────────────────────────

    /// Instantly kill a specific cell (identified by sigma).
    /// Its pixels are immediately replaced by medium.
    pub fn kill_cell(&mut self, sigma: u32) {
        if sigma == 0 || sigma as usize >= self.cells.len() { return; }
        for px in self.grid.iter_mut() {
            if *px == sigma { *px = 0; }
        }
        let c = &mut self.cells[sigma as usize];
        c.alive = false;
        c.dying = false;
        c.area = 0;
        c.perimeter = 0;
        c.target_area = 0;
        c.target_perimeter = 0;
        self.recompute_stats();
        self.boundary = build_boundary(&self.grid, self.p.grid_w, self.p.grid_h);
    }

    /// Count currently living cells (alive && !dying).
    pub fn n_living(&self) -> usize {
        self.cells.iter().filter(|c| c.id > 0 && c.alive && !c.dying).count()
    }

    /// Count currently dying cells.
    pub fn n_dying(&self) -> usize {
        self.cells.iter().filter(|c| c.id > 0 && c.dying).count()
    }
}
