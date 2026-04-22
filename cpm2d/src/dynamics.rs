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
//!
//! `step_demography` returns a `Vec<DemographyEvent>` — one entry per birth or
//! confirmed death that occurred this MCS.  The caller can use this to write an
//! event log without re-scanning the cell list.

use std::f64::consts::PI;
use rand::Rng;
use crate::boundary::build_boundary;
use crate::cellstate::CellState;
use crate::grid::Cpm2d;

// ── Event types ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, serde::Serialize)]
pub enum EventKind {
    Birth,
    Death,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct DemographyEvent {
    pub kind:         EventKind,
    pub sigma:        u32,
    pub mcs:          usize,
    /// Actual area of the cell when the event fired
    pub area_at_event: i64,
    /// MCS when this cell was born (same as mcs for births)
    pub birth_mcs:    usize,
    /// mcs - birth_mcs  (0 for births)
    pub lifetime_mcs: usize,
}

// ── Cpm2d impl ────────────────────────────────────────────────────────────────

impl Cpm2d {
    // ── Public entry point ────────────────────────────────────────────────────

    /// Run one demography step.  Call after `run_mcs`.
    /// Returns all birth and death events that occurred this step.
    pub fn step_demography(&mut self) -> Vec<DemographyEvent> {
        let mut events = Vec::new();
        self.grow_cells();
        self.maybe_kill_cells();
        self.sweep_dead_cells(&mut events);
        self.maybe_birth_cell(&mut events);
        events
    }

    // ── Growth ────────────────────────────────────────────────────────────────

    /// Increase target_area (and matching target_perimeter) for every live,
    /// non-dying cell by `growth_rate` pixels.
    fn grow_cells(&mut self) {
        let rate = self.p.growth_rates as i64;
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
    /// Appends a Death event for each cell removed.
    fn sweep_dead_cells(&mut self, events: &mut Vec<DemographyEvent>) {
        let threshold = self.p.death_area_threshold;
        let mcs = self.mcs;

        let to_clear: Vec<(u32, i64, usize)> = self.cells.iter()
            .filter(|c| c.id > 0 && c.dying && c.area < threshold)
            .map(|c| (c.id, c.area, c.birth_mcs))
            .collect();

        if to_clear.is_empty() { return; }

        for &(sigma, area, birth_mcs) in &to_clear {
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

            events.push(DemographyEvent {
                kind: EventKind::Death,
                sigma,
                mcs,
                area_at_event: area,
                birth_mcs,
                lifetime_mcs: mcs.saturating_sub(birth_mcs),
            });
        }

        self.recompute_stats();
        self.boundary = build_boundary(&self.grid, self.p.grid_w, self.p.grid_h);
    }

    // ── Birth ─────────────────────────────────────────────────────────────────

    /// With probability `birth_prob`, place a new 2×2 cell in a free medium
    /// spot.  Appends a Birth event on success.
    fn maybe_birth_cell(&mut self, events: &mut Vec<DemographyEvent>) {
        let prob = self.p.birth_prob;
        if prob <= 0.0 { return; }
        if !self.rng.gen_bool(prob.min(1.0)) { return; }

        let W = self.p.grid_w;
        let H = self.p.grid_h;

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
        let mcs = self.mcs;
        {
            let cell = &mut self.cells[sigma as usize];
            cell.alive = true;
            cell.dying = false;
            cell.target_area = ta;
            cell.target_perimeter = tp;
            cell.birth_mcs = mcs;
        }

        self.grid[br*W+bc]       = sigma;
        self.grid[(br+1)*W+bc]   = sigma;
        self.grid[br*W+bc+1]     = sigma;
        self.grid[(br+1)*W+bc+1] = sigma;

        self.recompute_stats();
        self.boundary = build_boundary(&self.grid, self.p.grid_w, self.p.grid_h);

        events.push(DemographyEvent {
            kind: EventKind::Birth,
            sigma,
            mcs,
            area_at_event: 4,
            birth_mcs: mcs,
            lifetime_mcs: 0,
        });
    }

    // ── Convenience ───────────────────────────────────────────────────────────

    /// Instantly kill a specific cell (identified by sigma).
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
