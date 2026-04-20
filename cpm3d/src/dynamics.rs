//! Cell demography: growth, death, and birth.
//!
//! Call [`Cpm3d::step_demography`] once per MCS (after `run_mcs`) to evolve
//! the cell population:
//!   - Living cells grow linearly (target_area += growth_rate each MCS).
//!   - Each living cell dies with probability `death_prob` per MCS; death is
//!     implemented by setting target_area = 0 so the CPM gradually erodes the
//!     cell.  Once its actual area drops below `death_volume_threshold` the slot
//!     is freed.
//!   - With probability `birth_probability` a new 2×2 cell is placed in a free
//!     medium region; dead slots are reused before the cells vec is extended.
//!
//! `step_demography` returns a `Vec<DemographyEvent>` — one entry per birth or
//! confirmed death that occurred this MCS.  The caller can use this to write an
//! event log without re-scanning the cell list.

use std::f64::consts::PI;
use rand::Rng;
use rand::prelude::*;
use rand_distr::Normal;
use crate::cellstate::CellState;
use crate::grid::{Cpm3d, compute_surface_from_volume};

// ── Event types ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, serde::Serialize)]
pub enum EventKind {
    Birth,
    Dying,
    Dead,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct DemographyEvent {
    pub kind:         EventKind,
    pub sigma:        u32,
    pub mcs:          usize,
    /// Actual area of the cell when the event fired
    pub volume_at_event: i64,
    /// MCS when this cell was born (same as mcs for births)
    pub birth_mcs:    usize,
    /// mcs - birth_mcs  (0 for births)
    pub lifetime_mcs: usize,
}

// ── Cpm2d impl ────────────────────────────────────────────────────────────────

impl Cpm3d {
    // ── Public entry point ────────────────────────────────────────────────────

    /// Run one demography step.  Call after `run_mcs`.
    /// Returns all birth and death events that occurred this step.
    pub fn step_demography(&mut self) -> Vec<DemographyEvent> {
        let mut events = Vec::new();
        self.grow_cells();
        self.maybe_kill_cells(&mut events);
        self.sweep_dead_cells(&mut events);
        self.maybe_birth_cell(&mut events);
        events
    }

    // ── Growth ────────────────────────────────────────────────────────────────

    /// Increase target_area (and matching target_perimeter) for every live,
    /// non-dying cell by `growth_rate` pixels.
    fn grow_cells(&mut self) {
        let rate = self.p.growth_rate as f64;
        //if rate == 0 { return; }
        for cell in self.cells.iter_mut() {
            if cell.id == 0 || !cell.alive || cell.dying { continue; }
            cell.lipid += rate;
            cell.target_volume  = cell.lipid as i64 + 4i64;
            cell.target_surface = compute_surface_from_volume(cell.target_volume as f64) as i64;
        }
    }

    // ── Death ─────────────────────────────────────────────────────────────────

    /// Mark each live cell as dying with probability `death_prob`.
    /// A dying cell has target_area = 0 so CPM pressure shrinks it.
    /// Appends a Dying event for each cell removed
    fn maybe_kill_cells(&mut self, events: &mut Vec<DemographyEvent>) {
        let prob = self.p.death_rate;
        if prob <= 0.0 { return; }
        let mcs = self.mcs;
        let to_kill: Vec<usize> = (1..self.cells.len())
            .filter(|&i| {
                let c = &self.cells[i];
                c.alive && !c.dying && self.rng.gen_bool(prob.min(1.0))
            })
            .collect();

        for i in to_kill {
            self.cells[i].dying = true;
            self.cells[i].target_volume = 0;
            self.cells[i].target_surface = 0;

            events.push(DemographyEvent {
                kind: EventKind::Dying,
                sigma: i as u32,
                mcs: mcs,
                volume_at_event: self.cells[i].volume,
                birth_mcs:self.cells[i].birth_mcs,
                lifetime_mcs: mcs.saturating_sub(self.cells[i].birth_mcs),
            });
        }
    }

    /// Remove dying cells whose actual area has dropped below
    /// `death_area_threshold`: zero all their grid pixels and free the slot.
    /// Appends a Dead event for each cell removed.
    fn sweep_dead_cells(&mut self, events: &mut Vec<DemographyEvent>) {
        let threshold = self.p.death_volume_threshold;
        let mcs = self.mcs;

        let to_clear: Vec<(u32, i64, usize)> = self.cells.iter()
            .filter(|c| c.id > 0 && c.dying && (c.volume as f64) < threshold)
            .map(|c| (c.id, c.volume, c.birth_mcs))
            .collect();

        if to_clear.is_empty() { return; }

        for &(sigma, area, birth_mcs) in &to_clear {
            for px in self.grid.iter_mut() {
                if *px == sigma { *px = 0; }
            }
            let c = &mut self.cells[sigma as usize];
            c.alive = false;
            c.dying = false;
            c.volume = 0;
            c.surface = 0;
            c.target_volume = 0;
            c.target_surface = 0;

            events.push(DemographyEvent {
                kind: EventKind::Dead,
                sigma,
                mcs,
                volume_at_event: area,
                birth_mcs,
                lifetime_mcs: mcs.saturating_sub(birth_mcs),
            });
        }

        self.recompute_stats();
    }

    // ── Birth ─────────────────────────────────────────────────────────────────

    /// With probability `birth_prob`, place a new 2×2 cell in a free medium
    /// spot.  Appends a Birth event on success.
    fn maybe_birth_cell(&mut self, events: &mut Vec<DemographyEvent>) {
        let prob = self.p.birth_rate;
        if prob <= 0.0 { return; }
        if !self.rng.gen_bool(prob.min(1.0)) { return; }

        let w = self.p.grid_w;
        let h = self.p.grid_h;
        let d = self.p.grid_d;

        let mut free: Vec<(usize, usize, usize)> = Vec::new();
        for r in 1..h.saturating_sub(1) {
            for c in 1..w.saturating_sub(1) {
                for b in 1..d.saturating_sub(1) {
                    if self.grid_get(r, c, b) == 0
                        && self.grid_get(r + 1, c, b) == 0
                        && self.grid_get(r, c + 1, b) == 0
                        && self.grid_get(r + 1, c + 1, b) == 0
                        && self.grid_get(r, c, b + 1) == 0
                        && self.grid_get(r + 1, c, b + 1) == 0
                        && self.grid_get(r, c + 1, b + 1) == 0
                        && self.grid_get(r + 1, c + 1, b + 1) == 0
                    {
                        free.push((r, c, b));
                    }
                }
            }
        }
        if free.is_empty() { return; }

        let idx = self.rng.gen_range(0..free.len());
        let (br, bc, bd) = free[idx];

        let sigma: u32 = self.cells.iter()
            .position(|c| c.id > 0 && !c.alive)
            .map(|i| i as u32)
            .unwrap_or_else(|| {
                let new_id = self.cells.len() as u32;
                self.cells.push(CellState::new(new_id, 0, 0));
                new_id
            });
        let normal_dist = if self.p.volume_sigma > 0.0 {
            Normal::new(self.p.target_volume as f64, self.p.volume_sigma).ok()
        } else {
            None
        };
        let ta = if let Some(ref dist) = normal_dist {
            let v = dist.sample(&mut self.rng).round() as i64;
            v.max(1)
        } else {
            self.p.target_volume
        };
        let tp = compute_surface_from_volume(ta as f64) as i64;
        let mcs = self.mcs;
        {
            let cell = &mut self.cells[sigma as usize];
            cell.alive = true;
            cell.dying = false;
            cell.target_volume = ta;
            cell.target_surface = tp;
            cell.birth_mcs = mcs;
        }
        self.grid_set(sigma, br, bc, bd);
        self.grid_set(sigma, br + 1, bc, bd);
        self.grid_set(sigma, br, bc + 1, bd);
        self.grid_set(sigma, br + 1, bc + 1, bd);
        self.grid_set(sigma, br, bc, bd + 1);
        self.grid_set(sigma, br + 1, bc, bd + 1);
        self.grid_set(sigma, br, bc + 1, bd + 1);
        self.grid_set(sigma, br + 1, bc + 1, bd + 1);


        self.recompute_stats();

        events.push(DemographyEvent {
            kind: EventKind::Birth,
            sigma,
            mcs,
            volume_at_event: 4,
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
        c.volume = 0;
        c.surface = 0;
        c.target_volume = 0;
        c.target_surface = 0;
        self.recompute_stats();
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
