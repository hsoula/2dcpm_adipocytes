use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Params {
    // Grid
    pub grid_w: usize,
    pub grid_h: usize,
    pub grid_d: usize,

    // Cells
    pub n_cells: usize,
    pub target_volume: i64,
    pub target_surface: i64,

    // Hamiltonian weights
    pub lambda_vol:   f64,
    pub lambda_surf:  f64,
    pub lambda_spher: f64,   // sphericity penalty (S³/36πV² → 1 for sphere)

    // Surface tensions
    pub j_cell_medium: f64,
    pub j_cell_cell:   f64,

    // Monte Carlo
    pub temperature: f64,
    /// None → grid_w * grid_h * grid_d attempts per MCS
    pub mcs_per_step: Option<usize>,
    pub total_steps: usize,

    // Output
    pub out_dir:    String,
    pub save_every: usize,
    pub png_every:  usize,   // slice PNGs every N MCS (0 = off)
}

impl Default for Params {
    fn default() -> Self {
        let tv = 125i64; // ≈ 5³ voxels
        // Surface of sphere with volume V: S = (36π)^{1/3} · V^{2/3}
        let ts = ((36.0 * PI).powf(1.0 / 3.0) * (tv as f64).powf(2.0 / 3.0)).round() as i64;
        Self {
            grid_w: 60,
            grid_h: 60,
            grid_d: 60,
            n_cells: 8,
            target_volume: tv,
            target_surface: ts,
            lambda_vol:   1.0,
            lambda_surf:  0.5,
            lambda_spher: 0.05,
            j_cell_medium: 8.0,
            j_cell_cell:  10.0,
            temperature: 2.0,
            mcs_per_step: None,
            total_steps: 2000,
            out_dir:    "data/sim3d/".to_string(),
            save_every: 100,
            png_every:  50,
        }
    }
}
