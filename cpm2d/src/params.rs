use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
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
    pub lambda_iso: f64,

    // Surface tensions
    pub j_cell_medium: f64,
    pub j_cell_cell: f64,
    // Wall param
    pub wall_inset: usize,
    // Voronoi init
    pub j_wall: f64,          // surface tension between cell and wall
    pub poisson_min_dist: f64,
    // Monte Carlo
    pub temperature: f64,
    /// None → grid_w * grid_h attempts per MCS
    pub mcs_per_step: Option<usize>,
    pub total_steps: usize,

    // Output
    pub png_dir: String,
    pub frames_dir: String,
    pub png_every: usize,
    pub console_every: usize,
    pub save_every: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            grid_w: 60,
            grid_h: 60,
            n_cells: 4*4,
            target_area: 100,
            target_perim: (100f64 * 4.0 * PI).sqrt() as i64,
            lambda_area: 1.0,
            lambda_perim: 0.50,
            lambda_iso: 0.1,
            j_cell_medium: 8.0,
            j_cell_cell: 10.0,
            j_wall: 20f64,          // surface tension between cell and wall
            poisson_min_dist: 10f64, //0.8 * (60f64 * 60f64 / (3 * 3) as f64).sqrt(),
            wall_inset: 5usize,
            temperature: 2.0,
            mcs_per_step: None,
            total_steps: 2000,
            png_dir: "frames".to_string(),
            frames_dir: "states".to_string(),
            png_every: 10,
            console_every: 10,
            save_every: 50,
        }
    }
}