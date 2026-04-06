use serde::{Deserialize, Serialize};

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
            n_cells: 9,
            target_area: 300,
            target_perim: 300f64.sqrt() as i64,
            lambda_area: 1.0,
            lambda_perim: 0.1,
            lambda_iso: 0.1,
            j_cell_medium: 10.0,
            j_cell_cell: 2.0,
            temperature: 20.0,
            mcs_per_step: None,
            total_steps: 1000,
            png_dir: "frames".to_string(),
            frames_dir: "states".to_string(),
            png_every: 10,
            console_every: 10,
            save_every: 50,
        }
    }
}