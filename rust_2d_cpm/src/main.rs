// dependencies in Cargo.toml:
// [dependencies]
// rand = "0.8"
// image = "0.24"
// serde = { version = "1.0", features = ["derive"] }
// bincode = "1.3"
// anyhow = "1.0"

use rand::Rng;
use rand::seq::SliceRandom;
use image::{ImageBuffer, Rgb};
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::path::Path;
use anyhow::Result;
use bincode::config;  // for configuration

// ------------------------------------------------------------
// Simulation parameters
// ------------------------------------------------------------
const GRID_SIZE: usize = 200;          // grid is GRID_SIZE x GRID_SIZE
const BLOCK_SIZE: usize = 10;          // each initial cell is BLOCK_SIZE x BLOCK_SIZE
const TEMPERATURE: f64 = 2.0;         // Monte Carlo temperature
const LAMBDA_AREA: f64 = 0.5;   // area constraint strength
const CELL_SIZE: usize = 30;           // each cell is CELL_SIZE x CELL_SIZE
const NUM_CELLS_X: usize = 2;          // number of cells horizontally
const NUM_CELLS_Y: usize = 2;           // number of cells vertically
const J_CC: f64 = 0.0;                 // adhesion between two cells (low → attractive)
const J_CM: f64 = 10.0;                // adhesion between cell and medium (high → repulsive)
const J_MM: f64 = 0.0;                 // adhesion between medium and medium
const MCS_PER_FRAME: usize = 1;      // Monte Carlo steps between outputs
const TOTAL_MCS: usize = 1;       // total simulation time

// ------------------------------------------------------------
// Helper: cell type (0 = medium, 1 = cell)
// ------------------------------------------------------------
fn cell_type(id: u32) -> usize {
    if id == 0 { 0 } else { 1 }
}

// ------------------------------------------------------------
// Main simulation structure
// ------------------------------------------------------------
#[derive(Serialize, Deserialize)]
struct CPM {
    grid: Vec<Vec<u32>>,
    cell_areas: Vec<u32>,
    target_areas: Vec<u32>,
    #[serde(skip)]
    j: [[f64; 2]; 2],
    lambda_area: f64,
    temperature: f64,
    #[serde(skip)]
    rng: rand::rngs::ThreadRng,
}

impl CPM {
    /// Create a new simulation with square cells packed in a grid.
    fn new(grid_size: usize, block_size: usize) -> Self {
        let num_cells_per_side = grid_size / block_size;
        let num_cells = num_cells_per_side * num_cells_per_side;
        let mut grid = vec![vec![0; grid_size]; grid_size];
        let mut cell_areas = vec![0; num_cells + 1]; // index 0 is medium, cells start at 1
        let mut target_areas = vec![0; num_cells + 1];

        // Assign cell IDs and count initial areas
        for i in 0..num_cells_per_side {
            for j in 0..num_cells_per_side {
                let cell_id = (i * num_cells_per_side + j + 1) as u32;
                let x0 = i * block_size;
                let y0 = j * block_size;
                for dx in 0..block_size {
                    for dy in 0..block_size {
                        grid[x0 + dx][y0 + dy] = cell_id;
                        cell_areas[cell_id as usize] += 1;
                    }
                }
                target_areas[cell_id as usize] = 1900 as u32; //(block_size * block_size) as u32;
            }
        }

        // Adhesion matrix: [medium, cell] x [medium, cell]
        let j = [
            [J_MM, J_CM],
            [J_CM, J_CC],
        ];

        CPM {
            grid,
            cell_areas,
            target_areas,
            j,
            lambda_area: LAMBDA_AREA,
            temperature: TEMPERATURE,
            rng: rand::thread_rng(),
        }
    }
    fn new_cluster(grid_size: usize, cell_size: usize, num_cells_x: usize, num_cells_y: usize) -> Self {
        let total_cells = num_cells_x * num_cells_y;
        let mut grid = vec![vec![0; grid_size]; grid_size];
        let mut cell_areas = vec![0; total_cells + 1];
        let mut target_areas = vec![0; total_cells + 1];

        // Compute offset to center the cluster
        let total_width = num_cells_x * cell_size;
        let total_height = num_cells_y * cell_size;
        let offset_x = (grid_size - total_width) / 2;
        let offset_y = (grid_size - total_height) / 2;

        // Place each cell as a solid square
        for i in 0..num_cells_x {
            for j in 0..num_cells_y {
                let cell_id = (j * num_cells_x + i + 1) as u32; // row‑major order, IDs start at 1
                let x0 = offset_x + i * cell_size;
                let y0 = offset_y + j * cell_size;
                for dx in 0..cell_size {
                    for dy in 0..cell_size {
                        grid[x0 + dx][y0 + dy] = cell_id;
                    }
                }
                let area = (cell_size * cell_size) as u32;
                cell_areas[cell_id as usize] = area;
                target_areas[cell_id as usize] = area;
            }
        }

        // Adhesion matrix: [medium, cell] x [medium, cell]
        let j = [
            [J_MM, J_CM],
            [J_CM, J_CC],
        ];

        CPM {
            grid,
            cell_areas,
            target_areas,
            j,
            lambda_area: LAMBDA_AREA,
            temperature: TEMPERATURE,
            rng: rand::thread_rng(),
        }
    }

    /// Compute the energy change if we copy the source cell into the target pixel.
    fn delta_energy(&self, src_id: u32, tgt_id: u32, sx: usize, sy: usize, tx: usize, ty: usize) -> f64 {
        let mut delta = 0.0;

        // --- Adhesion part ---
        // Remove old contributions from the target pixel
        for (nx, ny) in self.neighbors(tx, ty) {
            let neigh_id = self.grid[nx][ny];
            if neigh_id != tgt_id {
                let tgt_type = cell_type(tgt_id);
                let neigh_type = cell_type(neigh_id);
                delta -= self.j[tgt_type][neigh_type];
            }
        }
        // Add new contributions from the target pixel (now becomes src_id)
        for (nx, ny) in self.neighbors(tx, ty) {
            let neigh_id = self.grid[nx][ny];
            if neigh_id != src_id {
                let src_type = cell_type(src_id);
                let neigh_type = cell_type(neigh_id);
                delta += self.j[src_type][neigh_type];
            }
        }

        // --- Area part (only for cells, medium has no area) ---
        if src_id != 0 {
            let area_src = self.cell_areas[src_id as usize];
            let target_src = self.target_areas[src_id as usize];
            let new_area_src = area_src + 1;
            delta += self.lambda_area * (
                (new_area_src as f64 - target_src as f64).powi(2) -
                    (area_src as f64 - target_src as f64).powi(2)
            );
        }
        if tgt_id != 0 {
            let area_tgt = self.cell_areas[tgt_id as usize];
            let target_tgt = self.target_areas[tgt_id as usize];
            let new_area_tgt = area_tgt - 1;
            delta += self.lambda_area * (
                (new_area_tgt as f64 - target_tgt as f64).powi(2) -
                    (area_tgt as f64 - target_tgt as f64).powi(2)
            );
        }

        delta
    }

    /// Return an iterator over the 4‑neighbors of (x, y) that lie inside the grid.
    fn neighbors(&self, x: usize, y: usize) -> Vec<(usize, usize)> {
        let mut neigh = Vec::with_capacity(4);
        if x > 0 { neigh.push((x - 1, y)); }
        if x + 1 < GRID_SIZE { neigh.push((x + 1, y)); }
        if y > 0 { neigh.push((x, y - 1)); }
        if y + 1 < GRID_SIZE { neigh.push((x, y + 1)); }
        neigh
    }

    /// Perform one Monte Carlo step (N_pixel attempts).
    /// Returns the number of times a cell tried to expand into medium with positive ΔE (pressure events).
    fn monte_carlo_step(&mut self) -> usize {
        let n_pixels = GRID_SIZE * GRID_SIZE;
        // Create a vector of all pixel indices (flat index = x * GRID_SIZE + y)
        let mut indices: Vec<usize> = (0..n_pixels).collect();
        // Shuffle the indices randomly
        indices.shuffle(&mut self.rng);

        let mut pressure_events = 0;

        for flat_idx in indices {

            let sx = flat_idx / GRID_SIZE;
            let sy = flat_idx % GRID_SIZE;

            let neighbors = self.neighbors(sx, sy);

            if neighbors.is_empty() { continue; }
            let (tx, ty) = neighbors[self.rng.gen_range(0..neighbors.len())];
            let src_id = self.grid[sx][sy];
            let tgt_id = self.grid[tx][ty];
            if src_id == tgt_id { continue; }

            println!("sx: {:?}, sy: {:?}", sx, sy);
            println!("neighboiurs {:?}", neighbors);
            println!("tx: {:?}, ty: {:?}", tx, ty);
            println!("tgt_id: {:?}, src_id: {:?}", tgt_id, src_id);

            // Compute energy change
            let delta = self.delta_energy(src_id, tgt_id, sx, sy, tx, ty);
            println!("-----------------------delta: {:?}", delta);
            // Pressure measurement: cell trying to expand into medium with positive ΔE
            if src_id != 0 && tgt_id == 0 && delta > 0.0 {
                pressure_events += 1;
            }

            // Metropolis acceptance
            if delta <= 0.0 || self.rng.random::<f64>() < (-delta / self.temperature).exp() {
                // Accept: update grid and areas
                self.grid[tx][ty] = src_id;
                if src_id != 0 {
                    self.cell_areas[src_id as usize] += 1;
                }
                if tgt_id != 0 {
                    self.cell_areas[tgt_id as usize] -= 1;
                }
            }
        }
        pressure_events
    }

    // /// Save the current state (grid, areas, target areas) to a binary file.
    // fn save_state(&self, path: &Path) -> Result<()> {
    //     let file = File::create(path)?;
    //     let config = config::standard();
    //     bincode::encode_into_std_write(self, file, config)?;
    //     Ok(())
    // }
    //
    // /// Load a previously saved state.
    // fn load_state(path: &Path) -> Result<Self> {
    //     let file = File::open(path)?;
    //     let config = config::standard();
    //     let (state, _): (Self, _) = bincode::decode_from_std_read(file, config)?;
    //     // Re‑create the RNG (cannot be serialized)
    //     Ok(state)
    // }

    /// Render the grid as an RGB image.
    /// Uses a simple hash to assign a colour per cell ID.
    fn render(&self) -> ImageBuffer<Rgb<u8>, Vec<u8>> {
        let mut img = ImageBuffer::new(GRID_SIZE as u32, GRID_SIZE as u32);
        for x in 0..GRID_SIZE {
            for y in 0..GRID_SIZE {
                let id = self.grid[x][y];
                let color = if id == 0 {
                    Rgb([0, 0, 0]) // medium = black
                } else {
                    // Simple hash to produce a colour
                    let r = ((id * 123) % 255) as u8;
                    let g = ((id * 456) % 255) as u8;
                    let b = ((id * 789) % 255) as u8;
                    Rgb([r, g, b])
                };
                img.put_pixel(x as u32, y as u32, color);
            }
        }
        img
    }
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
fn main() -> Result<()> {
    //let mut sim = CPM::new(GRID_SIZE, BLOCK_SIZE);
    let mut sim = CPM::new_cluster(GRID_SIZE, CELL_SIZE, 1,1);
    // Output directories
    std::fs::create_dir_all("frames")?;
    std::fs::create_dir_all("snapshots")?;

    // Total pressure events over the whole simulation
    let mut total_pressure = 0;

    // Initial frame
    sim.render().save("frames/step_0000.png")?;
    //sim.save_state(Path::new("snapshots/step_0000.bin"))?;

    for step in 1..=TOTAL_MCS {
        // Perform one Monte Carlo step and accumulate pressure events
        let pressure_this_step = sim.monte_carlo_step();
        total_pressure += pressure_this_step;

        // Output every MCS_PER_FRAME steps
        if step % MCS_PER_FRAME == 0 {
            let filename = format!("frames/step_{:04}.png", step);
            sim.render().save(&filename)?;

            let snapshot = format!("snapshots/step_{:04}.bin", step);
         //   sim.save_state(Path::new(&snapshot))?;

            // Print progress and pressure rate
            let avg_pressure = total_pressure as f64 / step as f64;
            println!(
                "Step {}/{}: total pressure events = {}, avg per MCS = {:.2}",
                step, TOTAL_MCS, total_pressure, avg_pressure
            );
        }
    }

    // Final output
    sim.render().save("frames/final.png")?;
  //  sim.save_state(Path::new("snapshots/final.bin"))?;
    println!("Simulation finished. Total pressure events = {}", total_pressure);
    Ok(())
}