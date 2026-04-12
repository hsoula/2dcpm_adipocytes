//! cargo run --bin debug_mcs -- state_stable.json --out debug_energies.csv

use clap::Parser;
use cpm2d::{grid::Cpm2d, wall::Wall};
use serde::Serialize;
use std::fs;
use cpm2d::fliprecord::{FlipRecord};
#[derive(Parser)]
struct Cli {
    input: String,
    #[arg(long, default_value = "debug_energies.csv")]
    out: String,
    /// How many rings to shrink the wall inward
    #[arg(long, default_value = "1")]
    wall_shrink: usize,
}

fn main() {
    let cli = Cli::parse();
    fs::create_dir_all("debug").unwrap();

    let mut sim = Cpm2d::load_state(&cli.input);

    // Shrink wall inward to create confinement pressure
    // for _ in 0..cli.wall_shrink {
    //     let freed = sim.wall.expand_inward();  // shrink = expand inward
    //     for (r, c) in freed {
    //         // assign freed-from-wall pixel to nearest cell
    //         sim.grid[r * sim.p.grid_w + c] = sim.nearest_cell(r, c);
    //     }
    // }
    sim.recompute_stats();
    println!("Wall shrunk by {} rings. Running one instrumented MCS…", cli.wall_shrink);

    let mut wtr     = csv::Writer::from_path(&cli.out).unwrap();
    let mut attempt = 0usize;

    let mcs_size = sim.mcs_size;
    for _ in 0..mcs_size {
        if let Some(record) = sim.attempt_instrumented(attempt) {
            wtr.serialize(record).unwrap();
        }
        attempt += 1;
    }

    wtr.flush().unwrap();
    println!("Written {} flip records → {}", attempt, cli.out);
}