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

use std::env::join_paths;
use std::f64;
use std::fs;

use image::{ImageBuffer, Rgb};
use rand::prelude::*;
use serde::{Deserialize, Serialize};

// ─── Parameters ──────────────────────────────────────────────────────────────

use std::f64::consts::PI;
use cpm2d::params::Params;
use cpm2d::grid::Cpm2d;
use cpm2d::render::{save_png, console_print};
/// Represents a single cell in the CPM lattice
pub struct Cell {
    pub area: f64,          // current pixel count
    pub perimeter: f64,     // current boundary length
    pub target_area: f64,   // A₀
    pub target_perim: f64,  // P₀
    pub hull_area: f64,      // convex hull area (method 1 only)
    pub boundary_pixels: Vec<(i32, i32, u8)>, // (x, y, neighbor_count)
}

impl Cell {
    /// Number of same-cell neighbors of a boundary pixel (0..8)
    pub fn neighbor_count(pixels: &[(i32,i32)], x: i32, y: i32) -> u8 {
        let offsets = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)];
        offsets.iter().filter(|&(dx,dy)| pixels.contains(&(x+dx, y+dy))).count() as u8
    }
}


// ─── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut sim = if args.len() > 1 { let path = args.get(1).expect("--load requires a path argument");
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

    console_print(&sim);
    save_png(&sim, None);

    while sim.mcs < end_mcs {
        sim.run_mcs();

        if sim.mcs % p.console_every == 0 {
            console_print(&sim);
        }
        if sim.mcs % p.png_every == 0 {
            save_png(&sim, None);
        }
        if sim.mcs % p.save_every == 0 {
            sim.save_state(None);
        }
    }

    println!("\nDone. Final MCS={}", sim.mcs);
    sim.save_state(Some("states/state_final.json"));
    save_png(&sim, Some(&format!("{}/final.png", p.png_dir)));
}