//! Volume analysis over a 3-D CPM simulation run.
//!
//! Reads every `state_mcs*.json` file in the given directory (sorted by MCS),
//! and writes two CSV files next to the snapshots:
//!
//!   volume_timeseries.csv   — one row per alive cell per snapshot
//!                             columns: mcs, sigma, volume, surface, target_volume
//!
//!   volume_distribution.csv — one row per snapshot, summary statistics
//!                             columns: mcs, n_cells, mean_vol, std_vol, min_vol, max_vol
//!
//! Usage
//! -----
//!   cargo run --bin analyze -- --dir data/sim3d/

use std::fs;
use std::path::PathBuf;
use clap::Parser;
use cpm3d::grid::SaveState;

#[derive(Parser)]
struct Cli {
    /// Directory containing state_mcs*.json snapshots.
    #[arg(long, default_value = "data/sim3d/")]
    dir: PathBuf,
}

fn main() {
    let cli = Cli::parse();

    // ── Collect and sort snapshot files ──────────────────────────────────────
    let mut files: Vec<PathBuf> = fs::read_dir(&cli.dir)
        .unwrap_or_else(|e| panic!("cannot open {:?}: {e}", cli.dir))
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.starts_with("state_mcs") && n.ends_with(".json"))
                .unwrap_or(false)
        })
        .collect();

    files.sort();  // lexicographic sort == MCS order since filenames are zero-padded

    if files.is_empty() {
        eprintln!("No state_mcs*.json files found in {:?}", cli.dir);
        std::process::exit(1);
    }

    println!("Found {} snapshots in {:?}", files.len(), cli.dir);

    // ── Output file handles ───────────────────────────────────────────────────
    let ts_path   = cli.dir.join("volume_timeseries.csv");
    let dist_path = cli.dir.join("volume_distribution.csv");

    let mut ts_rows   = vec!["mcs,sigma,volume,surface,target_volume".to_string()];
    let mut dist_rows = vec!["mcs,n_cells,mean_vol,std_vol,min_vol,max_vol".to_string()];

    // ── Parse each snapshot ───────────────────────────────────────────────────
    for path in &files {
        let json = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {:?}: {e}", path));
        let state: SaveState = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("cannot parse {:?}: {e}", path));

        let mcs = state.mcs;

        // Alive, non-medium cells
        let alive: Vec<_> = state.cells.iter()
            .filter(|c| c.id > 0 && c.alive)
            .collect();

        // Per-cell rows
        for c in &alive {
            ts_rows.push(format!("{},{},{},{},{}", mcs, c.id, c.volume, c.surface, c.target_volume));
        }

        // Distribution summary
        if alive.is_empty() {
            dist_rows.push(format!("{},0,0,0,0,0", mcs));
        } else {
            let vols: Vec<f64> = alive.iter().map(|c| c.volume as f64).collect();
            let n    = vols.len() as f64;
            let mean = vols.iter().sum::<f64>() / n;
            let std  = (vols.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n).sqrt();
            let min  = vols.iter().cloned().fold(f64::INFINITY, f64::min) as i64;
            let max  = vols.iter().cloned().fold(f64::NEG_INFINITY, f64::max) as i64;
            dist_rows.push(format!("{},{},{:.2},{:.2},{},{}", mcs, alive.len(), mean, std, min, max));
        }

        print!("  MCS {:6}  cells={}\r", mcs, alive.len());
    }
    println!();

    // ── Write CSVs ────────────────────────────────────────────────────────────
    fs::write(&ts_path,   ts_rows.join("\n")   + "\n").expect("cannot write timeseries csv");
    fs::write(&dist_path, dist_rows.join("\n") + "\n").expect("cannot write distribution csv");

    println!("→ {:?}", ts_path);
    println!("→ {:?}", dist_path);
}
