//! Per-cell analysis over a 3-D CPM simulation run.
//!
//! Reads every `state_mcs*.json` file in the given directory (sorted by MCS)
//! and writes two CSV files:
//!
//!   cells.csv           — one row per alive cell per snapshot
//!                         mcs, sigma, volume, surface, target_volume, target_surface,
//!                         lipid, dying, birth_mcs, vol_ratio
//!
//!   population.csv      — one row per snapshot, population-level summary
//!                         mcs, n_cells, mean_vol, std_vol, min_vol, max_vol,
//!                         mean_lipid, total_lipid
//!
//! Usage
//! -----
//!   cargo run --bin analyze -- --dir data/sim3d_42_g1_d0.02_b0.01/

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

    files.sort();

    if files.is_empty() {
        eprintln!("No state_mcs*.json files found in {:?}", cli.dir);
        std::process::exit(1);
    }

    println!("Found {} snapshots in {:?}", files.len(), cli.dir);

    // ── Output buffers ────────────────────────────────────────────────────────
    let mut cell_rows = vec![
        "mcs,sigma,volume,surface,target_volume,target_surface,lipid,dying,birth_mcs,vol_ratio"
            .to_string(),
    ];
    let mut pop_rows = vec![
        "mcs,n_cells,mean_vol,std_vol,min_vol,max_vol,mean_lipid,total_lipid".to_string(),
    ];

    // ── Parse each snapshot ───────────────────────────────────────────────────
    for path in &files {
        let json = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {:?}: {e}", path));
        let state: SaveState = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("cannot parse {:?}: {e}", path));

        let mcs = state.mcs;

        // All non-medium cells that are alive (includes dying — they still exist)
        let alive: Vec<_> = state.cells.iter()
            .filter(|c| c.id > 0 && c.alive)
            .collect();

        // Per-cell rows
        for c in &alive {
            let vol_ratio = if c.target_volume > 0 {
                c.volume as f64 / c.target_volume as f64
            } else {
                0.0
            };
            cell_rows.push(format!(
                "{},{},{},{},{},{},{:.6},{},{},{:.4}",
                mcs, c.id,
                c.volume, c.surface,
                c.target_volume, c.target_surface,
                c.lipid,
                c.dying as u8,
                c.birth_mcs,
                vol_ratio,
            ));
        }

        // Population summary (non-dying cells for volume stats)
        let active: Vec<_> = alive.iter().filter(|c| !c.dying).collect();
        if active.is_empty() {
            pop_rows.push(format!("{},0,0,0,0,0,0,0", mcs));
        } else {
            let n         = active.len() as f64;
            let vols: Vec<f64> = active.iter().map(|c| c.volume as f64).collect();
            let mean_vol  = vols.iter().sum::<f64>() / n;
            let std_vol   = (vols.iter().map(|v| (v - mean_vol).powi(2)).sum::<f64>() / n).sqrt();
            let min_vol   = vols.iter().cloned().fold(f64::INFINITY,     f64::min) as i64;
            let max_vol   = vols.iter().cloned().fold(f64::NEG_INFINITY, f64::max) as i64;
            let total_lip = active.iter().map(|c| c.lipid).sum::<f64>();
            let mean_lip  = total_lip / n;
            pop_rows.push(format!(
                "{},{},{:.2},{:.2},{},{},{:.4},{:.4}",
                mcs, active.len(),
                mean_vol, std_vol, min_vol, max_vol,
                mean_lip, total_lip,
            ));
        }

        print!("  MCS {:6}  cells={}\r", mcs, alive.len());
    }
    println!();

    // ── Write CSVs ────────────────────────────────────────────────────────────
    let cell_path = cli.dir.join("cells.csv");
    let pop_path  = cli.dir.join("population.csv");
    fs::write(&cell_path, cell_rows.join("\n") + "\n").expect("cannot write cells.csv");
    fs::write(&pop_path,  pop_rows.join("\n")  + "\n").expect("cannot write population.csv");
    println!("→ {:?}  ({} cell rows)", cell_path,  cell_rows.len() - 1);
    println!("→ {:?}  ({} mcs rows)",  pop_path,   pop_rows.len()  - 1);
}
