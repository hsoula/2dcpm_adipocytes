//! Extract demography statistics from a simulate_live output directory.
//!
//! Reads all state JSON files to build a population time-series, and
//! copies the already-written events.csv through as-is.
//!
//! Outputs:
//!   <out>/population.csv   — one row per saved MCS snapshot
//!   (events.csv is written directly by simulate_live)
//!
//! Usage:
//!   cargo run --bin extract_demography -- --dir data/live/ --out results/live/

use std::fs;
use clap::Parser;
use glob::glob;
use cpm2d::grid::SaveState;

#[derive(Parser)]
struct Cli {
    /// Directory containing the simulate_live output (frames/ sub-dir with JSONs)
    #[arg(long, default_value = "data/live/")]
    dir: String,
    /// Output directory for CSVs
    #[arg(long, default_value = "results/live/")]
    out: String,
}

fn main() {
    let cli = Cli::parse();
    fs::create_dir_all(&cli.out).expect("cannot create output dir");

    // ── Population CSV ────────────────────────────────────────────────────────
    let pattern = format!("{}/frames/state_mcs*.json", cli.dir.trim_end_matches('/'));
    let mut files: Vec<String> = glob(&pattern)
        .expect("invalid glob pattern")
        .filter_map(|e| e.ok())
        .map(|p| p.to_string_lossy().into_owned())
        .collect();

    if files.is_empty() {
        eprintln!("No state JSON files found at: {}", pattern);
        std::process::exit(1);
    }
    files.sort();
    println!("Found {} state snapshots", files.len());

    let pop_path   = format!("{}/population.csv",  cli.out.trim_end_matches('/'));
    let areas_path = format!("{}/cell_areas.csv",  cli.out.trim_end_matches('/'));

    let mut pop_wtr = csv::Writer::from_path(&pop_path).expect("cannot create population.csv");
    pop_wtr.write_record(&["mcs", "n_living", "n_dying", "n_dead", "total_slots",
                            "mean_area_living", "mean_target_area_living"])
        .unwrap();

    // cell_areas.csv — one row per (mcs, cell) for every living or dying cell.
    // Grouping by mcs gives the area distribution at that snapshot.
    // Grouping by sigma gives the full area time-course up to and including death.
    let mut area_wtr = csv::Writer::from_path(&areas_path).expect("cannot create cell_areas.csv");
    area_wtr.write_record(&["mcs", "sigma", "area", "target_area", "status", "birth_mcs"])
        .unwrap();

    for path in &files {
        let json = fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
        let state: SaveState = serde_json::from_str(&json)
            .unwrap_or_else(|e| panic!("cannot parse {path}: {e}"));

        let living: Vec<_> = state.cells.iter()
            .filter(|c| c.id > 0 && c.alive && !c.dying)
            .collect();
        let dying: Vec<_> = state.cells.iter()
            .filter(|c| c.id > 0 && c.dying)
            .collect();
        let dead = state.cells.iter()
            .filter(|c| c.id > 0 && !c.alive && !c.dying)
            .count();

        let mean_area = if living.is_empty() { 0.0 } else {
            living.iter().map(|c| c.area as f64).sum::<f64>() / living.len() as f64
        };
        let mean_target = if living.is_empty() { 0.0 } else {
            living.iter().map(|c| c.target_area as f64).sum::<f64>() / living.len() as f64
        };

        pop_wtr.write_record(&[
            state.mcs.to_string(),
            living.len().to_string(),
            dying.len().to_string(),
            dead.to_string(),
            (state.cells.len().saturating_sub(1)).to_string(),
            format!("{:.2}", mean_area),
            format!("{:.2}", mean_target),
        ]).unwrap();

        // Per-cell area rows for living and dying cells
        for c in living.iter().chain(dying.iter()) {
            let status = if c.dying { "dying" } else { "living" };
            area_wtr.write_record(&[
                state.mcs.to_string(),
                c.id.to_string(),
                c.area.to_string(),
                c.target_area.to_string(),
                status.to_string(),
                c.birth_mcs.to_string(),
            ]).unwrap();
        }
    }

    pop_wtr.flush().unwrap();
    area_wtr.flush().unwrap();
    println!("Population CSV → {}", pop_path);
    println!("Cell areas CSV → {}", areas_path);

    // ── Copy / verify events.csv ──────────────────────────────────────────────
    let src_events = format!("{}/events.csv", cli.dir.trim_end_matches('/'));
    let dst_events = format!("{}/events.csv", cli.out.trim_end_matches('/'));
    if std::path::Path::new(&src_events).exists() {
        fs::copy(&src_events, &dst_events)
            .unwrap_or_else(|e| panic!("cannot copy events.csv: {e}"));
        // Count rows
        let content = fs::read_to_string(&dst_events).unwrap();
        let n_events = content.lines().count().saturating_sub(1); // minus header
        println!("Events CSV   → {}  ({} events)", dst_events, n_events);
    } else {
        eprintln!("Warning: no events.csv found at {}", src_events);
    }

    println!("Done.");
}
