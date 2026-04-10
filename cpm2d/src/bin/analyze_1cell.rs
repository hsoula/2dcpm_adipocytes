//! cargo run --bin analyse -- "sweep_states/states_mcs*.json" --out results/sweep.csv

use clap::Parser;
use glob::glob;
use cpm2d::{grid::Cpm2d, stats::compute_stats};
use serde::Serialize;
use std::fs;

#[derive(Parser)]
struct Cli {
    patterns: Vec<String>,
    #[arg(long, default_value = "results/sweep.csv")]
    out: String,
}

/// One CSV row = one cell at one MCS snapshot with its lambda context.
#[derive(Serialize)]
struct CsvRow {
    // provenance
    mcs:           usize,
    lambda_area:   f64,
    lambda_perimeter:  f64,
    lambda_iso:    f64,
    // cell stats  (mirrors CellStats fields, flattened)
    cell_id:       u32,
    area:          i64,
    perimeter:     i64,
    com_x:         f64,
    com_y:         f64,
}

/// Extract (lambda_area, lambda_perim, lambda_iso) from a filename like
/// "states_mcs000300_1.00_0.10_0.05.json"
fn parse_lambdas(path: &str) -> Option<(f64, f64, f64, i32)> {
    // Strip directory and extension, then take the last three '_'-separated tokens
    let stem = std::path::Path::new(path)
        .file_stem()?          // "states_mcs000300_1.00_0.10_0.05"
        .to_str()?;

    let parts: Vec<&str> = stem.split('_').collect();
    if parts.len() < 4 { return None; }

    let n = parts.len();
    let la = parts[n - 4].parse().ok()?;
    let lp = parts[n - 3].parse().ok()?;
    let li = parts[n - 2].parse().ok()?;
    let run = parts[n - 1].parse().ok()?;
    Some((la, lp, li, run))
}

fn main() {
    let cli = Cli::parse();
    fs::create_dir_all("results").unwrap();

    // Expand globs
    let mut files: Vec<String> = Vec::new();
    for pattern in &cli.patterns {
        for entry in glob(pattern).expect("invalid glob pattern") {
            match entry {
                Ok(p)  => files.push(p.to_string_lossy().into_owned()),
                Err(e) => eprintln!("glob error: {e}"),
            }
        }
    }
    files.sort();   // chronological order within each lambda combo

    if files.is_empty() {
        eprintln!("no files matched");
        return;
    }
    println!("Analysing {} files…", files.len());

    let mut wtr = csv::Writer::from_path(&cli.out).unwrap();
    let mut ok  = 0usize;
    let mut skip = 0usize;

    for path in &files {
        let (la, lp, li, run) = match parse_lambdas(path) {
            Some(v) => v,
            None    => {
                eprintln!("  skipping (cannot parse lambdas): {path}");
                skip += 1;
                continue;
            }
        };

        let sim   = Cpm2d::load_state(path);
        let (stats,edges) = compute_stats(&sim);

        for s in stats {
            wtr.serialize(CsvRow {
                mcs:          sim.mcs,
                lambda_area:  la,
                lambda_perimeter: lp,
                lambda_iso:   li,
                cell_id:      s.cell_id,
                area:         s.area,
                perimeter:    s.perimeter,
                com_x:        s.com_x,
                com_y:        s.com_y,
            }).unwrap();
        }
        ok += 1;
    }
    wtr.flush().unwrap();
    println!("Written {ok} snapshots ({skip} skipped) → {}", cli.out);
}