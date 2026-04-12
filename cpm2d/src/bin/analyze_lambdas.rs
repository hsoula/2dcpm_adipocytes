//! cargo run --bin analyse_lambdas -- "sweep_states/states_mcs*.json" --out results/sweep.csv

use clap::Parser;
use glob::glob;
use cpm2d::{grid::Cpm2d, stats::compute_stats};
use serde::Serialize;
use std::fs;
use std::path::Path;
use cpm2d::boundary::build_boundary;
use cpm2d::stats::{average_boundary_convexity, boundary_profiles, convex_hull_areas};

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
    target_area:  i32,
    n_cells:      i32,
    run :           i32,
    // cell stats  (mirrors CellStats fields, flattened)
    cell_id:       u32,
    area:          i64,
    perimeter:     i64,
    com_x:         f64,
    com_y:         f64,
    avg_convexity: f64,
    convex_hull_area: f64,
    solidity : f64,
}

/// Extract (lambda_area, lambda_perim, lambda_iso) from a filename like
/// "states_mcs000300_1.00_0.10_0.05.json"
fn parse_lambdas(path: &str) -> Option<(f64, f64, f64,i32, i32, i32,i32)> {
    // Strip directory and extension, then take the last three '_'-separated tokens
    let stem = std::path::Path::new(path)
        .file_stem()?          // "states_mcs000300_1.00_0.10_0.05"
        .to_str()?;

    let parts: Vec<&str> = stem.split('_').collect();
    if parts.len() < 8 { return None; }

    let n = parts.len();
    let la = parts[n - 7].parse().ok()?;
    let lp = parts[n - 6].parse().ok()?;
    let li = parts[n - 5].parse().ok()?;
    let vo :i32 = parts[n - 4].parse().ok()?;
    let area:i32 = parts[n - 3].parse().ok()?;
    let n_cells :i32 = parts[n - 2].parse().ok()?;
    let run = parts[n - 1].parse().ok()?;
    Some((la, lp, li, vo, area, n_cells,run))
}

fn main() {
    let cli = Cli::parse();
    let path = Path::new(&cli.out);
    //fs::create_dir_all(path).unwrap();

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

    let mut wtr = csv::Writer::from_path(format!("{}/out.csv",cli.out)).unwrap();
    let mut ok  = 0usize;
    let mut skip = 0usize;

    for path in &files {
        let (la, lp, li, vo, area, n_cells,run) = match parse_lambdas(path) {
            Some(v) => v,
            None    => {
                eprintln!("  skipping (cannot parse lambdas): {path}");
                skip += 1;
                continue;
            }
        };

        let sim   = Cpm2d::load_state(path);
        let (stats,edges) = compute_stats(&sim);
        let boundary = build_boundary(&sim.grid, sim.p.grid_w, sim.p.grid_h);
        let convexity   = average_boundary_convexity(&sim, &boundary);
        let hull_areas  = convex_hull_areas(&sim, &boundary);
        let profile = boundary_profiles(&sim);
        
        for (k,s) in stats.iter().enumerate() {
            wtr.serialize(CsvRow {
                mcs:          sim.mcs,
                lambda_area:  la,
                lambda_perimeter: lp,
                lambda_iso:   li,
                target_area:  area,
                n_cells:      n_cells,
                run:          run,
                cell_id:      s.cell_id,
                area:         s.area,
                perimeter:    s.perimeter,
                com_x:        s.com_x,
                com_y:        s.com_y,
                avg_convexity:    convexity[k],
                convex_hull_area: hull_areas[k],
                // solidity = actual area / hull area  (1.0 = perfectly convex)
                solidity: if hull_areas[k] > 0.0 {
                    sim.cells[k].area as f64 / hull_areas[k]
                } else { 0.0 },
            }).unwrap();
        }
        ok += 1;
    }
    wtr.flush().unwrap();
    println!("Written {ok} snapshots ({skip} skipped) → {}", cli.out);
}