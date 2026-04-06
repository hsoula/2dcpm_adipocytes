//! cargo run --bin analyse -- states/state_mcs*.json --out results/stats.csv

use std::fs;
use clap::Parser;
use cpm2d::{grid::Cpm2d, stats::compute_stats};
use glob::glob;
#[derive(Parser)]
struct Cli {
    /// Glob pattern or explicit files, e.g. "states/state_mcs*.json"
    patterns: Vec<String>,
    #[arg(long, default_value = "results/stats.csv")]
    out: String,
}

fn main() {
    let cli = Cli::parse();
    fs::create_dir_all("results").unwrap();
    let mut wtr_stats  = csv::Writer::from_path("results/stats.csv").unwrap();
    let mut wtr_edges  = csv::Writer::from_path("results/neighbours.csv").unwrap();


    // Expand each argument as a glob pattern
    let mut files: Vec<String> = Vec::new();
    for pattern in &cli.patterns {
        for entry in glob(pattern).expect("invalid glob pattern") {
            match entry {
                Ok(path) => files.push(path.to_string_lossy().into_owned()),
                Err(e)   => eprintln!("glob error: {e}"),
            }
        }
    }

    if files.is_empty() {
        eprintln!("no files matched");
        return;
    }

    // Sort so MCS order is preserved
    files.sort();
    println!("Analysing {} files…", files.len());

    for path in &files {
        let sim = Cpm2d::load_state(path);
        let (stats, edges) = compute_stats(&sim);
        for s in stats { wtr_stats.serialize(&s).unwrap(); }
        for e in edges { wtr_edges.serialize(&e).unwrap(); }
    }

    wtr_stats.flush().unwrap();
    wtr_edges.flush().unwrap();
    println!("Written → {}", cli.out);
}