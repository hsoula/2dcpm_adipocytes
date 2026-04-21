//! Aggregate sweep results across directories named `*_s{seed}_g{grow}_d{death}_b{birth}`.
//!
//! Reads `volume_distribution.csv` and `events.csv` from each matching subdir
//! and writes two combined CSVs into the base directory:
//!
//!   sweep_volumes.csv  — volume_distribution rows with param columns prepended
//!   sweep_events.csv   — events rows with param columns prepended
//!
//! Run `analyze` on each subdir first to generate volume_distribution.csv.
//!
//! Usage
//! -----
//!   cargo run --bin collect_results -- --base-dir data/

use std::fs;
use std::path::PathBuf;
use clap::Parser;

#[derive(Parser)]
struct Cli {
    /// Parent directory containing sweep subdirectories.
    #[arg(long, default_value = "data/")]
    base_dir: PathBuf,
}

/// Parse `_s{seed}_g{grow}_d{death}_b{birth}` suffix from a directory name.
fn parse_params(dirname: &str) -> Option<(String, f64, f64, f64)> {
    let s_pos = dirname.rfind("_s")?;
    let rest  = &dirname[s_pos + 2..];
    let parts: Vec<&str> = rest.splitn(4, '_').collect();
    if parts.len() < 4 { return None; }
    let seed  = parts[0].to_string();
    let grow  = parts[1].strip_prefix('g')?.parse::<f64>().ok()?;
    let death = parts[2].strip_prefix('d')?.parse::<f64>().ok()?;
    let birth = parts[3].strip_prefix('b')?.parse::<f64>().ok()?;
    Some((seed, grow, death, birth))
}

fn main() {
    let cli = Cli::parse();
    let base = &cli.base_dir;

    let mut subdirs: Vec<PathBuf> = fs::read_dir(base)
        .unwrap_or_else(|e| panic!("cannot open {:?}: {e}", base))
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.is_dir())
        .collect();
    subdirs.sort();

    let mut vol_rows = vec![
        "seed,grow_rate,death_rate,birth_rate,mcs,n_cells,mean_vol,std_vol,min_vol,max_vol"
            .to_string(),
    ];
    let mut evt_rows = vec![
        "seed,grow_rate,death_rate,birth_rate,kind,mcs,sigma,area_at_event,birth_mcs,lifetime_mcs"
            .to_string(),
    ];
    let mut n_matched = 0usize;

    for dir in &subdirs {
        let dirname = dir.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let Some((seed, grow, death, birth)) = parse_params(dirname) else { continue };
        n_matched += 1;
        let prefix = format!("{},{},{},{}", seed, grow, death, birth);

        let vol_path = dir.join("volume_distribution.csv");
        match fs::read_to_string(&vol_path) {
            Ok(content) => {
                for (i, line) in content.lines().enumerate() {
                    if i == 0 || line.trim().is_empty() { continue; }
                    vol_rows.push(format!("{},{}", prefix, line));
                }
            }
            Err(_) => eprintln!("  [skip] no volume_distribution.csv in {dirname}  (run analyze first)"),
        }

        let evt_path = dir.join("events.csv");
        match fs::read_to_string(&evt_path) {
            Ok(content) => {
                for (i, line) in content.lines().enumerate() {
                    if i == 0 || line.trim().is_empty() { continue; }
                    evt_rows.push(format!("{},{}", prefix, line));
                }
            }
            Err(_) => eprintln!("  [skip] no events.csv in {dirname}"),
        }
    }

    if n_matched == 0 {
        eprintln!(
            "No directories matching *_s{{seed}}_g{{grow}}_d{{death}}_b{{birth}} found in {:?}",
            base
        );
        std::process::exit(1);
    }

    println!("Matched {} sweep directories.", n_matched);

    let vol_out = base.join("sweep_volumes.csv");
    let evt_out = base.join("sweep_events.csv");
    fs::write(&vol_out, vol_rows.join("\n") + "\n").expect("cannot write sweep_volumes.csv");
    fs::write(&evt_out, evt_rows.join("\n") + "\n").expect("cannot write sweep_events.csv");
    println!("→ {:?}  ({} data rows)", vol_out, vol_rows.len() - 1);
    println!("→ {:?}  ({} events)",    evt_out, evt_rows.len() - 1);
}
