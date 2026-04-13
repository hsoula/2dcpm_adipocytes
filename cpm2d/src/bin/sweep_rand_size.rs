//! Lambda sweep with normally-distributed cell target areas.
//!
//! cargo run --release --bin sweep_rand_size -- \
//!     --la 0.5,2.0,0.5 --lp 0.05,0.2,0.05 --li 0.0,0.1,0.1 \
//!     --steps 500 --out-dir sweep_rand

use clap::Parser;
use cpm2d::{grid::Cpm2d, params::Params};
use std::fs;
use std::str::FromStr;
use rand::prelude::*;
use rand_distr::num_traits::Pow;
use cpm2d::cellstate::CellState;

use cpm2d::render::save_png;

/// Parses "min,max,step" from the command line.
#[derive(Debug, Clone)]
pub struct LambdaRange {
    pub min:  f64,
    pub max:  f64,
    pub step: f64,
}

impl LambdaRange {
    pub fn values(&self) -> Vec<f64> {
        let mut v = Vec::new();
        let mut x = self.min;
        while x <= self.max + 1e-9 {
            v.push(x);
            x += self.step;
        }
        v
    }
}

impl FromStr for LambdaRange {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 3 {
            return Err(format!("expected min,max,step — got '{s}'"));
        }
        let [min, max, step] = [parts[0], parts[1], parts[2]]
            .map(|p| p.parse::<f64>()
                .map_err(|e| format!("cannot parse '{}': {e}", p)));
        Ok(Self { min: min?, max: max?, step: step? })
    }
}

#[derive(Parser)]
struct Cli {
    #[arg(long)]
    la: LambdaRange,
    #[arg(long)]
    lp: LambdaRange,
    #[arg(long)]
    li: LambdaRange,
    #[arg(long, default_value = "sweep_states")]
    out_dir: String,
    #[arg(long, default_value = "300")]
    steps: usize,
    #[arg(long, default_value = "5000")]
    save_every: usize,
    #[arg(long, default_value = "10")]
    n_runs: usize,
    /// Mean cell size (side length)
    #[arg(long, default_value = "10")]
    size_mean: f64,
    /// Std-dev of cell size distribution
    #[arg(long, default_value = "7")]
    size_std: f64,
}

fn lambda_tag(la: f64, lp: f64, li: f64) -> String {
    format!("{:.2}_{:.2}_{:.2}", la, lp, li)
}

fn main() {
    let cli = Cli::parse();
    fs::create_dir_all(&cli.out_dir).unwrap();
    let la_vals = cli.la.values();
    let lp_vals = cli.lp.values();
    let li_vals = cli.li.values();

    let mut combos: Vec<(f64, f64, f64)> = Vec::new();
    for la in la_vals.iter().copied() {
        for lp in lp_vals.iter().copied() {
            for li in li_vals.iter().copied() {
                combos.push((la, lp, li));
            }
        }
    }
    println!("{} combinations × {} steps × {} runs", combos.len(), cli.steps, cli.n_runs);

    for (la, lp, li) in &combos {
        let tag = lambda_tag(*la, *lp, *li);

        for run in 0..cli.n_runs {
            println!("  running λ={} run {}", tag, run);
            let mut p = Params::default();
            p.grid_w = 100;
            p.grid_h = 100;
            p.temperature = 10f64;
            p.lambda_area = *la;
            p.lambda_perim = *lp;
            p.lambda_iso = *li;
            p.total_steps = cli.steps;
            p.n_cells = 64;
            p.png_every = usize::MAX;
            p.console_every = usize::MAX;
            p.save_every = usize::MAX;
            p.wall_inset = 0;

            let mut sim = Cpm2d::new(p, true);
            sim.cells = (0..sim.p.n_cells + 1)
                .map(|k| {
                    let size = sim.rng.gen_range(cli.size_mean - cli.size_std..cli.size_mean + cli.size_std);
                    CellState::new(
                        k as u32,
                        sim.cells[k].area,
                        sim.cells[k].perimeter,
                        size.pow(2) as i64,
                        (size.pow(2) as f64 * 4.0).sqrt() as i64,
                    )
                })
                .collect();

            for mcs in 0..cli.steps {
                sim.run_mcs();

                let should_save = cli.save_every > 0 && (mcs + 1) % cli.save_every == 0;
                let is_last = mcs + 1 == cli.steps;
                if should_save || is_last {
                    let fname = format!(
                        "{}/states_mcs{:06}_{}_{}.json",
                        cli.out_dir, sim.mcs, tag, run
                    );
                    let png_name = format!(
                        "{}/png/states_mcs{:06}_{}_{}.png",
                        cli.out_dir, sim.mcs, tag, run
                    );
                    sim.save_state(Some(&fname));
                    save_png(&sim, Some(&png_name));
                }
            }
        }
    }

    println!("Sweep done → {}", cli.out_dir);
}
