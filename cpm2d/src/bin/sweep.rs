//! Lambda parameter sweep — replaces simulate_1cell / simulate_4cells.
//!
//! cargo run --release --bin sweep -- \
//!     --la 0.5,2.0,0.5 --lp 0.05,0.2,0.05 --li 0.0,0.1,0.1 \
//!     --n-cells 10 --steps 2000 --save-every 1000 --out-dir sweep_out

use std::f64::consts::PI;
use clap::Parser;
use cpm2d::{grid::Cpm2d, params::Params};
use std::fs;
use std::str::FromStr;
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
    #[arg(long, default_value = "60")]
    grid_w: usize,
    #[arg(long, default_value = "60")]
    grid_h: usize,
    #[arg(long)]
    la: LambdaRange,
    #[arg(long)]
    lp: LambdaRange,
    #[arg(long)]
    li: LambdaRange,
    #[arg(long, default_value = "simulate/data")]
    out_dir: String,
    #[arg(long, default_value = "2000")]
    steps: usize,
    #[arg(long, default_value = "1000")]
    save_every: usize,
    #[arg(long, default_value = "10")]
    n_runs: usize,
    #[arg(long, default_value = "10")]
    n_cells: usize,
    #[arg(long, default_value = "false")]
    voronoi_start: bool,
    #[arg(long, default_value = "0")]
    wall_inset: usize,
    #[arg(long, default_value = "200")]
    target_area: i64,
}

fn lambda_tag(la: f64, lp: f64, li: f64) -> String {
    format!("{:.2}_{:.2}_{:.2}", la, lp, li)
}

fn name_tag(voronoi_start: bool, target_area: i64, nb_cells: usize) -> String {
    let s = if voronoi_start { "1" } else { "0" };
    format!("{}_{}_{}", s, target_area, nb_cells)
}

fn main() {
    let cli = Cli::parse();
    let la_vals = cli.la.values();
    let lp_vals = cli.lp.values();
    let li_vals = cli.li.values();
    let target_perimeter = (cli.target_area as f64 * 4.0 * PI).sqrt() as i64;

    fs::create_dir_all(&cli.out_dir).unwrap();
    fs::create_dir_all(format!("{}/frames/", cli.out_dir)).unwrap();
    fs::create_dir_all(format!("{}/png/", cli.out_dir)).unwrap();

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
            p.grid_w = cli.grid_w;
            p.grid_h = cli.grid_h;
            p.lambda_area = *la;
            p.lambda_perim = *lp;
            p.lambda_iso = *li;
            p.total_steps = cli.steps;
            p.n_cells = cli.n_cells;
            p.wall_inset = cli.wall_inset;
            p.target_area = cli.target_area;
            p.target_perim = target_perimeter;
            p.png_every = usize::MAX;
            p.console_every = usize::MAX;
            p.save_every = usize::MAX;

            let mut sim = Cpm2d::new(p, cli.voronoi_start);
            let name = name_tag(cli.voronoi_start, cli.target_area, cli.n_cells);

            for mcs in 0..cli.steps {
                sim.run_mcs();

                let should_save = cli.save_every > 0 && (mcs + 1) % cli.save_every == 0;
                let is_last = mcs + 1 == cli.steps;
                if should_save || is_last {
                    let file_name = format!(
                        "{}/frames/states_mcs{:06}_{}_{}_{}.json",
                        cli.out_dir, sim.mcs, tag, name, run
                    );
                    let png_name = format!(
                        "{}/png/states_mcs{:06}_{}_{}_{}.png",
                        cli.out_dir, sim.mcs, tag, name, run
                    );
                    sim.save_state(Some(&file_name));
                    save_png(&sim, Some(&png_name));
                }
            }
        }
    }

    println!("Sweep done → {}", cli.out_dir);
}
