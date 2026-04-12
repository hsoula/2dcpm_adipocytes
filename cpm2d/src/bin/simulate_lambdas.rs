//

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
        while x <= self.max + 1e-9 {   // epsilon avoids float rounding cutting last value
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
    grid_w : usize,
    #[arg(long, default_value = "60")]
    grid_h : usize,
    #[arg(long)]
    la: LambdaRange,   // --la 0.5,2.0,0.5
    #[arg(long)]
    lp: LambdaRange,   // --lp 0.05,0.2,0.05
    #[arg(long)]
    li: LambdaRange,   // --li 0.0,0.1,0.1
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
    voronoi_start : bool,
    #[arg(long, default_value = "")]
    load : String,
    #[arg(long, default_value = "0")]
    wall_inset: usize,
    #[arg(long, default_value = "200")]
    target_area : i64,
}
fn lambda_tag(la: f64, lp: f64, li: f64) -> String {
    // e.g. "1.00_0.10_0.05"  — fixed 2 decimal places, unambiguous in filenames
    format!("{:.2}_{:.2}_{:.2}", la, lp, li)
}

fn name_tag(voronoi_start: bool, target_area:i64, nb_cells: usize) -> String {
    let s = if voronoi_start {"1"} else {"0"};
    format!("{}_{}_{}",s, target_area, nb_cells)
}
fn main() {

    let cli = Cli::parse();
    let la_vals = cli.la.values();
    let lp_vals = cli.lp.values();
    let li_vals = cli.li.values();
    let n_runs = cli.n_runs;
    let n_cells = cli.n_cells;
    let steps = cli.steps;
    let save_every = cli.save_every;
    let wx = cli.grid_w;
    let hx = cli.grid_h;
    let out_dir = cli.out_dir;
    let wall_inset = cli.wall_inset;
    let target_area = cli.target_area;
    let target_perimeter = (target_area as f64 * 4.0 * PI).sqrt() as i64;
    let voronoi= cli.voronoi_start;
    fs::create_dir_all(out_dir.clone()).unwrap();
    fs::create_dir_all(format!("{}/frames/",out_dir.clone())).unwrap();
    fs::create_dir_all(format!("{}/png/",out_dir.clone())).unwrap();

    let mut combos: Vec<(f64, f64, f64)> = Vec::new();
    for la in la_vals.iter().copied() {
        for lp in lp_vals.iter().copied() {
            for li in li_vals.iter().copied() {
                combos.push((la, lp, li));
            }
        }
    }
    println!("{} combinations × {} steps", combos.len(), cli.steps);

    for (la, lp, li) in &combos {
        let tag = lambda_tag(*la, *lp, *li);

        for run in 0..n_runs {

            println!("  running λ={} run {}", tag, run);
            let mut p = Params::default();
            p.grid_w = wx;
            p.grid_h = hx;
            p.lambda_area = *la;
            p.lambda_perim = *lp;
            p.lambda_iso = *li;
            p.total_steps = steps;
            p.n_cells = n_cells;
            p.wall_inset = wall_inset;
            p.target_area = target_area;
            p.target_perim = target_perimeter;
            p.png_every = usize::MAX;   // suppress PNG during sweep
            p.console_every = usize::MAX;
            p.save_every = usize::MAX;   // we control saving ourselves

            let mut sim = Cpm2d::new(p, voronoi);
            let name = name_tag(voronoi, target_area, n_cells);

            for mcs in 0..cli.steps {

                sim.run_mcs();

                let should_save = save_every > 0 && (mcs + 1) % save_every == 0;
                let is_last = mcs + 1 == steps;
                if should_save || is_last {
                    // pattern: states_mcs000300_1.00_0.10_0.05.json
                    let file_name = format!(
                        "{}/frames/states_mcs{:06}_{}_{}_{}.json",
                        out_dir,
                        sim.mcs,
                        tag,
                        name,
                        run
                    );
                    let png_name = format!(
                        "{}/png/states_mcs{:06}_{}_{}_{}.png",
                        out_dir,
                        sim.mcs,
                        tag,
                        name,
                        run
                    );

                    sim.save_state(Some(&file_name));
                    save_png(&sim, Some(&png_name));
                }
            }
        }
    }

    println!("Sweep done → {}", out_dir);
}