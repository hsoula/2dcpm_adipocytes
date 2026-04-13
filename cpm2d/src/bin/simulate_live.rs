//! CPM simulation with cell demography: linear growth, random death, and birth.
//!
//! Usage example:
//!   cargo run --bin simulate_live -- \
//!     --grid-w 120 --grid-h 120 --n-cells 16 \
//!     --growth-rate 1 --death-prob 0.002 --birth-prob 0.05 \
//!     --steps 5000 --out-dir data/live/

use std::{f64, fs};
use std::path::Path;
use clap::Parser;
use std::f64::consts::PI;
use cpm2d::params::Params;
use cpm2d::grid::Cpm2d;
use cpm2d::render::{save_png, console_print};

#[derive(Parser)]
struct Cli {
    #[arg(long, default_value = "120")]
    grid_w: usize,
    #[arg(long, default_value = "120")]
    grid_h: usize,
    #[arg(long, default_value = "16")]
    n_cells: usize,
    #[arg(long, default_value = "1.0")]
    la: f64,
    #[arg(long, default_value = "0.2")]
    lp: f64,
    #[arg(long, default_value = "0.1")]
    li: f64,
    #[arg(long, default_value = "200")]
    target_area: i64,

    // Demography
    /// Target area added per MCS to each living cell (pixels/step)
    #[arg(long, default_value = "1")]
    growth_rate: f64,
    /// Per-cell probability of dying each MCS
    #[arg(long, default_value = "0.002")]
    death_prob: f64,
    /// Probability of one birth attempt per MCS
    #[arg(long, default_value = "0.05")]
    birth_prob: f64,
    /// Remove dying cell when its area falls below this threshold
    #[arg(long, default_value = "4")]
    death_area_threshold: i64,

    #[arg(long, default_value = "5000")]
    steps: usize,
    /// How often to save JSON state (0 = never)
    #[arg(long, default_value = "500")]
    save_every: usize,
    /// How often to write a PNG frame (0 = never); defaults to save_every if not set
    #[arg(long)]
    png_every: Option<usize>,
    /// How often to print a console summary (0 = never); defaults to save_every if not set
    #[arg(long)]
    console_every: Option<usize>,
    #[arg(long, default_value = "data/live/")]
    out_dir: String,
    #[arg(long, default_value = "0")]
    wall_inset: usize,
    #[arg(long, default_value = "false")]
    voronoi_start: bool,
    #[arg(long, default_value = "")]
    load: String,
}

fn main() {
    let cli = Cli::parse();
    let out_dir = cli.out_dir.clone();

    fs::create_dir_all(&out_dir).unwrap();
    fs::create_dir_all(format!("{}/frames/", out_dir)).unwrap();
    fs::create_dir_all(format!("{}/png/", out_dir)).unwrap();

    let target_perimeter = (cli.target_area as f64 * 4.0 * PI).sqrt() as i64;

    let mut sim: Cpm2d = if !cli.load.is_empty() && Path::new(&cli.load).exists() {
        Cpm2d::load_state(&cli.load)
    } else {
        let mut p = Params::default();
        p.grid_w           = cli.grid_w;
        p.grid_h           = cli.grid_h;
        p.n_cells          = cli.n_cells;
        p.lambda_area      = cli.la;
        p.lambda_perim     = cli.lp;
        p.lambda_iso       = cli.li;
        p.target_area      = cli.target_area;
        p.target_perim     = target_perimeter;
        p.wall_inset       = cli.wall_inset;
        p.total_steps      = cli.steps;
        p.save_every       = cli.save_every;
        p.console_every    = cli.save_every;
        p.png_every        = cli.save_every;
        p.frames_dir       = out_dir.clone();
        p.growth_rate      = cli.growth_rate;
        p.death_prob       = cli.death_prob;
        p.birth_prob       = cli.birth_prob;
        p.death_area_threshold = cli.death_area_threshold;
        Cpm2d::new(p, cli.voronoi_start)
    };

    let png_every     = cli.png_every.unwrap_or(cli.save_every);
    let console_every = cli.console_every.unwrap_or(cli.save_every);

    println!(
        "simulate_live  grid={}×{}  cells={}  T={:.2}  steps={}\n\
         growth_rate={:.2}  death_prob={:.4}  birth_prob={:.4}  death_threshold={}\n\
         save_every={}  png_every={}  console_every={}",
        sim.p.grid_w, sim.p.grid_h, sim.p.n_cells, sim.p.temperature, cli.steps,
        sim.p.growth_rate, sim.p.death_prob, sim.p.birth_prob, sim.p.death_area_threshold,
        cli.save_every, png_every, console_every,
    );

    // MCS 0 snapshot
    console_print(&sim);
    save_png(&sim, None);

    while sim.mcs < cli.steps {
        sim.run_mcs();
        sim.step_demography();

        let is_last = sim.mcs + 1 == cli.steps;

        if console_every > 0 && (sim.mcs % console_every == 0 || is_last) {
            println!(
                "MCS {:5}  living={:3}  dying={:2}",
                sim.mcs, sim.n_living(), sim.n_dying()
            );
            console_print(&sim);
        }

        if png_every > 0 && (sim.mcs % png_every == 0 || is_last) {
            save_png(&sim, None);
        }

        if cli.save_every > 0 && (sim.mcs % cli.save_every == 0 || is_last) {
            sim.save_state(None);
        }
    }

    println!("\nDone. Final MCS={}", sim.mcs);
}
