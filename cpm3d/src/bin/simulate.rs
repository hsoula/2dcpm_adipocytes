//! 3-D CPM simulation.
//! cargo run --bin simulate -- --n-cells 8 --steps 1000 --out-dir data/sim3d/

use std::fs;
use std::path::Path;
use clap::Parser;
use cpm3d::params::Params;
use cpm3d::grid::Cpm3d;

#[derive(Parser)]
struct Cli {
    #[arg(long, default_value = "60")]  grid_w: usize,
    #[arg(long, default_value = "60")]  grid_h: usize,
    #[arg(long, default_value = "60")]  grid_d: usize,
    #[arg(long, default_value = "8")]   n_cells: usize,
    #[arg(long, default_value = "125")] target_volume: i64,
    #[arg(long, default_value = "1.0")] lv: f64,
    #[arg(long, default_value = "0.5")] ls: f64,
    #[arg(long, default_value = "0.05")] li: f64,
    #[arg(long, default_value = "2000")] steps: usize,
    #[arg(long, default_value = "100")]  save_every: usize,
    #[arg(long, default_value = "50")]   png_every: usize,
    #[arg(long, default_value = "data/sim3d/")] out_dir: String,
    #[arg(long, default_value = "")] load: String,
}

fn main() {
    let cli = Cli::parse();
    fs::create_dir_all(&cli.out_dir).unwrap();

    let mut sim: Cpm3d = if !cli.load.is_empty() && Path::new(&cli.load).exists() {
        Cpm3d::load_state(&cli.load)
    } else {
        let mut p = Params::default();
        p.grid_w        = cli.grid_w;
        p.grid_h        = cli.grid_h;
        p.grid_d        = cli.grid_d;
        p.n_cells       = cli.n_cells;
        p.target_volume = cli.target_volume;
        p.lambda_vol    = cli.lv;
        p.lambda_surf   = cli.ls;
        p.lambda_spher  = cli.li;
        p.total_steps   = cli.steps;
        p.save_every    = cli.save_every;
        p.png_every     = cli.png_every;
        p.out_dir       = cli.out_dir.clone();
        Cpm3d::new(p)
    };

    println!(
        "cpm3d  grid={}×{}×{}  cells={}  T={:.2}  steps={}",
        sim.p.grid_w, sim.p.grid_h, sim.p.grid_d,
        sim.p.n_cells, sim.p.temperature, cli.steps
    );

    sim.print_stats();

    while sim.mcs < cli.steps {
        sim.run_mcs();
        let is_last = sim.mcs + 1 == cli.steps;

        if cli.save_every > 0 && (sim.mcs % cli.save_every == 0 || is_last) {
            sim.print_stats();
            sim.save_state(None);
        }

        if cli.png_every > 0 && (sim.mcs % cli.png_every == 0 || is_last) {
            let (w, h, d) = (sim.p.grid_w, sim.p.grid_h, sim.p.grid_d);
            let tag = |ax: &str| format!("{}/slice_{}_mcs{:06}.ppm",
                sim.p.out_dir.trim_end_matches('/'), ax, sim.mcs);
            sim.save_slice_png(0, d / 2, &tag("xy"));
            sim.save_slice_png(1, h / 2, &tag("xz"));
            sim.save_slice_png(2, w / 2, &tag("yz"));
        }
    }

    println!("Done. MCS={}", sim.mcs);
}
