//! 2-D Cellular Potts Model (CPM)
//!
//!
use std::{f64, fs};
use std::path::Path;
use clap::Parser;
use cpm2d::params::Params;
use cpm2d::grid::Cpm2d;
use cpm2d::render::{save_png, console_print};
use std::f64::consts::PI;

// ─── Main ─────────────────────────────────────────────────────────────────────
#[derive(Parser)]
struct Cli {
    #[arg(long, default_value = "60")]
    grid_w : usize,
    #[arg(long, default_value = "60")]
    grid_h : usize,
    #[arg(long, default_value = "1.0")]
    la: f64,   // --la 0.5,2.0,0.5
    #[arg(long, default_value = "0.2")]
    lp: f64,   // --lp 0.05,0.2,0.05
    #[arg(long, default_value = "0.1")]
    li: f64,   // --li 0.0,0.1,0.1
    #[arg(long, default_value = "16")]
    n_cells: usize,
    #[arg(long, default_value = "data/simulate/")]
    out_dir: String,
    #[arg(long, default_value = "1000")]
    steps: usize,
    #[arg(long, default_value = "100")]
    save_every: usize,
    #[arg(long, default_value = "false")]
    voronoi_start : bool,
    #[arg(long, default_value = "")]
    load : String,
    #[arg(long, default_value = "0")]
    wall_inset: usize,
    #[arg(long, default_value = "200")]
    target_area : i64,
}

fn main() {

    let cli = Cli::parse();
    let path = Path::new(&cli.load);
    let voronoi= cli.voronoi_start;
    let la = cli.la;
    let lp = cli.lp;
    let li = cli.li;
    let n_cells = cli.n_cells;
    let steps = cli.steps;
    let save_every = cli.save_every;
    let wx = cli.grid_w;
    let hx = cli.grid_h;
    let out_dir = cli.out_dir;
    let wall_inset = cli.wall_inset;
    let target_area = cli.target_area;
    let target_perimeter = (target_area as f64 * 4.0 * PI).sqrt() as i64;
    fs::create_dir_all(out_dir.clone()).unwrap();
    fs::create_dir_all(format!("{}/frames/",out_dir.clone())).unwrap();
    fs::create_dir_all(format!("{}/png/",out_dir.clone())).unwrap();

    let mut sim: Cpm2d = if path.exists() {
        Cpm2d::load_state(path.to_str().unwrap())
    }
    else {
        let mut p = Params::default();
        p.lambda_area = la;
        p.lambda_perim = lp;
        p.lambda_iso = li;
        p.n_cells = n_cells;
        p.total_steps = steps;
        p.grid_h = hx;
        p.grid_w = wx;
        p.frames_dir = out_dir.clone();
        p.save_every = save_every;
        p.console_every = save_every;
        p.png_every = save_every;
        p.wall_inset = wall_inset;
        p.target_area = target_area;
        p.target_perim = target_perimeter;
        Cpm2d::new(p, voronoi)
    };

    let p = sim.p.clone();
    let end_mcs = p.total_steps;

    println!(
        "Starting CPM2D  grid={}×{}  cells={}  T={}  MCS={} save_every={}",
        p.grid_w, p.grid_h, p.n_cells, p.temperature, p.total_steps, p.save_every
    );

    console_print(&sim);
    save_png(&sim, None);

    while sim.mcs < end_mcs {
        sim.run_mcs();

        let should_save = cli.save_every > 0 && sim.mcs % cli.save_every == 0;
        let is_last = sim.mcs + 1 == cli.steps;

        if should_save || is_last {
            console_print(&sim);
            save_png(&sim, None);
            sim.save_state(None);
        }
    }

    println!("\nDone. Final MCS={}", sim.mcs);

}