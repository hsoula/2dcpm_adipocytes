//! 3-D CPM simulation.
//! cargo run --bin simulate_life -- --n-cells 8 --steps 1000 --out-dir data/life

use std::fs;
use std::path::Path;
use clap::Parser;
use cpm3d::params::Params;
use cpm3d::grid::Cpm3d;
use std::io::Write as IoWrite;
use cpm3d::dynamics::EventKind;
#[derive(Parser)]
struct Cli {
    #[arg(long, default_value = "60")]  grid_w: usize,
    #[arg(long, default_value = "60")]  grid_h: usize,
    #[arg(long, default_value = "60")]  grid_d: usize,
    #[arg(long, default_value = "8")]   n_cells: usize,
    /// Mean target volume (mu) for all cells.
    #[arg(long, default_value = "125")] target_volume: i64,
    /// Std-dev of Gaussian noise on initial target volume per cell (0 = deterministic).
    #[arg(long, default_value = "0.0")] volume_sigma: f64,
    #[arg(long, default_value = "1.0")] lv: f64,
    #[arg(long, default_value = "0.5")] ls: f64,
    #[arg(long, default_value = "0.05")] li: f64,
    #[arg(long, default_value = "2000")] steps: usize,
    #[arg(long, default_value = "100")]  save_every: usize,
    #[arg(long, default_value = "0")]   png_every: usize,
    #[arg(long, default_value = "data/sim3d/")] out_dir: String,
    #[arg(long, default_value = "")] load: String,
    #[arg(long, default_value = "0.02")] death_rate: f64,
    #[arg(long, default_value = "0.01")] birth_rate: f64,
    #[arg(long, default_value = "1.0")] grow_rate: f64,
    /// RNG seed for reproducibility (omit for random).
    #[arg(long)] seed: Option<u64>,
}

fn main() {
    let cli = Cli::parse();
    let seed_str = cli.seed.map(|s| s.to_string()).unwrap_or_else(|| "rand".to_string());
    let dir = format!("{}_s{}_g{}_d{}_b{}",
        cli.out_dir.trim_end_matches('/'), seed_str,
        cli.grow_rate, cli.death_rate, cli.birth_rate);
    fs::create_dir_all(&dir).unwrap();
    fs::create_dir_all(format!("{}/events/", &dir)).unwrap();

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
        p.volume_sigma  = cli.volume_sigma;
        p.total_steps   = cli.steps;
        p.save_every    = cli.save_every;
        p.png_every     = cli.png_every;
        p.out_dir       = dir.clone();
        p.growth_rate   = cli.grow_rate;
        p.birth_rate    = cli.birth_rate;
        p.death_rate    = cli.death_rate;
        p.seed          = cli.seed;
        Cpm3d::new_empty(p)
    };

    println!(
        "cpm3d  grid={}×{}×{}  cells={}  T={:.2}  steps={}",
        sim.p.grid_w, sim.p.grid_h, sim.p.grid_d,
        sim.p.n_cells, sim.p.temperature, cli.steps
    );

    let events_path = format!("{}/events.csv", &dir);
    let mut events_csv = fs::File::create(&events_path).expect("cannot create events.csv");
    writeln!(events_csv, "kind,mcs,sigma,area_at_event,birth_mcs,lifetime_mcs")
        .expect("cannot write events.csv header");

    sim.print_stats();
    sim.save_state(None);
    while sim.mcs < cli.steps {
        sim.run_mcs();

        // -- do the demography
        let dem_events = sim.step_demography();

        // ── Handle demography events ─────────────────────────────────────────
        if !dem_events.is_empty() {
            for ev in &dem_events {
                let kind_str = match ev.kind {
                    EventKind::Birth => "birth",
                    EventKind::Dead => "death",
                    EventKind::Dying => "dying",
                };
                writeln!(
                    events_csv,
                    "{},{},{},{},{},{}",
                    kind_str, ev.mcs, ev.sigma,
                    ev.volume_at_event, ev.birth_mcs, ev.lifetime_mcs
                ).expect("cannot write event row");
                //
                // // PNG snapshot for every event
                // let tag = format!("{}/{}/event_{}_{}_mcs{:06}_s{}.png",
                //                   &cli.out_dir, "events", kind_str, ev.sigma, ev.mcs, ev.sigma);
                // save_png(&sim, Some(&tag));
            }
        }

        let is_last = sim.mcs + 1 == cli.steps;

        if cli.save_every > 0 && (sim.mcs % cli.save_every == 0 || is_last) {
            sim.print_stats();
            sim.save_state(None);
        }

        if cli.png_every > 0 && (sim.mcs % cli.png_every == 0 || is_last) {
            let (w, h, d) = (sim.p.grid_w, sim.p.grid_h, sim.p.grid_d);
            let tag = |ax: &str| format!("{}/slice_{}_mcs{:06}.png",
                                         sim.p.out_dir.trim_end_matches('/'), ax, sim.mcs);
            sim.save_slice_png(0, d / 2, &tag("xy"));
            sim.save_slice_png(1, h / 2, &tag("xz"));
            sim.save_slice_png(2, w / 2, &tag("yz"));
        }
    }

    println!("Done. MCS={}", sim.mcs);
}
