## Single interactive run
cargo run --release --bin simulate -- --n-cells 16 --steps 1000 --save-every 100 --out-dir data/simulate

## Parameter sweep (replaces simulate_1cell / simulate_4cells / simulate_lambdas)
# 1-cell sweep
cargo run --release --bin sweep -- --n-cells 1 --la 1.0,1.0,1.0 --lp 0.1,0.9,0.1 --li 0.0,0.0,0.1 --steps 2000 --save-every 1000 --out-dir sweep1

# 4-cell sweep
cargo run --release --bin sweep -- --n-cells 4 --grid-w 40 --grid-h 40 --la 1.0,1.0,1.0 --lp 0.1,0.9,0.1 --li 0.0,0.0,0.1 --steps 2000 --save-every 1000 --out-dir sweep4

# General sweep (e.g. 10 cells, 10 runs)
cargo run --release --bin sweep -- --n-cells 10 --la 0.5,2.0,0.5 --lp 0.05,0.2,0.05 --li 0.0,0.1,0.1 --steps 2000 --n-runs 10 --out-dir sweep_out

## Sweep with normally-distributed cell sizes (replaces simulate_diff_size)
cargo run --release --bin sweep_rand_size -- --la 0.5,2.0,0.5 --lp 0.05,0.2,0.05 --li 0.0,0.1,0.1 --steps 500 --size-mean 10 --size-std 7 --out-dir sweep_rand
    
## Analysis (replaces analyze_1cell / analyze_4cells / analyze_lambdas)
cargo run --release --bin analyze_sweep -- "sweep_out/frames/states_mcs*.json" --out results/sweep.csv

## Basic stats (no lambda metadata required)
cargo run --release --bin analyse -- "data/simulate/frames/state_mcs*.json" --out results/stats.csv

## Debug one MCS step
cargo run --release --bin debug_mcs -- state_stable.json --out debug_energies.csv
