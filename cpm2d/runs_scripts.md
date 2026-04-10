## To run the simulation with 4 cells
cargo run --release --bin simulate_4cells --  --la 1.0,1.0,1.0  --lp 0.1,0.9,0.1 --li 0.0,0.0,0.1  --steps 2000 --save-every 1000 --out-dir sweep4

## To run analysis on 4 cells
cargo run --release --bin analyze_4cells -- sweep4/states_mcs002000_1.00*.json --out results/sweep4.csv


