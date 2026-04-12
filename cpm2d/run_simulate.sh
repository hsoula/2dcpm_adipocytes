echo "Run Main Simulation"
rm data/simulate/frames/*.json
rm data/simulate/png/*.png
cargo run --release --bin simulate --  --voronoi-start  --la 1.0 --lp 0.3 --li 0.0 --steps 10000 --save-every 1000 --out-dir data/simulate