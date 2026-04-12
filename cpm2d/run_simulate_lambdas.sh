## To run the simulation with 4 cells
echo "Run lambda estimations"
rm data/lambdas/frames/*.json
rm data/lambdas/png/*.png
cargo run --release --bin simulate_lambdas --  --voronoi-start  --la 0.1,1.0,0.10  --lp 0.1,1.0,0.1 --li 0.0,1.0,0.1  --steps 5000 --save-every 5000 --out-dir data/lambdas/
