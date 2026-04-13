## To run the simulation with 4 cells
echo "Run lambda estimations"
rm data/lambdas/frames/*.json
rm data/lambdas/png/*.png
cargo run --release --bin sweep -- --n-cells 1 --la 1.0,1.0,1.0 --lp 0.1,0.9,0.1 --li 0.0,0.0,0.1 --steps 2000 --save-every 1000 --out-dir data/lambdas
