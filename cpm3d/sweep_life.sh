for seed in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
cargo run --release --bin simulate_life \
-- --grid-w 60 --grid-h 60 --grid-d 60 \
--target-volume 160 --volume-sigma 100  \
--n-cells 2450 --lv 2.0 --ls 0.01 --li 0.01  \
--steps 10 --out-dir data/sim3d --save-every 100 \
--death-rate 0.01 --birth-rate 0.2 --grow-rate 0.1 --seed $seed > data/out$seed.txt
done
