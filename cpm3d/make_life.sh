cargo run --release --bin simulate_life \
-- --grid-w 60 --grid-h 60 --grid-d 60 \
--target-volume 160 --volume-sigma 100  \
--n-cells 2450 --lv 2.0 --ls 0.01 --li 0.01  \
--steps 1000 --out-dir data/sim3d/growth --save-every 100 \
--death-rate 0.0 --birth-rate 0.0 --grow-rate 0.001 --seed 1917

cargo run --bin analyze -- --dir data/sim3d/growth
#python analysis/check_volumes.py
#python analysis/ripley.py
