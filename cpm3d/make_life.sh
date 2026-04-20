cargo run --release --bin simulate_life \
-- --grid-w 60 --grid-h 60 --grid-d 60 \
--target-volume 160 --volume-sigma 100  \
--n-cells 2450 --lv 2.0 --ls 0.01 --li 0.01  \
--steps 1000 --out-dir data/sim3d --save-every 100 \
--death-rate 0.02 --birth-rate 0.1 --seed 42

cargo run --bin analyze -- --dir data/sim3d
python analysis/check_volumes.py
python analysis/ripley.py
