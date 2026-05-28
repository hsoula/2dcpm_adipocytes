dir=data/growth_functional
rm -f $dir/*.*
cargo run --release --bin simulate_life \
-- --grid-w 20 --grid-h 20 --grid-d 20 \
--target-volume 100 --volume-sigma 0  \
--n-cells 450 --lv 2.0 --ls 0.01 --li 0.01  \
--steps 10 --out-dir $dir --save-every 1 \
--death-rate 0.0 --birth-rate 0.0 --grow-rate 0.0 --seed 1917

cargo run --bin analyze -- --dir $dir
python analysis/plot_cell_growth.py
#python analysis/ripley.py
