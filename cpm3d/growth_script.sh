dir=data/growth_functional
rm -f $dir/*.*
cargo run --release --bin simulate_life \
-- --grid-w 20 --grid-h 20 --grid-d 20 \
--target-volume 200 --volume-sigma 0  \
--n-cells 40  --lv 2.0 --ls 2.0 --li 0.0  \
--steps 100 --out-dir $dir --save-every 1 \
--death-rate 0.0 --birth-rate 0.0 --grow-rate 0.0 --seed 1917

cargo run --bin analyze -- --dir $dir
python analysis/plot_cell_growth.py $dir --mcs-range 0 100
#python analysis/ripley.py
