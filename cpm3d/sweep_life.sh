# 1. run combinations
for seed in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do
  for grow in 0.001 0.01 0.1; do
    for death in 0.01 0.02; do
      for birth in 0.005 0.01; do
        cargo run --release --bin simulate_life -- \
          --grid-w 60 --grid-h 60 --grid-d 60 \
          --seed $seed --grow-rate $grow --death-rate $death --birth-rate $birth \
          --n-cells 2450 --lv 2.0 --ls 0.01 --li 0.01  \
          --save-every 5000 \
          --out-dir data/sim3d --steps 5000
      done
    done
  done
done

# 2. generate volume CSVs per run
for d in data/sim3d_*; do
  cargo run --release --bin analyze -- --dir $d
done

# 3. aggregate everything
cargo run --release --bin collect_results -- --base-dir data/