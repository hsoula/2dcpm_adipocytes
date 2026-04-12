use crate::params::Params;
use rand::prelude::*;
use crate::wall::Wall;
/// Poisson disk sampling via Bridson's algorithm.
/// Returns at most `n` points guaranteed to be >= min_dist apart.

pub fn poisson_seeds(
    w: usize, h: usize,
    min_dist: f64,
    n: usize,
    wall: &Wall,
    rng: &mut impl Rng,
) -> Vec<(f64, f64)> {
    let cell_size = min_dist / std::f64::consts::SQRT_2;
    let grid_w    = (w as f64 / cell_size).ceil() as usize + 1;
    let grid_h    = (h as f64 / cell_size).ceil() as usize + 1;

    let mut bg: Vec<Option<usize>> = vec![None; grid_w * grid_h];
    let mut points: Vec<(f64, f64)> = Vec::new();
    let mut active: Vec<usize>       = Vec::new();

    let bg_idx = |px: f64, py: f64| -> usize {
        let gx = (px / cell_size) as usize;
        let gy = (py / cell_size) as usize;
        gy * grid_w + gx
    };

    // Compute the valid inner bounding box from the wall inset
    let lo = wall.inset as f64 + 1.0;  // just inside the wall ring
    let hi_c = w as f64 - lo;
    let hi_r = h as f64 - lo;

    if lo >= hi_c || lo >= hi_r {
        eprintln!("Warning: wall inset too large, no room for seeds");
        return vec![];
    }

    // First seed: random point strictly inside the wall
    let first = (rng.gen_range(lo..hi_c), rng.gen_range(lo..hi_r));
    bg[bg_idx(first.0, first.1)] = Some(0);
    points.push(first);
    active.push(0);

    const K: usize = 30;

    while !active.is_empty() && points.len() < n {
        let idx      = rng.gen_range(0..active.len());
        let (px, py) = points[active[idx]];
        let mut found = false;

        for _ in 0..K {
            let angle = rng.gen_range(0.0..std::f64::consts::TAU);
            let dist  = rng.gen_range(min_dist..2.0 * min_dist);
            let cx    = px + angle.cos() * dist;
            let cy    = py + angle.sin() * dist;

            // Must stay strictly inside the wall ring
            if cx < lo || cx >= hi_c || cy < lo || cy >= hi_r { continue; }

            // Also check the pixel itself isn't a wall pixel
            if wall.contains(cy as usize, cx as usize) { continue; }

            let gx = (cx / cell_size) as usize;
            let gy = (cy / cell_size) as usize;
            let mut ok = true;

            'outer: for dy in gy.saturating_sub(2)..=(gy + 2).min(grid_h - 1) {
                for dx in gx.saturating_sub(2)..=(gx + 2).min(grid_w - 1) {
                    if let Some(qi) = bg[dy * grid_w + dx] {
                        let (qx, qy) = points[qi];
                        if (cx-qx).powi(2) + (cy-qy).powi(2) < min_dist * min_dist {
                            ok = false;
                            break 'outer;
                        }
                    }
                }
            }

            if ok {
                let new_idx = points.len();
                bg[gy * grid_w + gx] = Some(new_idx);
                points.push((cx, cy));
                active.push(new_idx);
                found = true;
                break;
            }
        }

        if !found { active.swap_remove(idx); }
    }

    points
}
/// Fill the grid with a Voronoi diagram based on `seeds`.
/// The outermost ring of pixels is set to 0 (frozen wall).
/// Seeds are assigned sigma 1..=seeds.len().
pub fn voronoi_fill(
    grid: &mut Vec<u32>,
    w: usize, h: usize,
    seeds: &[(f64, f64)],
    wall: &Wall,
) {
    let lo = wall.inset  + 1;  // just inside the wall ring
    let hi_c = w  - lo;
    let hi_r = h  - lo;

    for r in lo..hi_r {
        for c in lo..hi_c {
            if wall.contains(r, c) {
                grid[r * w + c] = 0;  // frozen wall pixel
                continue;
            }
            // nearest seed — zero medium, complete fill
            let best = seeds.iter().enumerate()
                .min_by(|(_, a), (_, b)| {
                    let da = (c as f64 - a.0).powi(2) + (r as f64 - a.1).powi(2);
                    let db = (c as f64 - b.0).powi(2) + (r as f64 - b.1).powi(2);
                    da.partial_cmp(&db).unwrap()
                });
            grid[r * w + c] = best.map(|(i, _)| i as u32 + 1).unwrap_or(0);
        }
    }
}