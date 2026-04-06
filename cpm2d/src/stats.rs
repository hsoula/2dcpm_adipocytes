use std::collections::HashSet;
use crate::{params::Params, grid::Cpm2d};

#[derive(Debug, Clone, serde::Serialize)]
pub struct CellStats {
    pub mcs:       usize,
    pub cell_id:   u32,
    pub area:      i64,
    pub perimeter: i64,
    pub com_x:     f64,
    pub com_y:     f64,
}

#[derive(Debug, serde::Serialize)]
pub struct NeighbourEdge {
    pub mcs:  usize,
    pub cell_a: u32,
    pub cell_b: u32,
}

/// Compute centre of mass for every live cell.
/// Returns one CellStats per cell, sorted by cell_id.
pub fn compute_stats(sim: &Cpm2d) -> (Vec<CellStats>, Vec<NeighbourEdge>) {
    let W = sim.p.grid_w;
    let H = sim.p.grid_h;
    let n = sim.p.n_cells;

    // Accumulators: sum_x, sum_y, count per cell
    let mut sum_x = vec![0i64; n + 1];
    let mut sum_y = vec![0i64; n + 1];

    // Neighbour sets: for each cell, which other cells touch it?
    let mut neighbours: Vec<HashSet<u32>> = vec![HashSet::new(); n + 1];

    for r in 0..H {
        for c in 0..W {
            let s = sim.grid[r * W + c] as usize;
            if s == 0 { continue; }
            sum_x[s] += c as i64;
            sum_y[s] += r as i64;

            // Check Moore neighbourhood for contact with other cells
            for (dr, dc) in crate::grid::MOORE {
                let nr = r as i32 + dr;
                let nc = c as i32 + dc;
                if nr < 0 || nr >= H as i32 || nc < 0 || nc >= W as i32 { continue; }
                let nb = sim.grid[nr as usize * W + nc as usize];
                if nb > 0 && nb != s as u32 {
                    neighbours[s].insert(nb);
                }
            }
        }
    }

    let cell_stats =  (1..=n as u32).map(|k| {
        let a = sim.area[k as usize];
       CellStats {
            mcs:       sim.mcs,
            cell_id:   k,
            area:      a,
            perimeter: sim.perimeter[k as usize],
            com_x:     if a > 0 { sum_x[k as usize] as f64 / a as f64 } else { 0.0 },
            com_y:     if a > 0 { sum_y[k as usize] as f64 / a as f64 } else { 0.0 },
        }
    }).collect();
    let mut edges = Vec::new();
    for (a, set) in neighbours.iter().enumerate().skip(1) {
        for &b in set {
            if (a as u32) < b {
                edges.push(NeighbourEdge { mcs: sim.mcs, cell_a: a as u32, cell_b: b });
            }
        }
    }
    (cell_stats, edges)
}