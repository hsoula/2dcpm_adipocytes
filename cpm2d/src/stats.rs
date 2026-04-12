use std::collections::HashSet;
use crate::{grid::MOORE, grid::Cpm2d};
use crate::boundary::build_boundary;

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
        let a = sim.cells[k as usize].area;
       CellStats {
            mcs:       sim.mcs,
            cell_id:   k,
            area:      a,
            perimeter: sim.cells[k as usize].perimeter,
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

///Now the code. The rotation step is the key — for each boundary pixel,
/// find the first Von Neumann neighbour that differs,
/// compute the rotation angle to bring it to north, then apply that rotation to the full Moore
/// neighbourhood before building the 7-bit key.

// src/stats.rs

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryClass {
    FlatEdge,        // straight boundary, fully convex locally
    OuterCorner,     // convex corner (cell protrudes)
    InnerCorner,     // concave corner (bay opening)
    FilamentTip,     // 1-pixel-wide tip — pathological
    FilamentNeck,    // narrow bridge between two cell regions
    ConcaveBay,      // deep concavity
    Isolated,        // pixel surrounded on all sides by others (shouldn't happen)
    Other(u8),       // catch-all with raw key for inspection
}

/// Moore neighbourhood positions in order, starting from NW going clockwise:
/// NW=0 N=1 NE=2 E=3 SE=4 S=5 SW=6 W=7
const MOORE_CW: [(i32, i32); 8] = [
    (-1, -1), (-1, 0), (-1, 1),
    ( 0,  1),
    ( 1,  1), ( 1, 0), ( 1, -1),
    ( 0, -1),
];

/// Bit assignment after rotation (N is always 0, centre always 1):
/// bit 6=NW  5=NE  4=W  3=E  2=SW  1=S  0=SE
/// Positions in MOORE_CW after rotation to north:
/// index 0=NW 1=N(fixed 0) 2=NE 3=E 4=SE 5=S 6=SW 7=W
const BIT_FOR_POS: [u8; 8] = [6, 255, 5, 3, 0, 1, 2, 4];
// 255 = N slot, always 0, not stored

/// Rotate the 8-neighbour offsets so that `north_idx` becomes index 1 (north).
/// Returns the rotated array of (dr, dc).
fn rotate_to_north(north_idx: usize) -> [(i32, i32); 8] {
    let mut out = [(0i32, 0i32); 8];
    for i in 0..8 {
        out[i] = MOORE_CW[(north_idx + i) % 8];
    }
    out
}

/// For a boundary pixel (r, c) of cell `sigma`, compute the oriented 7-bit key.
/// Returns None if no Von Neumann neighbour differs (pixel is not really boundary).
fn oriented_key(
    grid: &[u32], w: usize, h: usize,
    r: usize, c: usize, sigma: u32,
) -> Option<u8> {
    // Find first Von Neumann neighbour that differs — this becomes "north"
    // Von Neumann positions in MOORE_CW: N=1, E=3, S=5, W=7
    let vn_indices = [1usize, 3, 5, 7];
    let north_idx = vn_indices.iter().find(|&&vi| {
        let (dr, dc) = MOORE_CW[vi];
        let nr = r as i32 + dr;
        let nc = c as i32 + dc;
        if nr < 0 || nr >= h as i32 || nc < 0 || nc >= w as i32 { return true; } // wall = different
        grid[nr as usize * w + nc as usize] != sigma
    })?;

    let rotated = rotate_to_north(*north_idx);
    let mut key = 0u8;

    for (pos, &(dr, dc)) in rotated.iter().enumerate() {
        let bit_slot = BIT_FOR_POS[pos];
        if bit_slot == 255 { continue; } // N slot, always 0
        let nr = r as i32 + dr;
        let nc = c as i32 + dc;
        let is_same = if nr < 0 || nr >= h as i32 || nc < 0 || nc >= w as i32 {
            false // out of bounds = different
        } else {
            grid[nr as usize * w + nc as usize] == sigma
        };
        if is_same {
            key |= 1 << bit_slot;
        }
    }

    Some(key)
}

/// Classify a 7-bit oriented pattern key into a BoundaryClass.
/// key bits: NW=6 NE=5 W=4 E=3 SW=2 S=1 SE=0   (N always 0, centre always 1)
pub fn classify_pattern(key: u8) -> BoundaryClass {
    let nw = (key >> 6) & 1;
    let ne = (key >> 5) & 1;
    let w  = (key >> 4) & 1;
    let e  = (key >> 3) & 1;
    let sw = (key >> 2) & 1;
    let s  = (key >> 1) & 1;
    let se = (key >> 0) & 1;

    let same_count = nw + ne + w + e + sw + s + se;

    match (w, e, s) {
        // Flat edge: both sides full, south full
        (1, 1, 1) => BoundaryClass::FlatEdge,

        // Outer (convex) corner: only south neighbour
        (0, 0, 1) => BoundaryClass::OuterCorner,

        // Filament tip: isolated from all sides
        _ if same_count == 0 => BoundaryClass::FilamentTip,
        _ if same_count == 1 && s == 1 => BoundaryClass::FilamentTip,

        // Concave bay: many same-cell neighbours but opens toward N
        _ if same_count >= 5 => BoundaryClass::ConcaveBay,

        // Inner corner: one lateral side missing
        (1, 0, 1) | (0, 1, 1) => BoundaryClass::InnerCorner,
        (1, 1, 0) => BoundaryClass::InnerCorner,

        // Narrow neck: connected only through diagonals
        _ if w == 0 && e == 0 && s == 0 && same_count > 0 => BoundaryClass::FilamentNeck,

        _ => BoundaryClass::Other(key),
    }
}

/// Per-cell distribution of boundary pattern classes.
#[derive(Debug, Clone, serde::Serialize)]
pub struct BoundaryProfile {
    pub cell_id:        u32,
    pub flat_edge:      u32,
    pub outer_corner:   u32,
    pub inner_corner:   u32,
    pub concave_bay:    u32,
    pub filament_tip:   u32,
    pub filament_neck:  u32,
    pub other:          u32,
    pub total_boundary: u32,
    /// Convexity score: (flat_edge + outer_corner) / total  → 1.0 = perfectly convex
    pub convexity:      f64,
}

pub fn boundary_profiles(sim: &Cpm2d) -> Vec<BoundaryProfile> {
    let W = sim.p.grid_w;
    let H = sim.p.grid_h;
    let n = sim.p.n_cells;

    let mut profiles: Vec<BoundaryProfile> = (0..=n).map(|k| BoundaryProfile {
        cell_id: k as u32,
        flat_edge: 0, outer_corner: 0, inner_corner: 0,
        concave_bay: 0, filament_tip: 0, filament_neck: 0, other: 0,
        total_boundary: 0, convexity: 0.0,
    }).collect();

    for r in 0..H {
        for c in 0..W {
            let s = sim.boundary[r * W + c] as usize;
            if s == 0 { continue; }

            let p = &mut profiles[s];
            p.total_boundary += 1;

            if let Some(key) = oriented_key(&sim.grid, W, H, r, c, s as u32) {
                match classify_pattern(key) {
                    BoundaryClass::FlatEdge     => p.flat_edge     += 1,
                    BoundaryClass::OuterCorner  => p.outer_corner  += 1,
                    BoundaryClass::InnerCorner  => p.inner_corner  += 1,
                    BoundaryClass::ConcaveBay   => p.concave_bay   += 1,
                    BoundaryClass::FilamentTip  => p.filament_tip  += 1,
                    BoundaryClass::FilamentNeck => p.filament_neck += 1,
                    BoundaryClass::Isolated     => p.other         += 1,
                    BoundaryClass::Other(_)     => p.other         += 1,
                }
            }
        }
    }

    // Compute convexity score
    for p in profiles.iter_mut().skip(1) {
        if p.total_boundary > 0 {
            p.convexity = (p.flat_edge + p.outer_corner) as f64
                / p.total_boundary as f64;
        }
    }

    profiles
}

// src/stats.rs

/// Average Moore-neighbour count for boundary pixels of each cell.
/// Returns a Vec indexed by cell id (index 0 unused).
/// A value close to 8 means a very convex compact cell;
/// lower values indicate concave/irregular boundaries.
pub fn average_boundary_convexity(sim: &Cpm2d, boundary: &Vec<u32>) -> Vec<f64> {
    let W = sim.p.grid_w;
    let H = sim.p.grid_h;
    let n = sim.p.n_cells;

    let mut sum   = vec![0u32; n + 1];
    let mut count = vec![0u32; n + 1];

    for r in 0..H {
        for c in 0..W {
            let s = boundary[r * W + c] as usize;
            if s == 0 { continue; }  // not a boundary pixel

            // Count Moore neighbours with the same sigma
            let same = MOORE.iter().filter(|(dr, dc)| {
                let nr = r as i32 + dr;
                let nc = c as i32 + dc;
                nr >= 0 && nr < H as i32 && nc >= 0 && nc < W as i32
                    && sim.grid[nr as usize * W + nc as usize] == s as u32
            }).count() as u32;

            sum[s]   += same;
            count[s] += 1;
        }
    }

    (0..=n).map(|k| {
        if count[k] == 0 { 0.0 }
        else { sum[k] as f64 / count[k] as f64 }
    }).collect()
}

/// Convex hull area for each cell using the shoelace formula
/// on the Andrew's monotone chain hull of boundary pixels.
/// Returns a Vec indexed by cell id (index 0 unused).
pub fn convex_hull_areas(sim: &Cpm2d, boundary: &Vec<u32>) -> Vec<f64> {
    let W = sim.p.grid_w;
    let H = sim.p.grid_h;
    let n = sim.p.n_cells;

    // Collect boundary pixel coordinates per cell
    let mut points: Vec<Vec<(i64, i64)>> = vec![Vec::new(); n + 1];
    for r in 0..H {
        for c in 0..W {
            let s = boundary[r * W + c] as usize;
            if s > 0 {
                points[s].push((c as i64, r as i64));
            }
        }
    }

    (0..=n).map(|k| convex_hull_area(&mut points[k])).collect()
}

/// Andrew's monotone chain — computes convex hull and returns its area
/// via the shoelace formula. Modifies the input (sorts it).
fn convex_hull_area(pts: &mut Vec<(i64, i64)>) -> f64 {
    let n = pts.len();
    if n < 3 { return 0.0; }

    pts.sort_unstable();
    pts.dedup();
    if pts.len() < 3 { return 0.0; }

    let cross = |o: (i64,i64), a: (i64,i64), b: (i64,i64)| -> i64 {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };

    let mut hull: Vec<(i64, i64)> = Vec::with_capacity(pts.len() * 2);

    // Lower hull
    for &p in pts.iter() {
        while hull.len() >= 2 && cross(hull[hull.len()-2], hull[hull.len()-1], p) <= 0 {
            hull.pop();
        }
        hull.push(p);
    }

    // Upper hull
    let lower_len = hull.len() + 1;
    for &p in pts.iter().rev() {
        while hull.len() >= lower_len && cross(hull[hull.len()-2], hull[hull.len()-1], p) <= 0 {
            hull.pop();
        }
        hull.push(p);
    }

    hull.pop(); // last point == first point

    // Shoelace formula
    let area2: i64 = hull.windows(2)
        .map(|w| w[0].0 * w[1].1 - w[1].0 * w[0].1)
        .sum::<i64>()
        .abs();
    // close the polygon
    let last  = *hull.last().unwrap();
    let first = hull[0];
    let close = (last.0 * first.1 - first.0 * last.1).abs();

    (area2 + close) as f64 / 2.0
}