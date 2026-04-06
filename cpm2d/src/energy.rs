use std::f64::consts::PI;

/// Convex hull energy: penalizes area "missing" from the convex hull.
/// hull_area must be recomputed externally (e.g. Andrew's monotone chain).
pub fn energy_convex_hull(hull_area:f64, area:f64, lambda: f64) -> f64 {
    let deficit = (hull_area - area).max(0.0);
    lambda * deficit
}

/// ΔH for a Metropolis copy attempt.
/// new_hull_area: recomputed hull after the proposed spin flip.
pub fn delta_convex_hull(
    old_area:f64,
    old_hull_area:f64,
    new_area: f64,
    new_hull_area: f64,
    lambda: f64,
) -> f64 {
    let h_old = energy_convex_hull(old_hull_area, old_area, lambda);
    energy_convex_hull(new_hull_area, new_area, lambda) - h_old
}

/// Shape index energy. p0 ≈ 3.54 (≈2√π) targets a circular/convex cell.
/// p0 > 2√π drives elongated shapes; p0 < 2√π is unphysical.
pub fn energy_shape_index(area: f64, perimeter:f64, lambda_s: f64, p0: f64) -> f64 {
    let target = p0 * area.sqrt();
    let diff   = perimeter - target;
    lambda_s * diff * diff
}

/// ΔH — O(1), uses only updated P and A.
pub fn delta_shape_index(
    old_area: f64,
    old_perimeter:f64,
    new_area: f64,
    new_perim: f64,
    lambda_s: f64,
    p0: f64,
) -> f64 {
    let h_old = energy_shape_index(old_area, old_perimeter, lambda_s, p0);
    energy_shape_index(new_area, new_perim, lambda_s, p0) - h_old
}
/// Isoperimetric ratio energy. A circle achieves the minimum 4π ≈ 12.57.
/// Any concavity or elongation raises P²/A above this value.
pub fn energy_isoperimetric(area: f64, perimeter:f64, lambda_iso: f64) -> f64 {
    if area < 1e-9 { return 0.0; }
    lambda_iso * perimeter * perimeter / area
}

/// Normalised variant: penalises excess above the circular minimum 4π.
pub fn energy_isoperimetric_normalised(area: f64, perimeter:f64, lambda_iso: f64) -> f64 {
    if area < 1e-9 { return 0.0; }
    let excess = perimeter * perimeter / area - 4.0 * PI;
    lambda_iso * excess.max(0.0)
}

/// ΔH — O(1).
pub fn delta_isoperimetric(
    old_area: f64,
    old_perimeter: f64,
    new_area: f64,
    new_perim: f64,
    lambda_iso: f64,
) -> f64 {
    let h_old = energy_isoperimetric_normalised(old_area,old_perimeter, lambda_iso);
    energy_isoperimetric_normalised(new_area, new_perim, lambda_iso) - h_old
}

/// Lattice curvature at pixel (x,y): κ ≈ n_neighbors - threshold.
/// n < threshold → convex protrusion (κ > 0, no penalty)
/// n > threshold → concave bay      (κ < 0, penalised)
/// threshold = 4 works well for Moore neighbourhood (8-connected).
pub fn local_curvature(neighbor_count: u8, threshold: f64) -> f64 {
    threshold - neighbor_count as f64  // positive = convex, negative = concave
}

/// Curvature energy: integrates concave contributions over the boundary.
// pub fn energy_curvature(cell: &Cell, lambda_curv: f64, threshold: f64) -> f64 {
//     cell.boundary_pixels
//         .iter()
//         .map(|&(_, _, n)| {
//             let kappa = local_curvature(n, threshold);
//             (-kappa).max(0.0)  // only concave pixels contribute
//         })
//         .sum::<f64>()
//         * lambda_curv
// }

/// ΔH for flipping pixel (x,y): recompute only affected boundary pixels.
pub fn delta_curvature(
    old_boundary: &[(i32, i32, u8)],
    new_boundary: &[(i32, i32, u8)],
    lambda_curv: f64,
    threshold: f64,
) -> f64 {
    let sum_concave = |pixels: &[(i32,i32,u8)]| -> f64 {
        pixels.iter().map(|&(_,_,n)| (-local_curvature(n,threshold)).max(0.0)).sum()
    };
    lambda_curv * (sum_concave(new_boundary) - sum_concave(old_boundary))
}
/// Compute total ΔH for a Metropolis copy attempt.
/// Call with pre-computed new_area / new_perim from the standard CPM update.
pub fn delta_total_convexity    (
    old_area: f64,
    old_perimeter: f64,
    new_area: f64,
    new_perim: f64,
    old_boundary: &[(i32, i32, u8)],
    new_boundary: &[(i32, i32, u8)],
    lambda_s: f64,    // shape index weight
    p0: f64,           // target shape index (~3.54)
    lambda_iso: f64,  // isoperimetric weight (optional, set 0 to disable)
    lambda_curv: f64, // curvature weight      (optional, set 0 to disable)
    threshold: f64,   // curvature threshold   (typically 4.0)
) -> f64 {
    delta_shape_index(old_area,old_perimeter , new_area, new_perim, lambda_s, p0)
        + delta_isoperimetric(old_area,old_perimeter , new_area, new_perim, lambda_iso)
        + delta_curvature(old_boundary, new_boundary, lambda_curv, threshold)
}

#[cfg(test)]
mod tests {
    use super::*;   // imports everything from energy.rs

    #[test]
    fn isoperimetric_non_negative() {
        assert!(energy_isoperimetric_normalised(100.0, 40.0, 1.0) >= 0.0);
    }

    #[test]
    fn delta_is_antisymmetric() {
        let fwd = delta_isoperimetric(100.0, 40.0, 101.0, 40.5, 1.0);
        let bwd = delta_isoperimetric(101.0, 40.5, 100.0, 40.0, 1.0);
        assert!((fwd + bwd).abs() < 1e-10);
    }
}