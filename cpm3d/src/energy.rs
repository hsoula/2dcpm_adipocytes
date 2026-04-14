use std::f64::consts::PI;

// ── Adhesion ──────────────────────────────────────────────────────────────────

#[inline]
pub fn j(s1: u32, s2: u32, j_cm: f64, j_cc: f64) -> f64 {
    if s1 == s2 { return 0.0; }
    if s1 == 0 || s2 == 0 { j_cm } else { j_cc }
}

// ── Volume (same as area in 2D) ───────────────────────────────────────────────

/// ΔH when s_old loses one voxel.
#[inline]
pub fn delta_volume_loss(vol: i64, target: i64, lambda: f64) -> f64 {
    lambda * ((vol - 1 - target).pow(2) - (vol - target).pow(2)) as f64
}

/// ΔH when s_new gains one voxel.
#[inline]
pub fn delta_volume_gain(vol: i64, target: i64, lambda: f64) -> f64 {
    lambda * ((vol + 1 - target).pow(2) - (vol - target).pow(2)) as f64
}

// ── Surface (same structure as perimeter in 2D) ───────────────────────────────

/// ΔH for a change of `ds` in the surface pixel count.
#[inline]
pub fn delta_surface(surf: i64, target: i64, ds: i64, lambda: f64) -> f64 {
    lambda * ((surf + ds - target).pow(2) - (surf - target).pow(2)) as f64
}

// ── Sphericity ────────────────────────────────────────────────────────────────
//
// Sphericity Ψ = (36π)^{1/3} · V^{2/3} / S  ∈ (0, 1].
// Energy: λ · (1 − Ψ)²  →  penalises non-spherical shapes.

pub fn sphericity_energy(vol: i64, surf: i64, lambda: f64) -> f64 {
    if vol < 1 || surf < 1 { return 0.0; }
    let psi = (36.0 * PI).powf(1.0 / 3.0)
        * (vol as f64).powf(2.0 / 3.0)
        / surf as f64;
    lambda * (1.0 - psi).powi(2)
}

pub fn delta_sphericity(
    old_vol:  i64, old_surf: i64,
    new_vol:  i64, new_surf: i64,
    lambda: f64,
) -> f64 {
    sphericity_energy(new_vol, new_surf, lambda)
        - sphericity_energy(old_vol, old_surf, lambda)
}
