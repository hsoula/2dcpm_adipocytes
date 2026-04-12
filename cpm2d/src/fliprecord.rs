use serde::Serialize;

#[derive(Serialize)]
pub struct FlipRecord {
    pub attempt:       usize,
    pub r:             usize,
    pub c:             usize,
    pub s_old:         u32,
    pub s_new:         u32,
    // per-term delta H
    pub dh_contact:    f64,
    pub dh_area_old:   f64,   // contribution from s_old losing a pixel
    pub dh_area_new:   f64,   // contribution from s_new gaining a pixel
    pub dh_perim:      f64,
    pub dh_iso:        f64,
    pub dh_total:      f64,
    // state before flip
    pub area_old:      i64,
    pub target_old:    i64,
    pub area_new:      i64,
    pub target_new:    i64,
    // outcome
    pub accepted:      bool,
    // energy balance flags
    pub area_dominated: bool,   // |dh_area| > |dh_contact + dh_perim|
    pub negative_runaway: bool, // dh_total < -threshold
}