#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CellState {
    pub id:               u32,
    pub volume:           i64,
    pub surface:          i64,
    pub target_volume:    i64,
    pub target_surface:   i64,
    pub alive:            bool,
    pub dying:            bool,
    pub lipid:            f64,
    #[serde(default)]
    pub birth_mcs:        usize,
}
pub const MIN_VOL: i64 = 4;
impl CellState {
    fn new_base(id: u32, target_volume: i64, target_surface: i64, lipid: f64) -> Self {
        Self {
            id,
            volume: 0,
            surface: 0,
            target_volume,
            target_surface,
            alive: true,
            dying: false,
            lipid,
            birth_mcs: 0,
        }
    }

    pub fn new(id: u32, target_volume: i64, target_surface: i64) -> Self {
        Self::new_base(id, target_volume, target_surface, 0.0)
    }

    /// Start with zero lipid and the minimum viable target volume (4 voxels).
    pub fn new_empty(id: u32) -> Self {
        let ts = ((36.0 * std::f64::consts::PI).powf(1.0 / 3.0)
            * (MIN_VOL as f64).powf(2.0 / 3.0)).round() as i64;
        Self::new_base(id, MIN_VOL, ts, 0.0)
    }
}
