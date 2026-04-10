#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CellState {
    pub id:           u32,
    /// Current pixel count (replaces sim.area[k])
    pub area:         i64,
    /// Current 4-connected perimeter (replaces sim.perim[k])
    pub perimeter:    i64,
    /// Target area at the current MCS
    pub target_area:  i64,
    /// Target perimeter at the current MCS
    pub target_perimeter: i64,
}

impl CellState {
    pub fn new(id: u32, target_area: i64, target_perimeter: i64 ) -> Self {
        Self { id, area: 0, perimeter: 0, target_area, target_perimeter }
    }
}