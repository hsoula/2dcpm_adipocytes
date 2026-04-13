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
    /// False = free slot (dead or never used), can be reused for birth
    pub alive: bool,
    /// True = shrinking toward removal (target_area set to 0)
    pub dying: bool,
}

impl CellState {
    pub fn new(id: u32, area: i64, perimeter:i64, target_area: i64, target_perimeter: i64 ) -> Self {
        Self { id, area, perimeter, target_area, target_perimeter, alive: true, dying: false }
    }
}