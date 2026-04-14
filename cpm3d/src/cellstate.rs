#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CellState {
    pub id:               u32,
    pub volume:           i64,
    pub surface:          i64,
    pub target_volume:    i64,
    pub target_surface:   i64,
    pub alive:            bool,
    pub dying:            bool,
    #[serde(default)]
    pub birth_mcs:        usize,
}

impl CellState {
    pub fn new(id: u32, target_volume: i64, target_surface: i64) -> Self {
        Self {
            id,
            volume: 0,
            surface: 0,
            target_volume,
            target_surface,
            alive: true,
            dying: false,
            birth_mcs: 0,
        }
    }
}
