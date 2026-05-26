#[derive(Debug, Clone, serde::Serialize)]
pub enum EventKind {
    Birth,
    Dying,
    Dead,
    /// Cell gained one voxel during a copy attempt.
    Grow,
    /// Cell lost one voxel during a copy attempt.
    Shrink,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct DemographyEvent {
    pub kind:            EventKind,
    pub sigma:           u32,
    pub mcs:             usize,
    pub volume_at_event: i64,
    pub birth_mcs:       usize,
    pub lifetime_mcs:    usize,
}
