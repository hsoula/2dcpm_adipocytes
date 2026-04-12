/// A closed ring of frozen pixels at a given inset from the grid edge.
/// sigma=0 pixels inside this ring are the "wall" — never flipped.
/// The ring can be expanded by pushing it outward (increasing the grid or shrinking inset).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Wall {
    /// Pixel coordinates of the wall ring, stored as a flat HashSet for O(1) lookup
    pub pixels: std::collections::HashSet<(usize, usize)>,
    /// Current inset from the grid edge (1 = outermost ring)
    pub inset: usize,
    pub grid_w: usize,
    pub grid_h: usize,
}

impl Wall {
    /// Build a rectangular ring at `inset` pixels from the grid boundary.
    pub fn new(grid_w: usize, grid_h: usize, inset: usize) -> Self {
        let mut pixels = std::collections::HashSet::new();
        let i = inset;
        // top and bottom rows
        for c in i..grid_w.saturating_sub(i) {
            pixels.insert((i,              c));
            pixels.insert((grid_h - 1 - i, c));
        }
        // left and right columns (excluding corners already inserted)
        for r in (i + 1)..grid_h.saturating_sub(i + 1) {
            pixels.insert((r, i));
            pixels.insert((r, grid_w - 1 - i));
        }
        Self { pixels, inset, grid_w, grid_h }
    }

    #[inline]
    pub fn contains(&self, r: usize, c: usize) -> bool {
        self.pixels.contains(&(r, c))
    }

    /// Expand by moving the wall ring one pixel outward (decrease inset).
    /// Call this when the tissue grows and needs more room.
    /// Returns the newly freed pixels (former wall, now available for cells).
    pub fn expand(&mut self) -> Vec<(usize, usize)> {
        if self.inset == 0 { return vec![]; }
        let old_pixels = self.pixels.clone();
        self.inset -= 1;
        *self = Wall::new(self.grid_w, self.grid_h, self.inset);
        // freed = pixels that were wall but are no longer
        old_pixels.difference(&self.pixels).copied().collect()
    }
}