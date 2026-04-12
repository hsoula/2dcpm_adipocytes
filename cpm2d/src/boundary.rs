const MOORE: [(i32, i32); 8] = [
    (-1, -1), (-1, 0), (-1, 1),
    ( 0, -1),           ( 0, 1),
    ( 1, -1),  ( 1, 0), ( 1, 1),
];

pub fn build_boundary(grid: &[u32], w: usize, h: usize) -> Vec<u32> {
    let mut boundary = vec![0u32; w * h];
    for r in 0..h {
        for c in 0..w {
            let s = grid[r * w + c];
            if s > 0 && pixel_is_boundary(grid, w, h, r, c, s) {
                boundary[r * w + c] = s;
            }
        }
    }
    boundary
}
#[inline]
pub fn pixel_is_boundary(grid: &[u32], w: usize, h: usize,
                     r: usize, c: usize, s: u32) -> bool {
    for (dr, dc) in MOORE {
        let nr = r as i32 + dr;
        let nc = c as i32 + dc;
        let nb = if nr < 0 || nr >= h as i32 || nc < 0 || nc >= w as i32 {
            0u32  // grid edge counts as medium
        } else {
            grid[nr as usize * w + nc as usize]
        };
        if nb != s { return true; }
    }
    false
}
