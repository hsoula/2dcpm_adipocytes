/// Place `n_cells` cubic blobs in a regular 3-D grid layout.
/// Each blob is a cube of side `blob_side` voxels centred in its tile.
pub fn place_cells_grid(
    grid: &mut Vec<u32>,
    w: usize, h: usize, d: usize,
    n_cells: usize,
    blob_side: usize,
) {
    // Number of tiles per axis (ceil(n^{1/3}))
    let n_side = (n_cells as f64).cbrt().ceil() as usize;
    let tile_w = w / n_side;
    let tile_h = h / n_side;
    let tile_d = d / n_side;
    let bs = blob_side
        .max(2)
        .min(tile_w.saturating_sub(2))
        .min(tile_h.saturating_sub(2))
        .min(tile_d.saturating_sub(2));

    let mut sigma = 1u32;
    'outer: for iz in 0..n_side {
        for iy in 0..n_side {
            for ix in 0..n_side {
                if sigma as usize > n_cells { break 'outer; }

                let ox = ix * tile_w + (tile_w - bs) / 2;
                let oy = iy * tile_h + (tile_h - bs) / 2;
                let oz = iz * tile_d + (tile_d - bs) / 2;

                for dz in 0..bs {
                    for dy in 0..bs {
                        for dx in 0..bs {
                            let x = ox + dx;
                            let y = oy + dy;
                            let z = oz + dz;
                            if x < w && y < h && z < d {
                                grid[z * w * h + y * w + x] = sigma;
                            }
                        }
                    }
                }
                sigma += 1;
            }
        }
    }
}

/// Place spherical blobs with per-cell radii derived from each cell's target_volume.
/// `cells[0]` is medium (skipped); cells 1..n get their own sphere.
pub fn place_cells_spheres_individual(
    grid: &mut Vec<u32>,
    w: usize, h: usize, d: usize,
    cells: &[crate::cellstate::CellState],
) {
    let n_cells = cells.len().saturating_sub(1); // exclude medium at index 0
    if n_cells == 0 { return; }

    let n_side = (n_cells as f64).cbrt().ceil() as usize;
    let tile_w = w / n_side;
    let tile_h = h / n_side;
    let tile_d = d / n_side;

    let mut sigma = 1usize;
    'outer: for iz in 0..n_side {
        for iy in 0..n_side {
            for ix in 0..n_side {
                if sigma > n_cells { break 'outer; }

                let tv = cells[sigma].target_volume.max(1) as f64;
                let radius = ((3.0 / (4.0 * std::f64::consts::PI)) * tv).powf(1.0 / 3.0);

                let cx = (ix * tile_w + tile_w / 2) as f64;
                let cy = (iy * tile_h + tile_h / 2) as f64;
                let cz = (iz * tile_d + tile_d / 2) as f64;
                let r2 = radius * radius;

                let lo_x = ((cx - radius).floor() as usize).max(0);
                let hi_x = ((cx + radius).ceil() as usize).min(w - 1);
                let lo_y = ((cy - radius).floor() as usize).max(0);
                let hi_y = ((cy + radius).ceil() as usize).min(h - 1);
                let lo_z = ((cz - radius).floor() as usize).max(0);
                let hi_z = ((cz + radius).ceil() as usize).min(d - 1);

                for z in lo_z..=hi_z {
                    for y in lo_y..=hi_y {
                        for x in lo_x..=hi_x {
                            let dx = x as f64 - cx;
                            let dy = y as f64 - cy;
                            let dz = z as f64 - cz;
                            if dx*dx + dy*dy + dz*dz <= r2 {
                                grid[z * w * h + y * w + x] = sigma as u32;
                            }
                        }
                    }
                }
                sigma += 1;
            }
        }
    }
}

/// Place `n_cells` spherical blobs of given radius.
pub fn place_cells_spheres(
    grid: &mut Vec<u32>,
    w: usize, h: usize, d: usize,
    n_cells: usize,
    radius: f64,
) {
    let n_side = (n_cells as f64).cbrt().ceil() as usize;
    let tile_w = w / n_side;
    let tile_h = h / n_side;
    let tile_d = d / n_side;

    let mut sigma = 1u32;
    'outer: for iz in 0..n_side {
        for iy in 0..n_side {
            for ix in 0..n_side {
                if sigma as usize > n_cells { break 'outer; }

                let cx = (ix * tile_w + tile_w / 2) as f64;
                let cy = (iy * tile_h + tile_h / 2) as f64;
                let cz = (iz * tile_d + tile_d / 2) as f64;
                let r2 = radius * radius;

                let lo_x = ((cx - radius).floor() as usize).max(0);
                let hi_x = ((cx + radius).ceil() as usize).min(w - 1);
                let lo_y = ((cy - radius).floor() as usize).max(0);
                let hi_y = ((cy + radius).ceil() as usize).min(h - 1);
                let lo_z = ((cz - radius).floor() as usize).max(0);
                let hi_z = ((cz + radius).ceil() as usize).min(d - 1);

                for z in lo_z..=hi_z {
                    for y in lo_y..=hi_y {
                        for x in lo_x..=hi_x {
                            let dx = x as f64 - cx;
                            let dy = y as f64 - cy;
                            let dz = z as f64 - cz;
                            if dx*dx + dy*dy + dz*dz <= r2 {
                                grid[z * w * h + y * w + x] = sigma;
                            }
                        }
                    }
                }
                sigma += 1;
            }
        }
    }
}
