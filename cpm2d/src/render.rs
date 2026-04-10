// ─── Helper: collect candidate copy-source sigmas ────────────────────────────

use image::{ImageBuffer, Rgb};
use crate::grid::Cpm2d;

pub fn arrayvec_neighbours(
    r: usize, c: usize,
    W: usize, H: usize,
    grid: &[u32],
    s_old: u32,
) -> Vec<u32> {
    let mut out = Vec::with_capacity(8);
    for (dr, dc) in crate::grid::MOORE {
        let nr = r as i32 + dr;
        let nc = c as i32 + dc;
        if nr >= 0 && nr < H as i32 && nc >= 0 && nc < W as i32 {
            let nb = grid[nr as usize * W + nc as usize];
            if nb != s_old {
                out.push(nb);
            }
        }
    }
    out
}
// ─── Procedural colormap ─────────────────────────────────────────────────────

/// Convert HSV (h in [0,360), s and v in [0,1]) to RGB u8 triple.
pub fn hsv_to_rgb(h: f32, s: f32, v: f32) -> [u8; 3] {
    let c = v * s;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let m = v - c;
    let (r1, g1, b1) = match h as u32 {
        0..=59   => (c, x, 0.0),
        60..=119 => (x, c, 0.0),
        120..=179 => (0.0, c, x),
        180..=239 => (0.0, x, c),
        240..=299 => (x, 0.0, c),
        _        => (c, 0.0, x),
    };
    [
        ((r1 + m) * 255.0).round() as u8,
        ((g1 + m) * 255.0).round() as u8,
        ((b1 + m) * 255.0).round() as u8,
    ]
}

/// Return the RGB colour for sigma `s` out of `n_cells` total cells.
/// sigma 0 (medium) is a neutral dark grey.
/// Cells are spread evenly around the hue wheel, with alternating
/// saturation (1.0 / 0.65) and value (0.9 / 1.0) so adjacent indices
/// are visually distinct even when n is large.
pub fn cell_colour(s: usize, n_cells: usize) -> [u8; 3] {
    if s == 0 {
        return [45, 45, 45]; // medium: dark grey
    }
    let k = s - 1; // 0-based cell index
    let n = n_cells.max(1);
    let hue = (k as f32 / n as f32) * 360.0;
    // Alternate saturation and brightness so neighbouring hues stay distinct
    let sat = if k % 2 == 0 { 0.85 } else { 0.55 };
    let val = if k % 2 == 0 { 1.00 } else { 0.90 };
    hsv_to_rgb(hue, sat, val)
}
pub fn console_print(grid : &Cpm2d) {
    // ANSI background colours: 0=reset, 1=red, 2=green, 3=yellow, 4=blue, …
    const BG: &[&str] = &[
        "\x1b[0m  ",   // 0 medium
        "\x1b[41m  ",  // 1 red
        "\x1b[42m  ",  // 2 green
        "\x1b[43m  ",  // 3 yellow
        "\x1b[44m  ",  // 4 blue
        "\x1b[45m  ",  // 5 magenta
        "\x1b[46m  ",  // 6 cyan
        "\x1b[47m  ",  // 7 white
        "\x1b[101m  ", // 8 bright red
    ];
    let reset = "\x1b[0m";
    let W = grid.p.grid_w;
    let H = grid.p.grid_h;
    let n = grid.p.n_cells;

    println!("\n─── MCS {} ───", grid.mcs);
    for r in 0..H {
        for c in 0..W {
            let s = grid.grid[r * W + c] as usize;
            let [red, grn, blu] = cell_colour(s, n);
            print!("\x1b[48;2;{red};{grn};{blu}m  ");
            // let idx = s.min(BG.len() - 1);
            //print!("{}", BG[idx]);
        }
        println!("{}", reset);
    }
    for k in 1..=grid.p.n_cells {
        println!(
            "  cell {}: area={:4} (tgt {:3})  perim={:4} (tgt {:3})",
            k, grid.cells[k].area, grid.p.target_area,
            grid.cells[k].perimeter, grid.p.target_perim
        );
    }
}

// ── PNG output ────────────────────────────────────────────────────────────
pub fn save_png(grid : &Cpm2d, path: Option<&str>) {
    let default = format!("{}/frame_{:06}.png", grid.p.png_dir, grid.mcs);
    let path = path.unwrap_or(&default);

    let scale = 8u32;
    let W = grid.p.grid_w as u32;
    let H = grid.p.grid_h as u32;
    let img_w = W * scale;
    let img_h = H * scale + 20;

    let mut img = ImageBuffer::<Rgb<u8>, _>::new(img_w, img_h);

    // Fill background (label bar)
    for px in img.pixels_mut() {
        *px = Rgb([30u8, 30, 30]);
    }

    let n = grid.p.n_cells;

    // ── Pass 1: cell interiors ────────────────────────────────────────────
    for r in 0..H {
        for c in 0..W {
            let s = grid.grid[(r * W + c) as usize] as usize;
            let rgb = Rgb(cell_colour(s, n));
            for dy in 0..scale {
                for dx in 0..scale {
                    img.put_pixel(c * scale + dx, r * scale + dy, rgb);
                }
            }
        }
    }

    // ── Pass 2: boundary overlay (darker shade) ───────────────────────────
    const DARKEN: f32 = 0.55; // 1.0 = no change, 0.0 = black
    for r in 0..H {
        for c in 0..W {
            let s = grid.boundary[(r * W + c) as usize] as usize;
            if s == 0 { continue; }
            let base = cell_colour(s, n);
            let dark = Rgb(base.map(|ch| (ch as f32 * DARKEN) as u8));
            for dy in 0..scale {
                for dx in 0..scale {
                    img.put_pixel(c * scale + dx, r * scale + dy, dark);
                }
            }
        }
    }

    img.save(path).unwrap_or_else(|e| eprintln!("PNG save error: {e}"));
    println!("  [png] saved → {path}");
}