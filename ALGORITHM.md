# Algorithm Documentation

## Table of Contents

1. [The Cellular Potts Model](#1-the-cellular-potts-model)
2. [Grid and State](#2-grid-and-state)
3. [Initialisation](#3-initialisation)
4. [The Hamiltonian](#4-the-hamiltonian)
5. [Monte Carlo Step](#5-monte-carlo-step)
6. [Boundary Detection](#6-boundary-detection)
7. [Boundary Profile Classification](#7-boundary-profile-classification)
8. [Parameters](#8-parameters)

---

## 1. The Cellular Potts Model

The Cellular Potts Model (CPM) is a lattice-based framework for simulating
multi-cell tissue dynamics.  Each lattice site (pixel) is assigned a
**spin** σ identifying which cell owns it.  Cells evolve by stochastically
copying spins from neighbouring pixels, with acceptance governed by a
Hamiltonian that encodes biological constraints (area, shape, adhesion).

The model lives on a 2-D rectangular grid of W × H pixels.

---

## 2. Grid and State

```
    col →   0   1   2   3   4   5   6
row        ┌───┬───┬───┬───┬───┬───┬───┐
 ↓  0      │ 0 │ 0 │ 1 │ 1 │ 0 │ 0 │ 0 │
           ├───┼───┼───┼───┼───┼───┼───┤
    1      │ 0 │ 1 │ 1 │ 1 │ 2 │ 0 │ 0 │
           ├───┼───┼───┼───┼───┼───┼───┤
    2      │ 0 │ 1 │ 1 │ 2 │ 2 │ 2 │ 0 │
           └───┴───┴───┴───┴───┴───┴───┘

  σ = 0  →  medium (extracellular space)
  σ = k  →  cell k  (k = 1 … n_cells)
```

Each cell k tracks:
- **area**      — number of pixels currently owned
- **perimeter** — number of 4-connected (Von Neumann) edges shared with a different σ
- **target_area**, **target_perimeter** — homeostatic targets

---

## 3. Initialisation

Two modes are available.

### Block initialisation (default)

Cells are placed as square blobs tiled in a regular grid.
Each blob has side ≈ √(target_area / 2), centred in its tile.

### Voronoi initialisation (`--voronoi-start`)

1. **Poisson-disk sampling** places n_cells seed points with a minimum
   separation distance, ensuring no two seeds are too close.
2. **Voronoi fill** assigns every non-wall pixel to its nearest seed,
   producing a partition into approximately equal, space-filling cells.

---

## 4. The Hamiltonian

The total energy is a sum of three terms:

```
H = H_contact + H_area + H_perim   [+ H_iso if λ_iso > 0]
```

### 4.1 Contact energy (adhesion)

Sums over all pairs of neighbouring (Moore-adjacent) pixels with different σ:

```
H_contact = Σ_{<i,j>} J(σ_i, σ_j)
```

| Interface         | Energy J          |
|-------------------|-------------------|
| cell – medium     | `j_cell_medium`   |
| cell – cell       | `j_cell_cell`     |
| same cell         | 0                 |

Higher J means less favourable contact → cells prefer to round up or
cluster depending on relative J values.

### 4.2 Area constraint

Quadratic penalty keeping each cell's area near its target A_t:

```
H_area = λ_A · Σ_k  (A_k − A_t)²
```

### 4.3 Perimeter constraint

Quadratic penalty keeping each cell's 4-connected perimeter near its
target P_t:

```
H_perim = λ_P · Σ_k  (P_k − P_t)²
```

The target perimeter is derived from the target area assuming a circle:
`P_t = √(4π · A_t)`.  Increasing λ_P drives cells towards rounder shapes.

### 4.4 Isoperimetric term (optional, λ_iso usually ≈ 0)

Penalises excess above the circular minimum 4π:

```
H_iso = λ_iso · Σ_k  max(P_k²/A_k − 4π, 0)
```

From the parameter sweep this term has negligible independent effect
once λ_P is non-zero; it can be set to 0.

---

## 5. Monte Carlo Step

One **Monte Carlo Step (MCS)** consists of W × H individual flip attempts.

### Single attempt

```
┌──────────────────────────────────────────────────────┐
│  1. Pick pixel (r, c) uniformly at random            │
│                                                      │
│  2. Skip if (r,c) is a wall pixel                    │
│                                                      │
│  3. s_old ← σ(r, c)                                  │
│                                                      │
│  4. Collect Moore neighbours with σ ≠ s_old          │
│     (these are the candidate new spins)              │
│     If none → abort                                  │
│                                                      │
│  5. Pick one candidate s_new at random               │
│                                                      │
│  6. Compute ΔH = ΔH_contact + ΔH_area               │
│                + ΔH_perim  + ΔH_iso                 │
│                                                      │
│  7. Metropolis acceptance:                           │
│       if ΔH ≤ 0  → accept always                    │
│       else        → accept with prob exp(−ΔH / T)   │
│                                                      │
│  8. If accepted:                                     │
│       σ(r,c) ← s_new                                │
│       update area counts for s_old, s_new            │
│       recompute perimeter for s_old, s_new           │
│       update boundary map locally                    │
└──────────────────────────────────────────────────────┘
```

### ΔH computation (O(1) per attempt)

**Contact:**  scan the 8 Moore neighbours of (r,c), sum the change in J.

**Area:**  only s_old loses one pixel, s_new gains one:
```
ΔH_area = λ_A · [(A_old − 1 − A_t)² − (A_old − A_t)²]
         + λ_A · [(A_new + 1 − A_t)² − (A_new − A_t)²]
```

**Perimeter:**  scan 4 Von Neumann neighbours of (r,c); count how many
edges are added/removed for s_old and s_new:
```
ΔH_perim = λ_P · [(P_old + ΔP_old − P_t)² − (P_old − P_t)²]
          + λ_P · [(P_new + ΔP_new − P_t)² − (P_new − P_t)²]
```

### Temperature T

T controls the noise level.  At T = 0 only energy-lowering moves are
accepted.  Typical values are T = 10–50 for adipocyte simulations.

---

## 6. Boundary Detection

A pixel (r, c) of cell σ is a **boundary pixel** if at least one of its
4 Von Neumann neighbours belongs to a different σ (including medium and wall).

The boundary map is updated locally after every accepted flip: only the
flipped pixel and its 8 Moore neighbours need to be re-evaluated.

---

## 7. Boundary Profile Classification

After simulation, the shape of each cell's boundary is quantified by
classifying every boundary pixel into one of six local pattern types.

### 7.1 Rotation to canonical orientation

For each boundary pixel ★ of cell σ:

**Step 1.** Find the first Von Neumann (N/E/S/W) neighbour that belongs
to a **different** cell.  Call its direction **north**.

```
Example: the pixel above ★ belongs to another cell → north = up
```

**Step 2.** Rotate the entire 3×3 Moore neighbourhood so that "north"
always points up.  This makes the classification rotation-invariant.

```
Before rotation          After rotation (north = up)
(north happened to       (north is always position N)
 be pointing East)

  ?  |  ★  |  ?            ?  |  N  |  ?
  N  |     |               ?  |  ★  |  ?
  ?  |     |  ?            ?  |     |  ?
```

**Step 3.** Build a **7-bit key** encoding which of the 7 non-north
positions belong to the same cell as ★ (1 = same, 0 = different).

```
Rotated 3×3 neighbourhood and bit assignments:

  ┌────────┬────────┬────────┐
  │ NW     │   N    │ NE     │
  │ bit 6  │ always │ bit 5  │
  │        │   0    │        │
  ├────────┼────────┼────────┤
  │  W     │   ★    │  E     │
  │ bit 4  │(centre)│ bit 3  │
  ├────────┼────────┼────────┤
  │ SW     │   S    │ SE     │
  │ bit 2  │ bit 1  │ bit 0  │
  └────────┴────────┴────────┘

  Key = bit6·NW + bit5·NE + bit4·W + bit3·E + bit2·SW + bit1·S + bit0·SE
  (N is always 0, not stored)
```

### 7.2 Pattern classification

The key is decoded using just the three cardinal non-north bits **W**, **E**, **S**,
plus `same_count` (total 1-bits in the key).

---

#### FlatEdge — straight boundary

Condition: `W = 1, E = 1, S = 1`

```
  ░ │ N │ ░        ░ = different cell / medium
  █ │ ★ │ █        █ = same cell
  ? │ █ │ ?
```

Both lateral neighbours and the south neighbour belong to the same cell.
The pixel sits on a smooth, straight section of the boundary.

---

#### OuterCorner — convex protrusion

Condition: `W = 0, E = 0, S = 1`

```
  ░ │ N │ ░
  ░ │ ★ │ ░
  ? │ █ │ ?
```

Only the south neighbour is same-cell.  The pixel is at the tip of a
convex corner sticking outward.

---

#### InnerCorner — concave indentation

Three sub-cases:

**Left open** (`W = 1, E = 0, S = 1`):
```
  ░ │ N │ ░
  █ │ ★ │ ░
  ? │ █ │ ?
```

**Right open** (`W = 0, E = 1, S = 1`):
```
  ░ │ N │ ░
  ░ │ ★ │ █
  ? │ █ │ ?
```

**South open** (`W = 1, E = 1, S = 0`):
```
  ░ │ N │ ░
  █ │ ★ │ █
  ? │ ░ │ ?
```

The pixel is on a concave bend; one cardinal direction has escaped the cell.

---

#### ConcaveBay — deep concavity

Condition: `same_count ≥ 5`

```
  █ │ N │ █
  █ │ ★ │ █
  ? │ █ │ ?
```

Most neighbours belong to the same cell.  The pixel lies inside a deep
bay that is almost fully enclosed.

---

#### FilamentTip — isolated spike tip

Condition: `same_count = 0`  or  `same_count = 1 and S = 1`

```
  ░ │ N │ ░        ░ │ N │ ░
  ░ │ ★ │ ░        ░ │ ★ │ ░
  ? │ ░ │ ?        ? │ █ │ ?
  (count=0)        (count=1, only S)
```

The pixel is essentially isolated — a 1-pixel-wide spike tip.
High filament_tip fraction indicates pathological cell shapes.

---

#### FilamentNeck — narrow bridge

Condition: `W = 0, E = 0, S = 0`  but  `same_count > 0`
(only diagonal neighbours are same-cell)

```
  █ │ N │ █
  ░ │ ★ │ ░
  █ │ ░ │ █
```

The cell is connected only through diagonal pixels — a narrow 1-pixel
isthmus.  Also indicates pathological morphology.

---

### 7.3 Convexity score

```
profile_convexity = (flat_edge + outer_corner) / total_boundary
```

A value near 1.0 indicates a smooth, convex cell.  Decreasing values
reflect increasing concavity and irregularity.

---

### 7.4 Summary of classification logic

```
same_count = number of bits set in the 7-bit key

W  E  S  │ same_count │  Class
─────────┼────────────┼──────────────
 1  1  1 │  any       │  FlatEdge
 0  0  1 │  any       │  OuterCorner
 *  *  * │  0         │  FilamentTip
 *  *  * │  1, S only │  FilamentTip
 *  *  * │  ≥ 5       │  ConcaveBay
 1  0  1 │  *         │  InnerCorner
 0  1  1 │  *         │  InnerCorner
 1  1  0 │  *         │  InnerCorner
 0  0  0 │  > 0       │  FilamentNeck  (diagonal-only)
 *  *  * │  *         │  Other
```

---

## 8. Parameters

| Parameter         | Symbol  | Default | Effect                                      |
|-------------------|---------|---------|---------------------------------------------|
| `lambda_area`     | λ_A     | 1.0     | Area constraint strength; binary on/off     |
| `lambda_perim`    | λ_P     | 0.2     | **Dominant shape parameter**                |
| `lambda_iso`      | λ_iso   | 0.1     | Isoperimetric; negligible — set to 0        |
| `target_area`     | A_t     | 200     | Homeostatic area per cell (pixels)          |
| `temperature`     | T       | 10.0    | Metropolis noise level                      |
| `j_cell_cell`     | J_cc    | 10.0    | Cell–cell adhesion cost                     |
| `j_cell_medium`   | J_cm    | 15.0    | Cell–medium adhesion cost                   |
| `wall_inset`      | —       | 0       | Frozen border thickness                     |
| `n_cells`         | n       | 16      | Number of cells                             |
| `grid_w / grid_h` | W, H   | 60      | Grid dimensions                             |

### Parameter sweep findings

From systematic sweeps of (λ_A, λ_P, λ_iso):

- **λ_P is the dominant parameter.**  Increasing λ_P progressively
  rounds cells and reduces filament/concave fractions.
- **λ_A acts as a binary switch.**  Any non-zero value prevents cell
  collapse; the specific magnitude has minimal effect.
- **λ_iso is negligible** once λ_P is non-zero.  It can safely be
  fixed at 0.
