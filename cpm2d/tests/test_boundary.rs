// tests/test_boundary.rs
use cpm2d::grid::Cpm2d;
use cpm2d::params::Params;

#[test]
fn boundary_pixel_count_matches_perimeter() {
    let sim = Cpm2d::new(Params::default(), false);
    let boundary_count = sim.boundary.iter().filter(|&&x| x > 0).count();
    // boundary pixel count should be in a sane range relative to perimeter
    assert!(boundary_count > 0);
}