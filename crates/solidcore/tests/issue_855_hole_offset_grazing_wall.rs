//! tangency-labs/solidcore#855 (probe): with #854's crossing loops resolved, a plate
//! whose hole's offset crosses the wall's offset by only 0.01 still fails in the SMLib
//! fork with 1001 (a grazing SSI piece at a triple point becomes a sliver ProtoEdge).
//! Not fixed in the fork: this probe records what the Rust port does, with and without
//! the #854 opt-in (`Kernel::fixes.offset_resolve_crossing_loops`).
//!
//! A 40 x 40 x 10 plate with an r2 through hole at (4.99, 20), shelled 1.5 inward with
//! the top open: the hole's offset (r 3.5) reaches x = 1.49, 0.01 past the side's offset
//! at x = 1.5. The cavity is the inner box 37 x 37 x 8.5 less the part of the hole
//! offset's disc inside it.
//!
//! In the port, with the opt-in, a crossing shallower than the protective stitch's
//! largest tolerance had the dropped lobe stitched shut and kept; the stitch now leaves
//! the resolution's open edges open
//! (`plan/findings/offset.offset_executive_p2--issue-855-stitch-closes-dropped-lobe.md`).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use core::f64::consts::PI;

use offset_fixture::{planar_faces_at, report, shell};
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::{Brep, Kernel};

const HOLE_X: f64 = 4.99;

fn plate(k: &mut Kernel, hole_x: f64) -> Brep {
    let plate = create_box(k, Point3::new(0.0, 0.0, 0.0), 40.0, 40.0, 10.0).unwrap();
    let hole = create_cylinder(k, Point3::new(hole_x, 20.0, -1.0), 2.0, 12.0).unwrap();
    boolean_difference(k, plate, hole).unwrap()
}

fn want(hole_x: f64) -> f64 {
    let t: f64 = 1.5;
    let r = 2.0 + t;
    let d = hole_x - t; // the hole offset's centre from the side's offset
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    let disc_inside = PI * r * r - segment;
    let cavity = ((40.0 - 2.0 * t) * (40.0 - 2.0 * t) - disc_inside) * (10.0 - t);
    40.0 * 40.0 * 10.0 - PI * 4.0 * 10.0 - cavity
}

fn run(resolve_crossing_loops: bool) -> bool {
    run_at(resolve_crossing_loops, HOLE_X).0
}

/// Shells the plate with its hole at `hole_x`; whether the result matches the oracle,
/// and its face count.
fn run_at(resolve_crossing_loops: bool, hole_x: f64) -> (bool, Option<usize>) {
    let mut k = Kernel::default();
    k.fixes.offset_resolve_crossing_loops = resolve_crossing_loops;
    let b = plate(&mut k, hole_x);
    let top = planar_faces_at(&k, &b, 2, 10.0);
    assert_eq!(top.len(), 1, "the top face");
    let result = shell(&mut k, b, -1.5, &top);
    let faces = result.as_ref().ok().map(|b| b.face_list().len());
    let ok = report(
        &k,
        &format!(
            "issue_855 open top, hole at x = {hole_x}, resolve_crossing_loops={resolve_crossing_loops}"
        ),
        result,
        Some(want(hole_x)),
        1e-6,
    );
    (ok, faces)
}

#[test]
fn issue_855_default() {
    // upstream behaviour: reported only (the loops cross, as in #854)
    let _ = run(false);
}

#[test]
fn issue_855_with_crossing_loops_resolved() {
    assert!(run(true));
}

/// The stitch threshold (#855): with the #854 opt-in, a hole offset crossing the side's
/// offset by less than the protective stitch's largest tolerance (about 0.26 here: 1/200
/// of the merged offset's box diagonal) used to have its dropped lobe stitched shut and
/// kept (11 faces, not a solid); the crossing-loop resolution's open edges now stay
/// open, so every depth gives the 15-face shell.
#[test]
fn issue_855_with_crossing_loops_resolved_by_depth() {
    for crossing in [0.5, 0.28, 0.25, 0.1, 0.001] {
        let hole_x = 1.5 + 3.5 - crossing;
        let (ok, faces) = run_at(true, hole_x);
        assert!(ok, "crossing {crossing}: not the expected shell");
        assert_eq!(faces, Some(15), "crossing {crossing}: face count");
    }
}

/// Narrows #855 (reported only): the same plate with the hole moved so its offset
/// crosses the side's offset by 0.5, 0.1, 0.03, 0.01 and 0.001, shelled open at the top
/// with the #854 opt-in; per case the result's faces (surface type, area, middle point)
/// and its non-manifold edges.
#[test]
#[ignore = "diagnostic probe for tangency-labs/solidcore#855, run on demand"]
fn issue_855_by_crossing_depth() {
    use solidcore::mass::integrator::properties::face_area;
    use solidcore::topology::edge::is_manifold;
    for crossing in [0.5, 0.28, 0.25, 0.1, 0.03, 0.01, 0.001] {
        let mut k = Kernel::default();
        k.fixes.offset_resolve_crossing_loops = true;
        let x = 1.5 + 3.5 - crossing;
        let plate_body = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 40.0, 40.0, 10.0).unwrap();
        let hole = create_cylinder(&mut k, Point3::new(x, 20.0, -1.0), 2.0, 12.0).unwrap();
        let b = boolean_difference(&mut k, plate_body, hole).unwrap();
        let top = planar_faces_at(&k, &b, 2, 10.0);
        let result = shell(&mut k, b, -1.5, &top);
        let mut shelled = match result {
            Ok(s) => s,
            Err(e) => {
                eprintln!("PROBE issue_855 crossing {crossing}: Err({e:?})");
                continue;
            }
        };
        let faces: Vec<_> = shelled.face_list().to_vec();
        eprintln!(
            "PROBE issue_855 crossing {crossing}: faces={} edges={} vertices={}",
            faces.len(),
            shelled.edge_list().len(),
            shelled.vertex_list().len()
        );
        for f in faces {
            let (code, mid) = {
                let face = shelled.face(f).unwrap();
                let s = face.surface.as_ref().unwrap();
                (s.type_code(), s.evaluate_point(&k, face.uv_domain.mid()).ok())
            };
            let area = face_area(&k, &mut shelled, f, 1e-6);
            eprintln!("  face {f:?}: type {code} area {area:?} mid {mid:?}");
        }
        let edges: Vec<_> = shelled.edge_list().to_vec();
        for e in edges {
            if !is_manifold(&shelled, e).unwrap_or(false) {
                let edge = shelled.edge(e).unwrap();
                let c = edge.curve.as_ref();
                let ends = c.map(|c| {
                    (
                        c.evaluate_point(&k, edge.interval.min()).ok(),
                        c.evaluate_point(&k, edge.interval.max()).ok(),
                    )
                });
                eprintln!("  non-manifold edge {e:?}: uses {} ends {ends:?}", edge.edgeuses.len());
            }
        }
    }
}
