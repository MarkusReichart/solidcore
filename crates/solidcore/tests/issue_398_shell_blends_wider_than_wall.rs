//! tangency-labs/solidcore#398 (probe): shelling 2 mm a box with r10 rounds, a 5 mm
//! chamfer and a 4 mm fillet (every blend larger than the wall) fails with 1001 in the
//! SMLib fork (rig case shell/1525). Not fixed in the fork: this probe records what the
//! Rust port does at each stage.
//!
//! The rig's part is not public, so this is a reconstruction with the three blend
//! kinds: an L block (80 x 60 x 20 with an 80 x 30 x 20 step on its back half), its
//! four vertical outer edges rounded r10, the concave step edge filleted r4 and the
//! front top edge chamfered 5. There is no closed form for the shell, so the probe
//! checks for a manifold solid lighter than its input.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use offset_fixture::{report, shell, volume};
use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::boolean::api_brep::boolean_union;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::{Brep, EdgeId, Kernel};

/// The edge whose curve's middle is within 1e-6 of `at`.
fn edge_at(k: &Kernel, b: &Brep, at: Point3) -> EdgeId {
    b.edges()
        .iter()
        .find(|(_, e)| {
            let Some(c) = e.curve.as_ref() else {
                return false;
            };
            let mid = c.evaluate_point(k, e.interval.evaluate(0.5)).unwrap();
            mid.distance(at) < 1e-6
        })
        .map(|(id, _)| id)
        .unwrap_or_else(|| panic!("no edge at {at:?}"))
}

/// The L block, rounded, filleted and chamfered; stops at `stage` (0 block, 1 rounds,
/// 2 fillet, 3 chamfer).
fn part(k: &mut Kernel, stage: u32) -> solidcore::Result<Brep> {
    let base = create_box(k, Point3::new(0.0, 0.0, 0.0), 80.0, 60.0, 20.0).unwrap();
    let step = create_box(k, Point3::new(0.0, 30.0, 20.0), 80.0, 30.0, 20.0).unwrap();
    let mut b = boolean_union(k, base, step).unwrap();
    if stage >= 1 {
        let rounds = [
            edge_at(k, &b, Point3::new(0.0, 0.0, 10.0)),
            edge_at(k, &b, Point3::new(80.0, 0.0, 10.0)),
            edge_at(k, &b, Point3::new(0.0, 60.0, 20.0)),
            edge_at(k, &b, Point3::new(80.0, 60.0, 20.0)),
        ];
        fillet_edges(k, &mut b, &rounds, 10.0, 1, 1, 1.0)?;
    }
    if stage >= 2 {
        let step_edge = edge_at(k, &b, Point3::new(40.0, 30.0, 20.0));
        fillet_edges(k, &mut b, &[step_edge], 4.0, 1, 1, 1.0)?;
    }
    if stage >= 3 {
        // the front top edge runs between the front rounds' rails, x 10 to 70
        let front = edge_at(k, &b, Point3::new(40.0, 0.0, 20.0));
        fillet_edges(k, &mut b, &[front], 5.0, 0, 1, 1.0)?;
    }
    Ok(b)
}

#[test]
fn issue_398_blended_input() {
    let mut k = Kernel::default();
    for stage in 0..=3 {
        match part(&mut k, stage) {
            Ok(mut b) => eprintln!(
                "PROBE issue_398 stage {stage}: faces={} volume={:?}",
                b.face_list().len(),
                volume(&k, &mut b)
            ),
            // a blend the fillet cannot build yet: the later stages need it
            Err(e) => {
                eprintln!("PROBE issue_398 stage {stage}: blend Err({e:?})");
                break;
            }
        }
    }
}

#[test]
fn issue_398_shell_each_stage() {
    let mut all = true;
    for stage in 0..=3 {
        let mut k = Kernel::default();
        let mut b = match part(&mut k, stage) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("PROBE issue_398 stage {stage}: blend Err({e:?}), later stages not shelled");
                break;
            }
        };
        let input = volume(&k, &mut b).unwrap();
        let result = shell(&mut k, b, -2.0, &[]);
        let manifold = report(&k, &format!("issue_398 stage {stage} closed"), result, None, 0.0);
        all &= manifold;
        eprintln!("PROBE issue_398 stage {stage}: input volume {input}");
    }
    assert!(all, "a stage's shell is not a manifold solid (see PROBE lines)");
}
