//! Trace (tangency-labs/solidcore#919, #920, #879): the stages of a constant-radius blend
//! of a radius-10 cylinder's top rim through the fillet session, for each cross section
//! and width. Run with `--ignored --nocapture`; it asserts nothing.

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_geom::p1::FilletSolverId;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_cylinder;
use solidcore::topology::{Brep, EdgeId};
use solidcore::Kernel;

fn rim(b: &Brep, height: f64) -> EdgeId {
    b.edges()
        .iter()
        .find(|(_, e)| {
            e.vertices[0] == e.vertices[1]
                && (b.vertex(e.vertices[0]).unwrap().point.z - height).abs() < 1e-9
        })
        .map(|(id, _)| id)
        .unwrap()
}

fn stages(width: f64, xsect: u32, continuity: u32) {
    let mut k = Kernel::default();
    let mut target = create_cylinder(&mut k, Point3::ZERO, 10.0, 20.0).unwrap();
    let top = rim(&target, 20.0);
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    let tag = format!("width {width} xsect {xsect} continuity {continuity}");
    let r = session.prepare_constant_radius_edges(&mut k, &[top], width, xsect, continuity, 1.0);
    println!("{tag}: prepare {r:?}");
    if r.is_err() {
        return;
    }
    let r = session.create_corners(&k, None);
    println!("{tag}: create_corners {r:?}");
    if r.is_err() {
        return;
    }
    macro_rules! stage {
        ($name:literal, $call:expr) => {{
            let r = $call;
            println!("{tag}: {} {:?}", $name, r.as_ref().map(|_| ()));
            if r.is_err() {
                println!("{tag}: error info {:?}", session.executive().error);
                return;
            }
        }};
    }
    stage!("create_corner_vertices", session.create_corner_vertices(&mut k));
    stage!("adjust_corner_topology", session.adjust_corner_topology(&mut k));
    // The inner error the executive stage reduces to SM_FS_SURF_TRACING_FAILURE.
    let r = session.calc_solver_geometry(&k, FilletSolverId(0));
    println!("{tag}: calc_solver_geometry(0) {r:?}");
    if let Ok(solver) = session.solver(FilletSolverId(0)) {
        println!("{tag}: solver 0 status {:?}", solver.base().status);
    }
    if r.is_err() {
        return;
    }
    stage!("calculate_solver_geometry", session.calculate_solver_geometry(&mut k));
    stage!("calculate_corner_edges", session.calculate_corner_edges(&mut k));
    stage!("trim_solver_rails", session.trim_solver_rails(&mut k));
    stage!("classify_solver_rails", session.classify_solver_rails(&mut k));
    stage!("build_fillet_faces", session.build_fillet_faces(&mut k));
    stage!("stitch_fillet_faces", session.stitch_fillet_faces(&mut k));
    if let Some(body) = session.intermediate_brep() {
        println!(
            "{tag}: intermediate body v/e/f {} {} {}",
            body.vertices().len(),
            body.edges().len(),
            body.faces().len()
        );
    }
}

#[test]
#[ignore = "trace for #919/#920; run on request"]
fn trace_cylinder_rim_blends() {
    for width in [9.99, 10.0] {
        stages(width, 0, 1); // linear (chamfer-like)
        stages(width, 1, 1); // circular (round)
        stages(width, 2, 2); // G2 blend curve
    }
}
