//! tangency-labs/solidcore#879: a blend as wide as the face it rolls on (a radius-10
//! cylinder's top rim blended 10 across the top). Upstream builds it; the chamfer's cone
//! then ends in a point at the top's centre. With
//! `Kernel::fixes.fillet_refuse_face_consuming_blends` on, a blend whose cross section meets
//! the face at an angle (the linear section, a chamfer) is refused with
//! `FilletErrorType::NoRoom`; a round's arc is tangent to the face and passes the rule
//! (the SolidCraft fork's `SmFilletExecutive::FacesOutlastTheirBlends` and its
//! amendment).

use solidcore::blend::fillet_executive::session::FilletSession;
use solidcore::blend::fillet_geom::p1::FilletErrorType;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_cylinder;
use solidcore::topology::{Brep, EdgeId};
use solidcore::Kernel;

/// The closed rim edge of a cylinder about Z at `height`.
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

/// Blends the top rim of a radius-10, height-20 cylinder at `width` with the constant
/// radius solver (`xsect`: 0 linear, 1 circular) up to the intermediate fillet body;
/// whether that was built, and the executive's error code.
fn blend_top_rim(width: f64, xsect: u32, refuse: bool) -> (bool, i32, String) {
    let mut k = Kernel::default();
    k.fixes.fillet_refuse_face_consuming_blends = refuse;
    let mut target = create_cylinder(&mut k, Point3::ZERO, 10.0, 20.0).unwrap();
    let top = rim(&target, 20.0);
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[top], width, xsect, 1, 1.0)
        .unwrap();
    session.create_corners(&k, None).unwrap();
    let result = session.create_fillet_brep(&mut k, None);
    (result.is_ok(), session.executive().error.code, format!("{result:?}"))
}

#[test]
fn issue_879_chamfer_taking_the_whole_top_is_refused() {
    let (built, code, result) = blend_top_rim(10.0, 0, true);
    assert!(!built, "the chamfer as wide as the top is refused");
    assert_eq!(code, FilletErrorType::NoRoom as i32, "{result}");
}

#[test]
fn issue_879_round_taking_the_whole_top_passes_the_rule() {
    // the round's arc is tangent to the top at its rail: not refused by this rule
    let (_, code, _) = blend_top_rim(10.0, 1, true);
    assert_ne!(code, FilletErrorType::NoRoom as i32);
}

#[test]
fn issue_879_a_chamfer_short_of_the_whole_top_passes_the_rule() {
    let (_, code, _) = blend_top_rim(9.99, 0, true);
    assert_ne!(code, FilletErrorType::NoRoom as i32);
}

#[test]
fn issue_879_default_keeps_upstream() {
    // without the opt-in nothing is refused for taking the whole of a face
    let (_, code, _) = blend_top_rim(10.0, 0, false);
    assert_ne!(code, FilletErrorType::NoRoom as i32);
}
