//! Chamfers (constant-distance fillets) on a [`FilletSession`]: `SmApiChamferFillet`.
//!
//! Upstream: `source/SM_API/src/SmApiFillets.cpp` (`SmApiChamferFillet`).
//!
//! | Rust               | Upstream             |
//! |--------------------|----------------------|
//! | [`chamfer_fillet`] | `SmApiChamferFillet` |
//!
//! Deviations:
//! - `SmApiChamferFillet` sets an `SmMakeSurfaceBlendSIH(1.0)` self-intersection
//!   handler; the port's executive has no handler objects yet, so the solvers' handler
//!   read stays [`SolidError::Unported`] (`SmFilletExecutive::m_pSelfIntersectionHandler`,
//!   unit `blend.fillet_executive_p1`) and a chamfer whose rail traces self-intersect
//!   defers there instead of being blended.
//! - No recursion.

use super::FilletSession;
use crate::blend::fillet_solver::p1::FilletSurfaceGenerator;
use crate::error::Result;
use crate::kernel::Kernel;
use crate::math::PI;
use crate::topology::Brep;

/// `SmApiChamferFillet(pBrep, dDistance)`: every edge of `brep` that has a blend edge
/// use gets an `SmConstantDistanceFS` (the B-rep's tolerance, 30 and 2 degrees,
/// `distance` the chamfer width across the chamfer face) with one linear
/// cross-section generator, loaded in edge order; then `CreateFilletCorners`, every
/// corner a bevel, and `DoFilleting` with the global merge.
///
/// # Errors
/// Those of the solver construction and loading, `CreateFilletCorners` and
/// `DoFilleting`; the mutations of the stages before a failing one are kept.
pub fn chamfer_fillet(k: &mut Kernel, brep: &mut Brep, distance: f64) -> Result<()> {
    let tolerance = brep.tolerance();
    let edges = crate::topology::brep::p1::edges(brep, None)?;
    let mut session = FilletSession::new(k, brep)?;

    // SmLinearCrossSectionFSG sFSGLinear (one generator shared by every solver)
    let generator = session.add_surface_generator(FilletSurfaceGenerator::Linear)?;
    session.executive.set_do_global_merge(true);

    for e in edges {
        let Some(eu) = crate::topology::edge::blend_edgeuse(k, &mut *session.objects.target, e)?
        else {
            continue;
        };
        let solver = session.new_constant_distance_solver(
            k,
            eu,
            distance,
            tolerance,
            30.0 * PI / 180.0,
            2.0 * PI / 180.0,
        )?;
        session.set_solver_generator(solver, generator)?;
        session.load_solver(k, solver)?;
    }

    session.create_corners(k, None)?;
    session.set_corners_bevel();
    session.do_filleting(k, None)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::attribute::{AttributeStore, ObjectKey};
    use crate::mass::integrator::properties::{brep_volume, PropertyAttributes};
    use crate::math::Point3;

    fn volume(k: &Kernel, b: &mut Brep) -> f64 {
        let store = AttributeStore::default();
        let key_of = |_| ObjectKey(0);
        let attributes = PropertyAttributes {
            store: &store,
            key_of: &key_of,
        };
        brep_volume(k, b, 1e-8, &attributes).unwrap()
    }

    /// One box edge chamfered through the session's constant-distance entry (linear
    /// cross section): seven faces, and the box less a right prism of legs
    /// `L / sqrt(2)` along the edge.
    #[test]
    fn chamfer_one_box_edge_end_to_end() {
        let mut k = Kernel::default();
        let mut b = crate::primitive::api_primitives::create_box(
            &mut k,
            Point3::new(0.0, 0.0, 0.0),
            10.0,
            20.0,
            30.0,
        )
        .unwrap();
        let e = b.edges().ids().next().unwrap();
        let length = {
            let edge = b.edge(e).unwrap();
            let curve = edge.curve.as_ref().unwrap();
            let ivl = edge.interval;
            let p0 = curve.evaluate_point(&k, ivl.min()).unwrap();
            let p1 = curve.evaluate_point(&k, ivl.max()).unwrap();
            p0.distance(p1)
        };
        let distance = 1.0;
        {
            let mut session = FilletSession::new(&mut k, &mut b).unwrap();
            session
                .prepare_constant_distance_edges(&mut k, &[e], distance, 0, 1, 1.0)
                .unwrap();
            session.do_fillet(&mut k, None).unwrap();
        }
        crate::topology::brep::p5::stitch_and_orient(
            &mut k,
            &mut b,
            &crate::topology::brep::p5::StitchAndOrientOptions::default(),
        )
        .unwrap();
        let faces = crate::topology::brep::p1::faces(&b, None).unwrap().len();
        assert_eq!(faces, 7);
        assert!(crate::topology::brep::p2::is_manifold_solid(&b).unwrap());
        let want = 6000.0 - distance * distance / 4.0 * length;
        let got = volume(&k, &mut b);
        assert!((got - want).abs() < 1e-6 * want, "volume {got}, want {want}");
    }
}
