//! The `prog_test` `smTestBrep.cpp` fixtures (v1.0.1) as B-rep blocks of the shared topology
//! fixture (`harness/cpp/shim/sc_topology_fixture.h`, built by
//! `harness/cpp/shim/support/test_brep_fixture.h`), and the `-4` `SmCrvOnSurf` edge curve
//! record their exports carry.
//!
//! | Block | Upstream builder |
//! |-------|------------------|
//! | `15 case` | `sm_MakeLoopOrientationCase(lCase)` |
//! | `16 case` | `sm_MakeTorusVertexLoopCase(lCase)` (fixture `lCase % 10`) |
//! | `17 size` | `sm_make_xy_sheet(dSize)` |
//!
//! The C++ side builds the B-rep with the upstream constructors; the Rust side imports the
//! export like every other block. A runner takes a [`TestBrepIn`] beside its usual B-rep
//! input and encodes it with [`block`].

use serde::{Deserialize, Serialize};
use solidcore::curve::curve::p2::CRV_ON_SURF_TYPE;
use solidcore::curve::curve::{Curve as SmCurve, CurveKind};
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::topology::database_io::{DatabaseReader, FileType, InputStream};

use crate::HarnessError;

/// A `smTestBrep.cpp` fixture.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub(crate) enum TestBrepIn {
    /// `sm_MakeLoopOrientationCase(case)`: a torus face bounded by one closed `SmCrvOnSurf`
    /// edge (cases 8 and 9 add a hole; 3 and 4 a stale reversed UV trim).
    LoopOrientation { case: u32 },
    /// `sm_MakeTorusVertexLoopCase(case)`: a torus face with a UV-square loop and 1..5
    /// vertex loops, after the builder's `HasSeamProblem` call.
    TorusVertexLoop { case: u32 },
    /// `sm_make_xy_sheet(size)`: `SmPrimitiveCreation::CreateRectangle(1e-5, size, size)`.
    XySheet { size: f64 },
}

/// The B-rep block of `t`, with the trailing `validate_only` and `init_gaps` flags of every
/// block.
pub(crate) fn block(t: &TestBrepIn, validate_only: bool, init_gaps: bool) -> Vec<f64> {
    let mut v = match t {
        TestBrepIn::LoopOrientation { case } => vec![15.0, f64::from(*case)],
        TestBrepIn::TorusVertexLoop { case } => vec![16.0, f64::from(*case)],
        TestBrepIn::XySheet { size } => vec![17.0, *size],
    };
    v.push(f64::from(u8::from(validate_only)));
    v.push(f64::from(u8::from(init_gaps)));
    v
}

fn harness(message: impl Into<String>) -> HarnessError {
    HarnessError::new(message)
}

/// Reads an `SmCrvOnSurf` record after its `-4` marker: `type version dim nbytes bytes...`,
/// the curve's binary database record, read back with `SmCrvOnSurf::ReadFromDB` (owner flag
/// 2), then `has_approx [matches max_gap dim n p m (x y z w)*(n+1) knot*(m+1)]`: the cached
/// B-spline approximation the record leaves out. The import rebuilds that cache as upstream
/// made it (`GetOrCreateBSApproxPointer`'s defaults: `SM_APPROX_TOL_3D`, the breaks of the
/// record, the exported parameterisation flag) and checks the rebuilt curve and gap bit for
/// bit against the exported ones. Counts, bytes, kind and complete consumption are checked
/// before success.
pub(crate) fn read_crv_on_surf(
    num: &mut dyn FnMut() -> Result<f64, HarnessError>,
    remaining: usize,
) -> Result<SmCurve, HarnessError> {
    let integer = |num: &mut dyn FnMut() -> Result<f64, HarnessError>,
                   max: f64|
     -> Result<u64, HarnessError> {
        let n = num()?;
        if !n.is_finite() || !(0.0..=max).contains(&n) || n.fract() != 0.0 {
            return Err(harness("export: bad curve-on-surface head"));
        }
        Ok(n as u64)
    };
    let ty = integer(num, 1e6)? as i64;
    if ty != CRV_ON_SURF_TYPE {
        return Err(harness("export: -4 record is not an SmCrvOnSurf"));
    }
    let version = integer(num, 1e4)?;
    let dim = integer(num, 3.0)?;
    let n = integer(num, 1e8)? as usize;
    if n > remaining.saturating_sub(4) {
        return Err(harness("export too short for curve-on-surface record"));
    }
    let mut bytes = Vec::new();
    for _ in 0..n {
        bytes.push(integer(num, 255.0)? as u8);
    }
    let mut db = DatabaseReader::with_stream(FileType::Binary, InputStream::new(bytes));
    db.set_version(version);
    let k = Kernel::default();
    let mut curve = SmCurve::read_from_db(&k, ty, &mut db, dim, version)
        .map_err(|e| harness(format!("export: curve-on-surface record: {e}")))?;
    if db.in_stream().is_none_or(|s| s.remaining() != 0) {
        return Err(harness("export: curve-on-surface record not read whole"));
    }
    let CurveKind::CrvOnSurf(cos) = curve.kind_mut() else {
        return Err(harness("export: curve-on-surface reader changed its kind"));
    };
    if integer(num, 1.0)? == 0 {
        return Ok(curve);
    }
    // the cached approximation
    let matches = integer(num, 1.0)? != 0;
    let max_gap = num()?;
    let mut bits = vec![integer(num, 3.0)?];
    let n = integer(num, 1e7)?;
    let p = integer(num, 64.0)?;
    let m = integer(num, 1e7)?;
    let needed = (n + 1) * 4 + m + 1;
    if needed as usize > remaining {
        return Err(harness(
            "export too short for curve-on-surface approximation",
        ));
    }
    bits.extend([1, n, p, m]);
    for _ in 0..needed {
        bits.push(num()?.to_bits());
    }
    let rebuilt = cos
        .get_or_create_bs_approx(&k, k.tolerances.approx_tol_3d, None, Some(matches))
        .map_err(|e| harness(format!("export: curve-on-surface approximation: {e}")))?;
    if approx_bits(&rebuilt) != bits || cos.bs_approx_max_gap_3d().to_bits() != max_gap.to_bits() {
        return Err(harness(
            "export: curve-on-surface approximation not rebuilt bit for bit",
        ));
    }
    Ok(curve)
}

/// The dimension and raw NURB of a B-spline (`dim 1 n p m (x y z w)* knot*`, as bits).
fn approx_bits(c: &BSplineCurve) -> Vec<u64> {
    let mut bits = vec![c.dim()];
    match c.gw_nurb() {
        None => bits.push(0),
        Some(n) => {
            bits.extend([1, n.n() as u64, n.degree() as u64, n.m() as u64]);
            for p in n.cpts() {
                bits.extend([p.x, p.y, p.z, p.w].map(f64::to_bits));
            }
            bits.extend(n.knots().iter().map(|v| v.to_bits()));
        }
    }
    bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blocks() {
        let t: TestBrepIn =
            serde_json::from_str(r#"{"kind":"loop_orientation","case":3}"#).unwrap();
        assert_eq!(block(&t, false, true), vec![15.0, 3.0, 0.0, 1.0]);
        let t: TestBrepIn = serde_json::from_str(r#"{"kind":"xy_sheet","size":20}"#).unwrap();
        assert_eq!(block(&t, false, false), vec![17.0, 20.0, 0.0, 0.0]);
    }

    #[test]
    fn rejects_other_types() {
        let v = [1.0, 1.0, 3.0, 0.0];
        let mut it = v.iter().copied();
        let mut num = || it.next().ok_or_else(|| harness("short"));
        assert!(read_crv_on_surf(&mut num, v.len()).is_err());
    }
}
