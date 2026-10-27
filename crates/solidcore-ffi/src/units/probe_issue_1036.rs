//! FFI wrapper for the one-off tangency-labs/solidcore#1036 probe (not a unit).
//!
//! Shim: `harness/cpp/shim/probe_issue_1036.cpp` (`sc_probe_issue_1036_run`). Upstream:
//! `source/SM_API/src/SmApiPrimitives.cpp` (`SmApiCreateSkinPrimitive`),
//! `source/SM_API/src/SmApiBrep.cpp` (`SmApiBooleanDifference`, `SmApiBooleanIntersection`,
//! `SmApiBooleanUnion`), `source/SM_API/src/SmApiQueries.cpp` (`SmApiBrepComputeVolume`).
//!
//! The program is `op skinA skinB`. A skin is `nc`, `nc` curves, `degree cap` (one profile
//! curve per section, capped). A curve is `dim n p m`, `n+1` x `x y z w`, `m+1` knots
//! (`SmBSplineCurve(gw_CURVE*, dim, FALSE)`). The output is, per skin,
//! `status [faces edges vertices volume_status volume]` (bracketed only on success), then
//! for a Boolean op its `status [has faces edges vertices volume_status volume]`.
//!
//! [`run`] accepts only clamped, finite, non-rational 3D NURB curves (`m = n+p+1`,
//! `1 <= p <= min(n, 17)`, non-decreasing knots, `w` = `NOW`), 2 to 64 sections per skin
//! and a loft degree of 1 to 17.

use crate::validate::{bad, NOW, NOZ};
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_probe_issue_1036_run(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// What to run on the two skins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Build and measure both skins only.
    Measure = 0,
    /// `SmApiBooleanDifference(A, B)`.
    Difference = 1,
    /// `SmApiBooleanIntersection(A, B)`.
    Intersection = 2,
    /// `SmApiBooleanUnion(A, B)`.
    Union = 3,
}

/// A non-rational 3D NURB section curve: degree, weighted control points
/// (`[x, y, z, w]`, `w` = `NOW`) and the full knot vector.
#[derive(Debug, Clone, PartialEq)]
pub struct Nurb {
    /// Degree `p`.
    pub p: usize,
    /// Control points.
    pub cpts: Vec<[f64; 4]>,
    /// Knots (`cpts.len() + p + 1` of them).
    pub knots: Vec<f64>,
}

/// A skin solid: its sections, loft degree and whether its ends are capped.
#[derive(Debug, Clone, PartialEq)]
pub struct Skin {
    /// One curve per section.
    pub sections: Vec<Nurb>,
    /// Loft degree.
    pub degree: usize,
    /// Cap the first and last sections.
    pub cap: bool,
}

fn encode_nurb(v: &mut Vec<f64>, c: &Nurb) -> Result<(), FfiError> {
    let n1 = c.cpts.len();
    if n1 < 2 || c.p < 1 || c.p > 17 || c.p > n1 - 1 || c.knots.len() != n1 + c.p + 1 {
        return Err(bad("probe_issue_1036: malformed NURB sizes"));
    }
    if c.knots.iter().any(|k| !k.is_finite()) || c.knots.windows(2).any(|w| w[0] > w[1]) {
        return Err(bad("probe_issue_1036: knots not finite and non-decreasing"));
    }
    if c.cpts.iter().any(|p| {
        !(p[0].is_finite() && p[1].is_finite() && p[2].is_finite()) || p[2] == NOZ || p[3] != NOW
    }) {
        return Err(bad(
            "probe_issue_1036: control points must be finite, 3D, non-rational",
        ));
    }
    v.extend([3.0, (n1 - 1) as f64, c.p as f64, (c.knots.len() - 1) as f64]);
    for p in &c.cpts {
        v.extend(p);
    }
    v.extend(&c.knots);
    Ok(())
}

fn encode_skin(v: &mut Vec<f64>, s: &Skin) -> Result<(), FfiError> {
    if !(2..=64).contains(&s.sections.len()) || !(1..=17).contains(&s.degree) {
        return Err(bad("probe_issue_1036: 2 to 64 sections, degree 1 to 17"));
    }
    v.push(s.sections.len() as f64);
    for c in &s.sections {
        encode_nurb(v, c)?;
    }
    v.extend([s.degree as f64, f64::from(u8::from(s.cap))]);
    Ok(())
}

/// Runs `op` on skins `a` and `b` in the C++ kernel; the raw output record (module docs).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for inputs outside the accepted range; the shim's status.
pub fn run(op: Op, a: &Skin, b: &Skin) -> Result<Vec<f64>, FfiError> {
    let mut prog = vec![op as i32 as f64];
    encode_skin(&mut prog, a)?;
    encode_skin(&mut prog, b)?;
    let n_in = i64::try_from(prog.len()).map_err(|_| bad("program too long"))?;
    call_with_buffer(64, |buf: &mut [f64], len| {
        // SAFETY: `prog` is valid for `n_in` reads; `buf` for `buf.len()` writes; `len`
        // is a valid i64; nothing is retained.
        unsafe {
            sys::sc_probe_issue_1036_run(
                prog.as_ptr(),
                n_in,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}
