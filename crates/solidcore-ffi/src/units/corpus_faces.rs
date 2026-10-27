//! FFI wrapper for the diagnostic shim `harness/cpp/shim/corpus_faces.cpp`
//! (`sc_corpus_faces_measure`): per-face measures of one upstream test model read with
//! `SmBrep::ReadFromFile`, to locate a whole-B-rep disagreement of `corpus.parity` on a
//! face.
//!
//! | Function    | Output |
//! |-------------|--------|
//! | [`measure`] | `read_status n`, then per face `type area area_same volume_same mid_x mid_y mid_z status` |
//! | [`drops`]   | `status`, then per edge use of the face's loops `status degree n t0 t1 e0 e1 x0 y0 xn yn` (`CreateUVTrimCurve` called directly) |
//! | [`fixsteps`] | `status poles`, then per corner of the face's loops 20 values (`SmFace::CreateUVTrimCurves` replayed, see the shim) |
//! | [`trimstep`] | `status`, the made trim `n nk knots..`, `Trim`'s status and interval, the shape after it, `EditParameterization`'s status and the shape after it |
//! | [`trims`]   | `status loops counts.. `, then per trim curve `degree n t0 t1` and `n` control points `x y w` |
//!
//! See the shim for the definition of each value. The objects read are leaked: call it
//! from a short-lived process.

use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_corpus_faces_measure(
            path: *const std::ffi::c_char,
            accuracy: f64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_corpus_faces_drops(
            path: *const std::ffi::c_char,
            face_index: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_corpus_faces_fixsteps(
            path: *const std::ffi::c_char,
            face_index: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_corpus_faces_trimstep(
            path: *const std::ffi::c_char,
            face_index: i64,
            eu_index: i64,
            t0: f64,
            t1: f64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_corpus_faces_trims(
            path: *const std::ffi::c_char,
            face_index: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Values per face in [`measure`]'s output.
pub const PER_FACE: usize = 8;

/// Reads B-rep file `path` with upstream's reader and measures every face at
/// integration accuracy `accuracy` (see the module docs for the layout).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a path with a NUL byte or a non-positive or
/// non-finite accuracy; [`FfiError::Status`] for a shim failure.
pub fn measure(path: &str, accuracy: f64) -> Result<Vec<f64>, FfiError> {
    if !(accuracy.is_finite() && accuracy > 0.0) {
        return Err(FfiError::InvalidArgument("accuracy".into()));
    }
    let cpath = std::ffi::CString::new(path)
        .map_err(|_| FfiError::InvalidArgument("path has a NUL byte".into()))?;
    call_with_buffer(256, |buf, len| {
        // SAFETY: cpath is a valid NUL-terminated string; buf is valid for buf.len()
        // writes and the shim honours the capacity; len is valid.
        unsafe {
            sys::sc_corpus_faces_measure(
                cpath.as_ptr(),
                accuracy,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Reads B-rep file `path` with upstream's reader and returns the trim boundary face
/// `face_index` (of `GetFaces`) is integrated over (`CreateOrientedTrimBoundary(ctx,
/// TRUE, TRUE, ..)`; see the module docs for the layout).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a path with a NUL byte or a negative face index;
/// [`FfiError::Status`] for a shim failure.
pub fn trims(path: &str, face_index: i64) -> Result<Vec<f64>, FfiError> {
    if face_index < 0 {
        return Err(FfiError::InvalidArgument("face_index".into()));
    }
    let cpath = std::ffi::CString::new(path)
        .map_err(|_| FfiError::InvalidArgument("path has a NUL byte".into()))?;
    call_with_buffer(4096, |buf, len| {
        // SAFETY: cpath is a valid NUL-terminated string; buf is valid for buf.len()
        // writes and the shim honours the capacity; len is valid.
        unsafe {
            sys::sc_corpus_faces_trims(
                cpath.as_ptr(),
                face_index,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Reads B-rep file `path` with upstream's reader and calls `SmEdgeuse::CreateUVTrimCurve`
/// directly on every edge use of face `face_index`'s loops (no trim attached, no
/// `sm_FixUVTrimCurves`; see the module docs for the layout).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a path with a NUL byte or a negative face index;
/// [`FfiError::Status`] for a shim failure.
pub fn drops(path: &str, face_index: i64) -> Result<Vec<f64>, FfiError> {
    if face_index < 0 {
        return Err(FfiError::InvalidArgument("face_index".into()));
    }
    let cpath = std::ffi::CString::new(path)
        .map_err(|_| FfiError::InvalidArgument("path has a NUL byte".into()))?;
    call_with_buffer(512, |buf, len| {
        // SAFETY: cpath is a valid NUL-terminated string; buf is valid for buf.len()
        // writes and the shim honours the capacity; len is valid.
        unsafe {
            sys::sc_corpus_faces_drops(
                cpath.as_ptr(),
                face_index,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Reads B-rep file `path` with upstream's reader and replays `SmFace::CreateUVTrimCurves`
/// on face `face_index`: the trims made, then each corner with `sm_FixUVTrimCurves`'
/// inputs, its `DropPoint` and the trims' control point counts before and after (see
/// the shim for the layout).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a path with a NUL byte or a negative face index;
/// [`FfiError::Status`] for a shim failure.
pub fn fixsteps(path: &str, face_index: i64) -> Result<Vec<f64>, FfiError> {
    if face_index < 0 {
        return Err(FfiError::InvalidArgument("face_index".into()));
    }
    let cpath = std::ffi::CString::new(path)
        .map_err(|_| FfiError::InvalidArgument("path has a NUL byte".into()))?;
    call_with_buffer(512, |buf, len| {
        // SAFETY: cpath is a valid NUL-terminated string; buf is valid for buf.len()
        // writes and the shim honours the capacity; len is valid.
        unsafe {
            sys::sc_corpus_faces_fixsteps(
                cpath.as_ptr(),
                face_index,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Reads B-rep file `path` with upstream's reader, makes the trim of edge use
/// `eu_index` of face `face_index` (`CreateUVTrimCurve`) and applies `sm_FixUVTrimCurves`'
/// trim step to it: `Trim([t0, t1])`, then `EditParameterization` back over its interval
/// (see the module docs for the layout).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a path with a NUL byte, a negative index or a
/// non-finite parameter; [`FfiError::Status`] for a shim failure.
pub fn trimstep(
    path: &str,
    face_index: i64,
    eu_index: i64,
    t0: f64,
    t1: f64,
) -> Result<Vec<f64>, FfiError> {
    if face_index < 0 || eu_index < 0 || !t0.is_finite() || !t1.is_finite() {
        return Err(FfiError::InvalidArgument("face_index, eu_index, t0 or t1".into()));
    }
    let cpath = std::ffi::CString::new(path)
        .map_err(|_| FfiError::InvalidArgument("path has a NUL byte".into()))?;
    call_with_buffer(1024, |buf, len| {
        // SAFETY: cpath is a valid NUL-terminated string; buf is valid for buf.len()
        // writes and the shim honours the capacity; len is valid.
        unsafe {
            sys::sc_corpus_faces_trimstep(
                cpath.as_ptr(),
                face_index,
                eu_index,
                t0,
                t1,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}
