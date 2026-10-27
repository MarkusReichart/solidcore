//! Face-merge cases of the part-5 runner: a classification curve classified on a face
//! (`SmFace::CurveOnClassify`) and merged into it (`SmCurveClassification::MergeIntervals`),
//! as prog_test `my_test_separating_closed_edges` (`tests/prog_test/src/smTopoTest.cpp`,
//! v1.0.1) does (shim `harness/cpp/shim/intersect_curve_class_p5_face.cpp`, layouts in
//! `solidcore_ffi::units::intersect_curve_class_p5`).
//!
//! The C++ side builds the source B-rep (a shared-topology-fixture block such as
//! `sm_make_xy_sheet`, or `SmPrimitiveCreation::CreateCone`), enables its editing and runs
//! the classification and the merge. The Rust side imports the shim's export of the same
//! B-rep (unit `topology.vertexuse`'s importer) and runs [`curve_on_classify_mut`] and
//! [`CurveClassification::merge_intervals_mut`] on it over the face's own B-rep
//! ([`BrepSeamObjects`]), with the classification curve as a temporary handle. Both sides
//! compare the classification after `CurveOnClassify` (intervals and the start, mid and
//! end classes and objects), the merge status, the new face, edge and vertex counts, the
//! B-rep's counts after it and whether the B-rep is valid (C++: `ValidatePointers` and
//! `AssertValid`; Rust: [`solidcore::topology::integrity::check`]).

use serde::{Deserialize, Serialize};

use solidcore::curve::crv_on_surf::CrvOnSurf;
use solidcore::curve::curve::Curve;
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::curve::line::Line;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::object::ObjectRef;
use solidcore::heal::face_props::LoopObject;
use solidcore::intersect::curve_class::p4::ClassificationInit;
use solidcore::intersect::curve_class::p5::ClassGeometryMut;
use solidcore::intersect::curve_class::p6::{
    MergeContext, MergeContextWithIdentity, MergeIntervalBrep,
};
use solidcore::intersect::curve_class::CurveClassification;
use solidcore::intersect::point_class::{
    ClassGeometry, ClassObject, ClassObjectKind, PointClassification,
};
use solidcore::kernel::Kernel;
use solidcore::math::{Point3, Vec3};
use solidcore::nurbs::b_spline_surface::SurfParamType;
use solidcore::query::face::p2::{
    curve_on_classify_mut, CurveOnClassifyBrep, CurveOnClassifyObjects, FaceBrep,
};
use solidcore::surface::surface::Surface;
use solidcore::surface::volume::Volume;
use solidcore::topology::face::p4::BrepSeamObjects;
use solidcore::topology::{Brep, EdgeId, EdgeuseId, FaceId, LoopId, VertexId, VertexuseId};
use solidcore::{EntityRef, SolidError};
use solidcore_ffi::units::intersect_curve_class_p5::{self as ffi, FaceOp};

use crate::units::topology_vertexuse::{import, sel};
use crate::util::test_brep::{self, TestBrepIn};
use crate::{HarnessError, Output};

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

fn yes() -> bool {
    true
}

/// `SmPrimitiveCreation::CreateCone(height, base_radius, top_radius, start_deg, end_deg,
/// SmAxis2Placement())` into a new B-rep.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(super) struct ConeIn {
    height: f64,
    base_radius: f64,
    top_radius: f64,
    start_deg: f64,
    end_deg: f64,
}

/// The classification curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub(super) enum CurveIn {
    /// An `SmCircle` (dimension 3).
    Circle {
        center: [f64; 3],
        x_axis: [f64; 3],
        y_axis: [f64; 3],
        start_deg: f64,
        end_deg: f64,
        radius: f64,
    },
    /// An `SmCrvOnSurf` (owner flag 2) on the face surface over the UV line between the
    /// fractions `(u0, v0)` and `(u1, v1)` of its natural UV domain.
    UvLine { u0: f64, v0: f64, u1: f64, v1: f64 },
}

/// A face selector by geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(super) enum FaceBy {
    /// The first face whose surface is periodic in U.
    UPeriodic,
}

/// A face-merge case: exactly one of `test_brep` and `cone`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(super) struct FaceMergeIn {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    test_brep: Option<TestBrepIn>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    cone: Option<ConeIn>,
    /// The face, modulo the face count in export order (unused with `face_by`).
    #[serde(default)]
    face: i64,
    /// `u_periodic`: the test's pick, the face whose surface is periodic in U over its
    /// natural domain (`IsPeriodic(GetNaturalUVDomain(), SM_SP_U)`); none or several is a
    /// harness error.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    face_by: Option<FaceBy>,
    curve: CurveIn,
    #[serde(default = "yes")]
    do_point_classify: bool,
    #[serde(default)]
    has_big_gaps: bool,
    #[serde(default = "yes")]
    check_3d: bool,
    #[serde(default = "yes")]
    uv_space_okay: bool,
}

impl FaceMergeIn {
    /// The source block of the shim input.
    fn source(&self) -> Result<Vec<f64>, HarnessError> {
        match (&self.test_brep, &self.cone) {
            (Some(t), None) => {
                let mut v = vec![0.0];
                v.extend(test_brep::block(t, false, false));
                Ok(v)
            }
            (None, Some(c)) => Ok(vec![
                1.0,
                c.height,
                c.base_radius,
                c.top_radius,
                c.start_deg,
                c.end_deg,
            ]),
            _ => Err(harness(
                "face_merge: give exactly one of test_brep and cone",
            )),
        }
    }

    /// The op arguments after the source.
    fn args(&self) -> Vec<f64> {
        let mut v = vec![
            f64::from(u8::from(self.face_by.is_some())),
            self.face as f64,
        ];
        match &self.curve {
            CurveIn::Circle {
                center,
                x_axis,
                y_axis,
                start_deg,
                end_deg,
                radius,
            } => {
                v.push(0.0);
                v.extend(center);
                v.extend(x_axis);
                v.extend(y_axis);
                v.extend([*start_deg, *end_deg, *radius]);
            }
            CurveIn::UvLine { u0, v0, u1, v1 } => v.extend([1.0, *u0, *v0, *u1, *v1]),
        }
        for f in [
            self.do_point_classify,
            self.has_big_gaps,
            self.check_3d,
            self.uv_space_okay,
        ] {
            v.push(f64::from(u8::from(f)));
        }
        v
    }
}

/// The shim's output (after a successful build) as the case output; a failed build is
/// reported identically by both sides.
fn shim_output(o: &[f64]) -> Result<Output, HarnessError> {
    match o.get(1..) {
        Some([-4.0]) => return Err(harness("face_by u_periodic: no U-periodic face")),
        Some([-5.0]) => return Err(harness("face_by u_periodic: several U-periodic faces")),
        _ => {}
    }
    match o.first() {
        Some(&0.0) => Ok(Output::new().values("out", o.get(1..).unwrap_or(&[]).to_vec())),
        Some(&-1.0) => Ok(Output::err(format!(
            "build_failed: {}",
            o.get(1).copied().unwrap_or(-9.0)
        ))),
        _ => Err(harness("unexpected c++ output")),
    }
}

/// The C++ side.
pub(super) fn cpp(i: &FaceMergeIn) -> Result<Output, HarnessError> {
    let mut input = i.source()?;
    input.extend(i.args());
    let o = ffi::face(FaceOp::Merge, "", &input).map_err(|e| harness(format!("c++: {e}")))?;
    shim_output(&o)
}

// ---- the Rust side ----------------------------------------------------------------------

/// The classification handles: the face's topology, its surface, edge and trim curves,
/// and temporaries (the classification curve is temporary 0).
#[derive(Debug, Clone, Copy, PartialEq)]
enum H {
    Topo(EntityRef),
    Surface(FaceId),
    EdgeCurve(EdgeId),
    Trim(EdgeuseId),
    Temp(usize),
}

/// The classification curve's handle.
const CURVE: H = H::Temp(0);

impl ClassObject for H {
    fn class_kind(&self) -> ClassObjectKind {
        match self {
            Self::Topo(e) => e.class_kind(),
            Self::Surface(_) => ClassObjectKind::Surface,
            Self::EdgeCurve(_) | Self::Trim(_) | Self::Temp(_) => ClassObjectKind::Curve,
        }
    }

    fn temporary(index: usize, kind: ClassObjectKind) -> Option<Self> {
        (kind == ClassObjectKind::Curve).then_some(Self::Temp(index))
    }
}

impl LoopObject for H {
    fn from_loop(l: LoopId) -> Self {
        Self::Topo(EntityRef::Loop(l))
    }

    fn as_loop(&self) -> Option<LoopId> {
        match self {
            Self::Topo(e) => e.as_loop(),
            _ => None,
        }
    }

    fn as_vertex(&self) -> Option<VertexId> {
        match self {
            Self::Topo(e) => e.as_vertex(),
            _ => None,
        }
    }

    fn from_face(f: FaceId) -> Option<Self> {
        Some(Self::Topo(EntityRef::Face(f)))
    }

    fn from_face_surface(f: FaceId) -> Option<Self> {
        Some(Self::Surface(f))
    }

    fn from_topology(t: EntityRef) -> Option<Self> {
        Some(Self::Topo(t))
    }

    fn as_topology(&self) -> Option<EntityRef> {
        match self {
            Self::Topo(t) => Some(*t),
            _ => None,
        }
    }

    fn from_edge_curve(e: EdgeId) -> Option<Self> {
        Some(Self::EdgeCurve(e))
    }

    fn as_edge_curve(&self) -> Option<EdgeId> {
        match self {
            Self::EdgeCurve(e) => Some(*e),
            _ => None,
        }
    }

    fn from_trim_curve(eu: EdgeuseId) -> Option<Self> {
        Some(Self::Trim(eu))
    }

    fn as_trim_curve(&self) -> Option<EdgeuseId> {
        match self {
            Self::Trim(eu) => Some(*eu),
            _ => None,
        }
    }
}

/// The face's own objects ([`BrepSeamObjects`]) with the classification curve as
/// temporary 0 (as `topology::face::p4`'s seam wrapper adds the seam curves).
struct WithCurve<'a, R> {
    inner: R,
    curve: &'a Curve,
}

impl<R: core::ops::Deref> ClassGeometry<H> for WithCurve<'_, R>
where
    R::Target: ClassGeometry<H>,
{
    fn curve(&self, o: H) -> Option<&Curve> {
        if o == CURVE {
            return Some(self.curve);
        }
        self.inner.curve(o)
    }
    fn surface(&self, o: H) -> Option<&Surface> {
        self.inner.surface(o)
    }
    fn volume(&self, o: H) -> Option<&Volume> {
        self.inner.volume(o)
    }
    fn point(&self, o: H) -> Option<Point3> {
        self.inner.point(o)
    }
    fn classification_vertex_point(&self, vertex: H) -> Option<Point3> {
        self.inner.classification_vertex_point(vertex)
    }
    fn topology_zone_tol_3d(&self, object: H) -> Option<f64> {
        self.inner.topology_zone_tol_3d(object)
    }
    fn vertex_face_gap(&self, vertex: H, surface: &Surface, face: Option<H>) -> Option<f64> {
        self.inner.vertex_face_gap(vertex, surface, face)
    }
    fn topology_entity(&self, object: H) -> Option<(&Brep, EntityRef)> {
        self.inner.topology_entity(object)
    }
    fn topology_handle(&self, entity: EntityRef) -> Option<H> {
        self.inner.topology_handle(entity)
    }
    fn topology_handle_in(&self, brep: &Brep, entity: EntityRef) -> Option<H> {
        self.inner.topology_handle_in(brep, entity)
    }
    fn owner_of(&self, owner: ObjectRef) -> Option<ObjectRef> {
        self.inner.owner_of(owner)
    }
}

impl<R: core::ops::Deref> CurveOnClassifyObjects<H> for WithCurve<'_, R>
where
    R::Target: CurveOnClassifyObjects<H>,
{
    fn face(&self, f: FaceId) -> H {
        self.inner.face(f)
    }
    fn vertex(&self, v: VertexId) -> H {
        self.inner.vertex(v)
    }
    fn edge(&self, e: EdgeId) -> H {
        self.inner.edge(e)
    }
    fn vertexuse(&self, vu: VertexuseId) -> H {
        self.inner.vertexuse(vu)
    }
    fn edgeuse(&self, eu: EdgeuseId) -> H {
        self.inner.edgeuse(eu)
    }
    fn edge_curve(&self, e: EdgeId) -> H {
        self.inner.edge_curve(e)
    }
    fn trim_curve(&self, eu: EdgeuseId) -> H {
        self.inner.trim_curve(eu)
    }
    fn face_surface(&self, f: FaceId) -> Option<H> {
        self.inner.face_surface(f)
    }
    fn loop_id(&self, tgt_loop: H) -> Option<LoopId> {
        self.inner.loop_id(tgt_loop)
    }
}

impl<R: core::ops::DerefMut> ClassGeometryMut<H> for WithCurve<'_, R>
where
    R::Target: ClassGeometryMut<H>,
{
    fn curve_mut(&mut self, o: H) -> Option<&mut Curve> {
        if o == CURVE {
            // the caller's curve (upstream's classification does not own it)
            return None;
        }
        self.inner.curve_mut(o)
    }
    fn topology_entity_mut(&mut self, o: H) -> Option<(&mut Brep, EntityRef)> {
        self.inner.topology_entity_mut(o)
    }
}

impl<R: core::ops::Deref> FaceBrep for WithCurve<'_, R>
where
    R::Target: FaceBrep,
{
    fn brep(&self) -> &Brep {
        self.inner.brep()
    }
}

impl<R: core::ops::DerefMut> CurveOnClassifyBrep<H> for WithCurve<'_, R>
where
    R::Target: CurveOnClassifyBrep<H>,
{
    fn brep_mut(&mut self) -> &mut Brep {
        self.inner.brep_mut()
    }
    fn classify_surface(&self, f: FaceId) -> solidcore::Result<Option<&Surface>> {
        self.inner.classify_surface(f)
    }
    fn keeps_made_trim(&self, eu: EdgeuseId) -> solidcore::Result<bool> {
        self.inner.keeps_made_trim(eu)
    }
}

impl<R: core::ops::DerefMut> MergeIntervalBrep<H> for WithCurve<'_, R>
where
    R::Target: MergeIntervalBrep<H>,
{
    fn merge_context(&mut self, object: H) -> Option<MergeContext<'_, H>> {
        self.inner.merge_context(object)
    }
    fn merge_context_with_curve_identity(
        &mut self,
        object: H,
    ) -> Option<MergeContextWithIdentity<'_, H>> {
        self.inner.merge_context_with_curve_identity(object)
    }
}

/// A Rust error as the shim's status (`Unported` defers the case).
enum Stop {
    Status(f64),
    Output(Box<Output>),
}

fn stop(e: &SolidError) -> Stop {
    if let Some(o) = Output::unported(e) {
        return Stop::Output(Box::new(o));
    }
    Stop::Status(match e {
        SolidError::Upstream { code, .. } => f64::from(*code),
        SolidError::InvalidInput { .. } => 1007.0,
        _ => -9.0,
    })
}

fn put_object(o: &mut Vec<f64>, h: Option<H>) {
    o.extend(match h {
        None => [-1.0, -1.0],
        Some(H::Temp(0)) => [100.0, -1.0],
        Some(H::Surface(_)) => [101.0, -1.0],
        Some(H::Topo(EntityRef::Face(f))) => [3.0, f64::from(f.index())],
        Some(H::Topo(EntityRef::Edge(e))) => [7.0, f64::from(e.index())],
        Some(H::Topo(EntityRef::Vertex(v))) => [9.0, f64::from(v.index())],
        Some(_) => [-2.0, -1.0],
    });
}

fn put_pc(o: &mut Vec<f64>, pc: &PointClassification<H>) {
    o.push(f64::from(pc.point_class().code()));
    put_object(o, pc.object());
}

/// The Rust side.
pub(super) fn rust(i: &FaceMergeIn) -> Result<Output, HarnessError> {
    let source = i.source()?;
    let ex =
        ffi::face(FaceOp::Export, "", &source).map_err(|e| harness(format!("c++ export: {e}")))?;
    if ex.first() != Some(&0.0) {
        // a failed build: the C++ side reports the same
        return shim_output(&ex);
    }
    if ex.get(1) == Some(&-3.0) {
        return Ok(Output::new().values("out", vec![-3.0]));
    }
    let (mut b, _exact) = import(&ex)?;
    b.set_editing_enabled(true);
    let n_faces = solidcore::topology::brep::p1::faces(&b, None)
        .map_err(|e| harness(format!("rust import: {e}")))?
        .len();
    let mut k = Kernel::default();
    let f = match i.face_by {
        None => FaceId::from_index(
            u32::try_from(sel(i.face, n_faces)).map_err(|_| harness("bad face selector"))?,
        ),
        Some(FaceBy::UPeriodic) => {
            // the test's pick: the face whose surface is periodic in U (exactly one, else a
            // harness error, as the C++ side)
            let mut found = None;
            let mut matches = 0usize;
            for index in 0..n_faces {
                let id =
                    FaceId::from_index(u32::try_from(index).map_err(|_| harness("face index"))?);
                let Some(s) = b
                    .face(id)
                    .map_err(|e| harness(e.to_string()))?
                    .surface
                    .as_ref()
                else {
                    continue;
                };
                let d = s
                    .natural_uv_domain_with_nurb(&k)
                    .map_err(|e| harness(format!("domain: {e}")))?;
                if s.is_periodic(&k, &d, SurfParamType::U)
                    .map_err(|e| harness(format!("is_periodic: {e}")))?
                {
                    found = found.or(Some(id));
                    matches += 1;
                }
            }
            match (found, matches) {
                (Some(id), 1) => id,
                (_, 0) => return Err(harness("face_by u_periodic: no U-periodic face")),
                _ => return Err(harness("face_by u_periodic: several U-periodic faces")),
            }
        }
    };

    // the classification curve
    let curve = match &i.curve {
        CurveIn::Circle {
            center,
            x_axis,
            y_axis,
            start_deg,
            end_deg,
            radius,
        } => Curve::from(
            Ellipse::new_circle(
                &k,
                Point3::new(center[0], center[1], center[2]),
                Vec3::new(x_axis[0], x_axis[1], x_axis[2]),
                Vec3::new(y_axis[0], y_axis[1], y_axis[2]),
                Extent1d::from_min_max(*start_deg, *end_deg),
                *radius,
                &EllipseOptions::default(),
            )
            .map_err(|e| harness(format!("circle: {e}")))?,
        ),
        CurveIn::UvLine { u0, v0, u1, v1 } => {
            let surface = b
                .face(f)
                .map_err(|e| harness(e.to_string()))?
                .surface
                .clone()
                .ok_or_else(|| harness("face without surface"))?;
            let d = surface
                .natural_uv_domain_with_nurb(&k)
                .map_err(|e| harness(format!("domain: {e}")))?;
            let (lo, hi) = (d.min(), d.max());
            let p0 = Point3::new(lo.x + u0 * (hi.x - lo.x), lo.y + v0 * (hi.y - lo.y), 0.0);
            let p1 = Point3::new(lo.x + u1 * (hi.x - lo.x), lo.y + v1 * (hi.y - lo.y), 0.0);
            let uv = Curve::from(Line::from_points(p0, p1, 2));
            Curve::from(
                CrvOnSurf::new(&k, uv, surface, None, 2)
                    .map_err(|e| harness(format!("crv on surf: {e}")))?,
            )
        }
    };

    let mut out = Vec::new();
    let tol = b.tolerance();
    let init = ClassificationInit {
        curve: Some(CURVE),
        interval: curve.natural_interval(),
        uv_curve: None,
        src_zone_tol_3d: tol,
        opt_data_size: 0,
        is_curve_owned: false,
        is_curve_current: true,
        surface1: Some(H::Surface(f)),
        surface2: None,
        uv_curve2: None,
        ti: None,
    };
    let mut cc = match CurveClassification::with_curve(&k, &init) {
        Ok(cc) => cc,
        Err(e) => return Err(harness(format!("classification: {e}"))),
    };

    // CurveOnClassify
    let classified = {
        let mut inner = BrepSeamObjects::<H, &mut Brep>::new(&k, &mut b, f)
            .map_err(|e| harness(format!("objects: {e}")))?;
        let mut objects = WithCurve {
            inner: &mut inner,
            curve: &curve,
        };
        curve_on_classify_mut(
            &k,
            f,
            &mut objects,
            i.do_point_classify,
            &mut cc,
            i.has_big_gaps,
            i.check_3d,
            i.uv_space_okay,
        )
    };
    let status = match classified {
        Ok(_) => 0.0,
        Err(e) => match stop(&e) {
            Stop::Status(s) => s,
            Stop::Output(o) => return Ok(*o),
        },
    };
    out.push(status);
    out.push(cc.size() as f64);
    for ivl in cc.intervals() {
        out.push(ivl.interval.min());
        out.push(ivl.interval.max());
        put_pc(&mut out, &ivl.start);
        put_pc(&mut out, &ivl.mid);
        put_pc(&mut out, &ivl.end);
    }
    // the test SERs CurveOnClassify before merging
    if status != 0.0 {
        return Ok(Output::new().values("out", out));
    }

    // MergeIntervals
    let merged = {
        let mut inner = BrepSeamObjects::<H, &mut Brep>::for_merge(&mut b, f)
            .map_err(|e| harness(format!("merge objects: {e}")))?;
        let mut objects = WithCurve {
            inner: &mut inner,
            curve: &curve,
        };
        cc.merge_intervals_mut(&mut k, &mut objects, false)
    };
    let counts = match merged {
        Ok(m) => [m.faces.len(), m.edges.len(), m.vertices.len()],
        Err(e) => match stop(&e) {
            // the test SERs MergeIntervals too
            Stop::Status(s) => {
                out.push(s);
                return Ok(Output::new().values("out", out));
            }
            Stop::Output(o) => return Ok(*o),
        },
    };
    out.push(0.0);
    out.extend(counts.map(|n| n as f64));
    let count =
        |n: solidcore::Result<usize>| n.map(|n| n as f64).map_err(|e| harness(e.to_string()));
    out.push(count(
        solidcore::topology::brep::p1::faces(&b, None).map(|v| v.len()),
    )?);
    out.push(count(
        solidcore::topology::brep::p1::edges(&b, None).map(|v| v.len()),
    )?);
    out.push(count(
        solidcore::topology::brep::p1::vertices(&b, None).map(|v| v.len()),
    )?);
    out.push(f64::from(u8::from(
        solidcore::topology::integrity::check(&b).is_ok(),
    )));
    Ok(Output::new().values("out", out))
}
