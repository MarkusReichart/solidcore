// SPDX-License-Identifier: Apache-2.0
// Diagnostic shim for the `corpus.parity` runner: per-face measures of one upstream
// test model (SmBrep::ReadFromFile, B-rep files only), to locate a whole-B-rep area or
// volume disagreement on a face. Not a unit runner: its probe test is
// crates/difftest/tests/probe_sphere_diff_faces.rs.
//
// The objects read are leaked on purpose (the probe runs once per process).
#include "sc_shim.h"
#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <SmBrep.h>
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmExtent3d.h>
#include <SmEdge.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmGap.h>
#include <SmLoop.h>
#include <SmLoopuse.h>
#include <SmSurface.h>
#include <SmVertex.h>
#include <SmVertexuse.h>
#include <cmath>
#include <limits>
#include <vector>

namespace {
const SmContext& faces_context()
{
    static SmContext* p = new SmContext();
    return *p;
}
} // namespace

// Output: `read_status n`, then per face of `GetFaces` (in its order) 8 values:
// surface type (`GetType`), area (`SmApiFaceComputeArea`'s integration without its
// accuracy clamp: face box midpoint origin, SM_PPF_AREA, made positive), the signed
// area and volume of `ComputePreciseProperties(SM_OT_SAME, accuracy, brep box mid,
// brep estimate, 0, .., SM_PPF_ALL)` (the origin and estimate of
// `SmBrep::ComputePreciseProperties`), the face box midpoint (x y z) and the status of
// the two integrations (0 or the first failing SmStatus).
SC_EXPORT int32_t sc_corpus_faces_measure(const char* path, double accuracy, double* out,
                                          int64_t cap, int64_t* out_len)
{
    if (path == nullptr || !(accuracy > 0.0) || !std::isfinite(accuracy))
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = faces_context();
    SmBrep* b = new (ctx) SmBrep();
    SmStatus st = b->ReadFromFile(ctx, path, SM_ASCII);
    std::vector<double> o;
    if (st != SM_SUCCESS)
    {
        o.push_back((double)st);
        o.push_back(0.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    SmExtent3d bbox;
    SmStatus bs = b->CalculateBoundingBox(bbox);
    SmVector3d bsize = bbox.GetSize();
    double est = smos_Max(bsize.z, smos_Max(bsize.x, bsize.y));
    est = est * est / 100.0;
    SmTArray<SmFace*> faces;
    b->GetFaces(faces);
    o.push_back(bs == SM_SUCCESS ? 0.0 : (double)bs);
    o.push_back((double)faces.GetSize());
    for (ULONG i = 0; i < faces.GetSize(); i++)
    {
        const SmFace* f = faces[i];
        SmStatus fs = SM_SUCCESS;
        SmExtent3d fbox;
        SmStatus s = f->CalculateBoundingBox(fbox);
        if (s != SM_SUCCESS)
            fs = s;
        SmVector3d fsize = fbox.GetSize();
        double d = smos_Max(fsize.z, smos_Max(fsize.x, fsize.y));
        double area = 0.0, dv = 0.0;
        SmTArray<SmVector3d> moments;
        s = f->ComputePreciseProperties(SM_OT_SAME, accuracy, fbox.GetMid(), d * d / 100.0, 0.0,
                                        area, dv, moments, SM_PPF_AREA);
        if (s != SM_SUCCESS && fs == SM_SUCCESS)
            fs = s;
        if (area < 0.0)
            area = -area;
        double a2 = 0.0, v2 = 0.0;
        SmTArray<SmVector3d> m2;
        s = f->ComputePreciseProperties(SM_OT_SAME, accuracy, bbox.GetMid(), est, 0.0, a2, v2,
                                        m2, SM_PPF_ALL);
        if (s != SM_SUCCESS && fs == SM_SUCCESS)
            fs = s;
        const SmSurface* srf = f->GetSurface();
        SmPoint3d mid = fbox.GetMid();
        o.insert(o.end(), {srf != nullptr ? (double)srf->GetType() : -1.0, area, a2, v2, mid.x,
                           mid.y, mid.z, (double)fs});
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}

// The trim boundary face `face_index` (of `GetFaces`) integrates over:
// `CreateOrientedTrimBoundary(ctx, TRUE, TRUE, ..)` as `SmFace::ComputePreciseProperties`
// calls it. Output: `read_status`, the loop count, each loop's curve count, then per
// curve `degree n t0 t1` and its `n` control points as `x y w` (cartesian, weight).
SC_EXPORT int32_t sc_corpus_faces_trims(const char* path, int64_t face_index, double* out,
                                        int64_t cap, int64_t* out_len)
{
    if (path == nullptr || face_index < 0)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = faces_context();
    SmBrep* b = new (ctx) SmBrep();
    SmStatus st = b->ReadFromFile(ctx, path, SM_ASCII);
    std::vector<double> o;
    SmTArray<SmFace*> faces;
    if (st == SM_SUCCESS)
        b->GetFaces(faces);
    if (st != SM_SUCCESS || (ULONG)face_index >= faces.GetSize())
    {
        o.push_back(st != SM_SUCCESS ? (double)st : -1.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    SmTArray<ULONG> loops;
    SmTArray<SmBSplineCurve*> curves;
    st = faces[(ULONG)face_index]->CreateOrientedTrimBoundary(ctx, TRUE, TRUE, loops, curves);
    o.push_back((double)st);
    o.push_back((double)loops.GetSize());
    for (ULONG i = 0; i < loops.GetSize(); i++)
        o.push_back((double)loops[i]);
    for (ULONG i = 0; i < curves.GetSize(); i++)
    {
        SmTArray<SmPoint3d> cps;
        SmTArray<double> w;
        curves[i]->GetControlPolygon(cps, w);
        SmExtent1d ivl = curves[i]->GetNaturalInterval();
        o.insert(o.end(), {(double)curves[i]->GetDegree(), (double)cps.GetSize(), ivl.GetMin(),
                           ivl.GetMax()});
        for (ULONG j = 0; j < cps.GetSize(); j++)
            o.insert(o.end(), {cps[j].x, cps[j].y, j < w.GetSize() ? w[j] : 1.0});
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}

// Each edge use of face `face_index`'s loops (`GetLoops`, each loop's `GetEdgeuses`)
// with `SmEdgeuse::CreateUVTrimCurve` called on it directly, before any
// `sm_FixUVTrimCurves` of `SmFace::CreateUVTrimCurves` (the face is read afresh and no
// trim is attached). Output: `read_status`, then per edge use `status degree n t0 t1
// e0 e1 x0 y0 xn yn` (the curve's interval, the edge's interval, its first and last
// control points; zeros when the call fails).
SC_EXPORT int32_t sc_corpus_faces_drops(const char* path, int64_t face_index, double* out,
                                        int64_t cap, int64_t* out_len)
{
    if (path == nullptr || face_index < 0)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = faces_context();
    SmBrep* b = new (ctx) SmBrep();
    SmStatus st = b->ReadFromFile(ctx, path, SM_ASCII);
    std::vector<double> o;
    SmTArray<SmFace*> faces;
    if (st == SM_SUCCESS)
        b->GetFaces(faces);
    if (st != SM_SUCCESS || (ULONG)face_index >= faces.GetSize())
    {
        o.push_back(st != SM_SUCCESS ? (double)st : -1.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    o.push_back(0.0);
    SmTArray<SmLoop*> loops;
    faces[(ULONG)face_index]->GetLoops(loops);
    for (ULONG l = 0; l < loops.GetSize(); l++)
    {
        SmTArray<SmEdgeuse*> eus;
        loops[l]->GetEdgeuses(eus);
        for (ULONG i = 0; i < eus.GetSize(); i++)
        {
            double dist = 0.0;
            SmBSplineCurve* c = nullptr;
            SmStatus cs = eus[i]->CreateUVTrimCurve(dist, c, nullptr, nullptr);
            SmExtent1d eivl = eus[i]->GetEdge()->GetInterval();
            if (cs != SM_SUCCESS || c == nullptr)
            {
                o.insert(o.end(), {(double)cs, 0.0, 0.0, 0.0, 0.0, eivl.GetMin(), eivl.GetMax(),
                                   0.0, 0.0, 0.0, 0.0});
                continue;
            }
            SmTArray<SmPoint3d> cps;
            SmTArray<double> w;
            c->GetControlPolygon(cps, w);
            SmExtent1d ivl = c->GetNaturalInterval();
            SmPoint3d p0 = cps.GetSize() ? cps[0] : SmPoint3d(0, 0, 0);
            SmPoint3d pn = cps.GetSize() ? cps[cps.GetSize() - 1] : SmPoint3d(0, 0, 0);
            o.insert(o.end(), {(double)cs, (double)c->GetDegree(), (double)cps.GetSize(),
                               ivl.GetMin(), ivl.GetMax(), eivl.GetMin(), eivl.GetMax(), p0.x,
                               p0.y, pn.x, pn.y});
        }
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}

// upstream's corner fix of SmFace::CreateUVTrimCurves (SmFace.cpp, not in a header)
SmStatus sm_FixUVTrimCurves(SmEdgeuse* pEUCurrent, SmEdgeuse* pEUPrevious, SmVertex& crVertex,
                            double dTol3d, double dUVDistBetween, const SmExtent2d& crUVDomain,
                            SmBoolean& rbFixedCurves);

namespace {
double trim_points(const SmEdgeuse* eu)
{
    SmBSplineCurve* c = eu->GetUVTrimCurve();
    if (c == nullptr)
        return -1.0;
    SmTArray<SmPoint3d> cps;
    SmTArray<double> w;
    c->GetControlPolygon(cps, w);
    return (double)cps.GetSize();
}
} // namespace

// The steps of SmFace::CreateUVTrimCurves on face `face_index`, replayed: every edge use
// of the upward loop uses gets its trim (GetOrCreateUVTrimCurve), then the corner loop
// (CCW gap, sm_FixUVTrimCurves when above the scaled zero; the gap-at-pole test is not
// replayed: `poles` reports the surface singularities so a reader can tell). Output:
// `read_status poles`, then per corner 20 values: loop, index, gap, called, the vertex
// use's UV (x y: VertexFaceGap's face UV, else ComputeUVPoint), the current edge use's
// end parameter tCurr, its DropPoint (success, parameter, distance) as
// sm_FixUVTrimCurves calls it, the control point counts of the current and previous
// trims before and after, whether the fix reports a change, and the previous edge
// use's tPrev and DropPoint (success, parameter, distance).
SC_EXPORT int32_t sc_corpus_faces_fixsteps(const char* path, int64_t face_index, double* out,
                                           int64_t cap, int64_t* out_len)
{
    if (path == nullptr || face_index < 0)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = faces_context();
    SmBrep* b = new (ctx) SmBrep();
    SmStatus st = b->ReadFromFile(ctx, path, SM_ASCII);
    std::vector<double> o;
    SmTArray<SmFace*> faces;
    if (st == SM_SUCCESS)
        b->GetFaces(faces);
    if (st != SM_SUCCESS || (ULONG)face_index >= faces.GetSize())
    {
        o.push_back(st != SM_SUCCESS ? (double)st : -1.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    SmFace* f = faces[(ULONG)face_index];
    SmExtent2d uv_domain = f->GetUVDomain();
    o.push_back(0.0);
    o.push_back((double)f->GetSurface()->GetSingularities());
    SmTArray<SmLoopuse*> lus;
    f->GetUpwardLoopuses(lus);
    for (ULONG l = 0; l < lus.GetSize(); l++)
    {
        SmTArray<SmEdgeuse*> eus;
        lus[l]->GetEdgeuses(eus);
        for (ULONG j = 0; j < eus.GetSize(); j++)
        {
            if (eus[j]->GetEdge() == nullptr || eus[j]->GetUVTrimCurve() != nullptr)
                continue;
            SmBSplineCurve* c = nullptr;
            double d = 0.0;
            eus[j]->GetOrCreateUVTrimCurve(c, &d);
        }
    }
    for (ULONG l = 0; l < lus.GetSize(); l++)
    {
        SmTArray<SmEdgeuse*> eus;
        lus[l]->GetEdgeuses(eus);
        for (ULONG j = 0; j < eus.GetSize(); j++)
        {
            SmEdgeuse* eu = eus[j];
            SmEdgeuse* cw = eu->GetCWEdgeuse();
            SmVertex* v = eu->GetVertexuse()->GetVertex();
            SmPoint3d p0, p1;
            eu->NormalizedEvaluate(0.0, TRUE, p0);
            cw->NormalizedEvaluate(1.0, TRUE, p1);
            double gap = p0.DistanceBetween(p1);
            double scaled = SM_EFF_ZERO * (1.0 + p0.GetMaxDimension());
            const double nan = std::numeric_limits<double>::quiet_NaN();
            if (!(gap > scaled))
            {
                o.insert(o.end(), {(double)l, (double)j, gap, 0.0, nan, nan, nan, nan, nan, nan,
                                   trim_points(eu), trim_points(cw), nan, nan, nan, nan, nan,
                                   nan, nan, nan});
                continue;
            }
            // the inputs sm_FixUVTrimCurves computes before its DropPoint
            SmBSplineCurve* cur = eu->GetUVTrimCurve();
            SmExtent1d ivl = cur->GetNaturalInterval();
            double t_curr = eu->GetOrientation() == SM_OT_SAME ? ivl.GetMin() : ivl.GetMax();
            SmPoint3d uv_curr;
            cur->EvaluatePoint(t_curr, uv_curr);
            SmPoint2d uv_vertex;
            SmVertexuse* vu = eu->GetVertexuse();
            SmVertexFaceGap* fg = vu->GetVertexFaceGap();
            if (fg)
                uv_vertex = fg->GetFaceUV();
            else
                vu->ComputeUVPoint(uv_vertex, TRUE);
            SmPoint3d uv_vertex3d(uv_vertex.x, uv_vertex.y, 0.0);
            double edge_gap = uv_vertex3d.DistanceBetween(uv_curr);
            SmBoolean ok = FALSE;
            double t_new = nan, dist = nan;
            cur->DropPoint(ivl, uv_vertex3d, NULL, 1.1 * (SM_EFF_ZERO + edge_gap), &t_curr, ok,
                           t_new, dist, SM_SO_INTERSECT);
            SmBSplineCurve* prev = cw->GetUVTrimCurve();
            SmExtent1d pivl = prev->GetNaturalInterval();
            double t_prev = cw->GetOrientation() == SM_OT_SAME ? pivl.GetMax() : pivl.GetMin();
            SmPoint3d uv_prev;
            prev->EvaluatePoint(t_prev, uv_prev);
            double prev_gap = uv_vertex3d.DistanceBetween(uv_prev);
            SmBoolean pok = FALSE;
            double pt_new = nan, pdist = nan;
            prev->DropPoint(pivl, uv_vertex3d, NULL, 1.1 * (SM_EFF_ZERO + prev_gap), &t_prev, pok,
                            pt_new, pdist, SM_SO_INTERSECT);
            double n_cur = trim_points(eu), n_cw = trim_points(cw);
            double tol = smos_4Max((double)b->GetTolerance(), (double)v->GetTolerance(),
                                   (double)eu->GetEdge()->GetTolerance(),
                                   (double)cw->GetEdge()->GetTolerance());
            SmBoolean fixed = FALSE;
            sm_FixUVTrimCurves(eu, cw, *v, tol, gap, uv_domain, fixed);
            o.insert(o.end(), {(double)l, (double)j, gap, 1.0, uv_vertex.x, uv_vertex.y, t_curr,
                               ok ? 1.0 : 0.0, t_new, dist, n_cur, n_cw, trim_points(eu),
                               trim_points(cw), fixed ? 1.0 : 0.0, 0.0, t_prev, pok ? 1.0 : 0.0,
                               pt_new, pdist});
        }
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}

namespace {
// `n nk knots..` of a B-spline: control point count and the full knot vector.
void push_shape(const SmBSplineCurve* c, std::vector<double>& o)
{
    SmTArray<SmPoint3d> cps;
    SmTArray<double> w;
    c->GetControlPolygon(cps, w);
    SmTArray<double> knots;
    c->GetKnotsAll(knots);
    o.push_back((double)cps.GetSize());
    o.push_back((double)knots.GetSize());
    for (ULONG i = 0; i < knots.GetSize(); i++)
        o.push_back(knots[i]);
}
} // namespace

// The trim step of sm_FixUVTrimCurves on one freshly made trim: edge use `eu_index` of
// face `face_index` (its upward loop uses' edge uses in order) gets
// `CreateUVTrimCurve`, then `Trim([t0, t1])` and `EditParameterization` back over the
// curve's interval, as sm_FixUVTrimCurves does. Output: `read_status`, the made curve
// `n nk knots..`, `Trim`'s status and its (possibly snapped) interval `t0 t1`, the
// shape after `Trim`, `EditParameterization`'s status and the shape after it.
SC_EXPORT int32_t sc_corpus_faces_trimstep(const char* path, int64_t face_index,
                                           int64_t eu_index, double t0, double t1,
                                           double* out, int64_t cap, int64_t* out_len)
{
    if (path == nullptr || face_index < 0 || eu_index < 0)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = faces_context();
    SmBrep* b = new (ctx) SmBrep();
    SmStatus st = b->ReadFromFile(ctx, path, SM_ASCII);
    std::vector<double> o;
    SmTArray<SmFace*> faces;
    if (st == SM_SUCCESS)
        b->GetFaces(faces);
    if (st != SM_SUCCESS || (ULONG)face_index >= faces.GetSize())
    {
        o.push_back(st != SM_SUCCESS ? (double)st : -1.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    SmTArray<SmLoopuse*> lus;
    faces[(ULONG)face_index]->GetUpwardLoopuses(lus);
    std::vector<SmEdgeuse*> all;
    for (ULONG l = 0; l < lus.GetSize(); l++)
    {
        SmTArray<SmEdgeuse*> eus;
        lus[l]->GetEdgeuses(eus);
        for (ULONG j = 0; j < eus.GetSize(); j++)
            all.push_back(eus[j]);
    }
    if ((size_t)eu_index >= all.size())
    {
        o.push_back(-2.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    double dist = 0.0;
    SmBSplineCurve* c = nullptr;
    st = all[(size_t)eu_index]->CreateUVTrimCurve(dist, c, nullptr, nullptr);
    if (st != SM_SUCCESS || c == nullptr)
    {
        o.push_back((double)st);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    o.push_back(0.0);
    push_shape(c, o);
    SmExtent1d full = c->GetNaturalInterval();
    SmExtent1d ivl;
    ivl.SetMinMax(t0, t1);
    SmStatus ts = c->Trim(ivl);
    o.push_back((double)ts);
    o.push_back(ivl.GetMin());
    o.push_back(ivl.GetMax());
    push_shape(c, o);
    SmStatus es = c->EditParameterization(full);
    o.push_back((double)es);
    push_shape(c, o);
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
