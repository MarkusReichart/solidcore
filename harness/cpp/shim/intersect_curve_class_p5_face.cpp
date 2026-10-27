// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `intersect.curve_class_p5`, face-merge cases: a classification
// curve classified on a face (SmFace::CurveOnClassify) and merged into it
// (SmCurveClassification::MergeIntervals), as prog_test my_test_separating_closed_edges
// (tests/prog_test/src/smTopoTest.cpp, v1.0.1) does.
//
// One entry point, sc_intersect_curve_class_p5_face(op, path, in, n_in, out, cap, out_len).
// Layouts are documented in crates/solidcore-ffi/src/units/intersect_curve_class_p5.rs,
// which validates every input. Every op starts with a source:
//   `0 <B-rep block>`               a B-rep block of the shared topology fixture
//                                   (sc_topology_fixture.h; e.g. block 12, sm_make_xy_sheet)
//   `1 height base top a0 a1`       SmPrimitiveCreation(infinite region).CreateCone(height,
//                                   base, top, a0, a1, SmAxis2Placement()) into a new B-rep
// then the B-rep's editing is enabled (the test sets m_bEditingEnabled for both cases) and
// its non-exported analytic geometry replaced by the plain B-splines the Rust import
// rebuilds (as topology_face_p4.cpp's normalise). Op 0 exports it (the fixture's layout,
// imported by unit topology.vertexuse's importer); op 1 runs the merge.
#include "sc_topology_fixture.h"
#include <SmAssertArray.h>
#include <SmCircle.h>
#include <SmCrvOnSurf.h>
#include <SmCurveClass.h>
#include <SmLine.h>
#include <SmMemory.h>
#include <SmPrimitiveCreation.h>

namespace
{
enum : int32_t
{
    P5F_OP_EXPORT = 0,
    P5F_OP_MERGE = 1,
};

struct P5fEdgePeek : SmEdge
{
    static SmCurve* SmEdge::*curve() { return &P5fEdgePeek::m_pCurve; }
};

// The Rust import rebuilds NURBs as plain B-splines (and lines, planes) and the analytic
// surfaces the export sends exactly (sc_analytic_export.h) as their own kind: give every
// other face and edge of the C++ B-rep the same geometry (topology_face_p4.cpp's
// fp4_normalise). The analytic originals are leaked (they may still be referenced).
void p5f_normalise(Graph& g)
{
    for (SmFace* f : g.faces)
    {
        SmSurface* s = f->*FPeek::srf();
        if (s == nullptr || !s->IsKindOf(SmBSplineSurface_TYPE) || s->GetType() == SmBSplineSurface_TYPE ||
            s->GetType() == SmPlane_TYPE)
            continue;
        std::vector<double> head;
        if (sc_analytic_export::put_exact_head(head, s, the_context()))
            continue;
        SmBSplineSurface* b = (SmBSplineSurface*)s;
        gw_SURFACE* nurb = b->GetGwNurbPointer();
        if (nurb == nullptr)
            continue;
        SmBSplineSurface* nb = new (the_context()) SmBSplineSurface(sm_AllocateAndCopyNurbSurface(nurb), FALSE, NULL);
        nb->GetOutOfBoundsEnabled() = b->GetOutOfBoundsEnabled();
        nb->SetBSplineSurfaceForm(b->GetBSplineSurfaceForm());
        bool owned = s->GetOwner() == f;
        f->*FPeek::srf() = nb;
        if (owned)
        {
            s->SetOwner(NULL);
            nb->SetOwner(f);
        }
    }
    for (SmEdge* e : g.edges)
    {
        SmCurve* c = e->*P5fEdgePeek::curve();
        if (c == nullptr || !c->IsKindOf(SmBSplineCurve_TYPE) || c->GetType() == SmBSplineCurve_TYPE ||
            c->GetType() == SmLine_TYPE)
            continue;
        SmBSplineCurve* b = (SmBSplineCurve*)c;
        gw_CURVE* nurb = b->GetGwNurbPointer();
        if (nurb == nullptr)
            continue;
        SmBSplineCurve* nc = new (the_context()) SmBSplineCurve(sm_AllocateAndCopyNurbCurve(nurb), b->GetDim(), FALSE);
        bool owned = c->GetOwner() == e;
        e->*P5fEdgePeek::curve() = nc;
        if (owned)
        {
            c->SetOwner(NULL);
            nc->SetOwner(e);
        }
    }
}

// Builds the source B-rep (nullptr with the status on failure).
SmBrep* p5f_source(Rd& rd, const char* path, SmStatus& status)
{
    int64_t kind = rd.idx();
    if (kind == 0)
        return build(rd, path, status);
    // the cone of my_test_separating_closed_edges' second case
    double height = rd.num(), base = rd.num(), top = rd.num(), a0 = rd.num(), a1 = rd.num();
    crv_on_surf_records() = false;
    SmBrep* brep = new (the_context()) SmBrep();
    brep->m_bEditingEnabled = TRUE;
    SmPrimitiveCreation create(brep->GetInfiniteRegion());
    status = create.CreateCone(height, base, top, a0, a1, SmAxis2Placement());
    if (status != SM_SUCCESS)
    {
        delete brep;
        return nullptr;
    }
    return brep;
}

// A point classification's object: 3 face, 7 edge, 9 vertex (with the export index of the
// B-rep before the merge), 100 the classification curve, 101 a surface, -1 none, -2 other.
void p5f_put_object(std::vector<double>& o, Graph& g, const SmObject* p, const SmCurve* curve)
{
    if (p == nullptr)
    {
        o.push_back(-1.0);
        o.push_back(-1.0);
        return;
    }
    if (p == curve)
    {
        o.push_back(100.0);
        o.push_back(-1.0);
        return;
    }
    if (p->IsKindOf(SmSurface_TYPE))
    {
        o.push_back(101.0);
        o.push_back(-1.0);
        return;
    }
    for (int k : {3, 7, 9})
    {
        int64_t i = g.find(k, p);
        if (i >= 0)
        {
            o.push_back((double)k);
            o.push_back((double)i);
            return;
        }
    }
    o.push_back(-2.0);
    o.push_back(-1.0);
}

void p5f_put_pc(std::vector<double>& o, Graph& g, const SmPointClassification& pc, const SmCurve* curve)
{
    o.push_back((double)pc.GetPointClass());
    p5f_put_object(o, g, pc.GetObject(), curve);
}

void p5f_merge(Rd& rd, SmBrep* brep, Graph& g, std::vector<double>& o)
{
    // `face_by face curve_kind curve... do_point_classify has_big_gaps check_3d uv_space_okay`
    const SmContext& ctx = the_context();
    int64_t face_by = rd.idx();
    int64_t fsel = rd.idx();
    SmFace* face = nullptr;
    if (face_by == 0)
        face = g.faces[(size_t)fsel % g.faces.size()];
    else
    {
        // the test's pick: the face whose surface is periodic in U over its natural domain
        // (the harness requires exactly one: -4 none, -5 several)
        int64_t matches = 0;
        for (SmFace* f : g.faces)
            if (f->GetSurface() != nullptr && f->GetSurface()->IsPeriodic(f->GetSurface()->GetNaturalUVDomain(), SM_SP_U))
            {
                if (face == nullptr)
                    face = f;
                matches++;
            }
        if (matches != 1)
        {
            o.push_back(matches == 0 ? -4.0 : -5.0);
            return;
        }
    }
    int64_t curve_kind = rd.idx();
    SmCurve* curve = nullptr;
    if (curve_kind == 0)
    {
        SmPoint3d c = rd.vec();
        SmVector3d x = rd.vec();
        SmVector3d y = rd.vec();
        double a0 = rd.num(), a1 = rd.num(), r = rd.num();
        curve = new (ctx) SmCircle(c, x, y, SmExtent1d(a0, a1), r, 3, &ctx);
    }
    else
    {
        // the UV line between two fractions of the face surface's natural domain
        double fu0 = rd.num(), fv0 = rd.num(), fu1 = rd.num(), fv1 = rd.num();
        SmExtent2d d = face->GetSurface()->GetNaturalUVDomain();
        SmPoint2d p0(d.GetMin().x + fu0 * (d.GetMax().x - d.GetMin().x),
                     d.GetMin().y + fv0 * (d.GetMax().y - d.GetMin().y));
        SmPoint2d p1(d.GetMin().x + fu1 * (d.GetMax().x - d.GetMin().x),
                     d.GetMin().y + fv1 * (d.GetMax().y - d.GetMin().y));
        SmLine uv(p0, p1, 2, &ctx);
        curve = new (ctx) SmCrvOnSurf(uv, *face->GetSurface(), NULL, 2, &ctx);
    }
    SmObjDelete clean_curve(curve);
    SmBoolean do_pc = rd.bflag();
    SmBoolean big = rd.bflag();
    SmBoolean check_3d = rd.bflag();
    SmBoolean uv_ok = rd.bflag();

    SmCurveClassification cls(curve, curve->GetNaturalInterval(), NULL, brep->GetTolerance(), 0, NULL, FALSE, TRUE,
                              face->GetSurface());
    SmStatus s = face->CurveOnClassify(do_pc, cls, big, check_3d, NULL, uv_ok);
    o.push_back(st(s));
    o.push_back((double)cls.GetSize());
    for (ULONG i = 0; i < cls.GetSize(); i++)
    {
        o.push_back(cls[i].m_vInterval.GetMin());
        o.push_back(cls[i].m_vInterval.GetMax());
        p5f_put_pc(o, g, cls[i].m_vStart, curve);
        p5f_put_pc(o, g, cls[i].m_vMid, curve);
        p5f_put_pc(o, g, cls[i].m_vEnd, curve);
    }
    // the test SERs CurveOnClassify before merging
    if (s != SM_SUCCESS)
        return;
    SmTArray<SmFace*> faces;
    SmTArray<SmEdge*> edges;
    SmTArray<SmVertex*> vertices;
    s = cls.MergeIntervals(&faces, &edges, &vertices, face);
    o.push_back(st(s));
    // the test SERs MergeIntervals too
    if (s != SM_SUCCESS)
        return;
    o.push_back((double)faces.GetSize());
    o.push_back((double)edges.GetSize());
    o.push_back((double)vertices.GetSize());
    Graph after;
    walk(brep, after);
    o.push_back((double)after.faces.size());
    o.push_back((double)after.edges.size());
    o.push_back((double)after.vertices.size());
    SmAssertArray pointers, valid;
    bool pointers_ok = brep->ValidatePointers(&pointers, SM_LEVEL_2, NULL) == SM_SUCCESS && pointers.GetSize() == 0;
    bool valid_ok = brep->AssertValid(&valid, SM_LEVEL_2, SM_WALK) && valid.GetSize() == 0;
    o.push_back(pointers_ok && valid_ok ? 1.0 : 0.0);
}
} // namespace

SC_EXPORT int32_t sc_intersect_curve_class_p5_face(int32_t op, const char* path, const double* in, int64_t n_in,
                                                   double* out, int64_t cap, int64_t* out_len)
{
    if (in == NULL || n_in < 1 || out_len == NULL || path == NULL || op < P5F_OP_EXPORT || op > P5F_OP_MERGE)
        return SC_ERR_ARG;
    SC_TRY
    Rd rd{in, 0};
    std::vector<double> o;
    SmStatus status = SM_SUCCESS;
    Built built;
    built.brep = p5f_source(rd, path, status);
    if (built.brep == nullptr)
    {
        o.push_back(-1.0);
        o.push_back(st(status));
    }
    else
    {
        built.brep->m_bEditingEnabled = TRUE;
        Graph g;
        walk(built.brep, g);
        p5f_normalise(g);
        o.push_back(0.0);
        if (g.faces.empty())
            o.push_back(-3.0);
        else if (op == P5F_OP_EXPORT)
            export_graph(built.brep, g, o);
        else
            p5f_merge(rd, built.brep, g, o);
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
