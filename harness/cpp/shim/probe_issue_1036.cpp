// SPDX-License-Identifier: Apache-2.0
// One-off probe for tangency-labs/solidcore#1036 (not a unit shim): builds two capped
// skin solids with SmApiCreateSkinPrimitive from exact NURB section curves and runs one
// two-B-rep Boolean on them (SmApiBooleanDifference / SmApiBooleanIntersection /
// SmApiBooleanUnion), or measures the first solid alone.
//
// One entry point, sc_probe_issue_1036_run(in, n_in, out, cap, out_len). `in` is
// `op skinA skinB`, a skin `nc curve*nc degree cap`, a curve `dim n p m (n+1) x [x y z w]
// (m+1) x knot` (as sweep_api_surfaces.cpp's CK_NURB block, without the kind). Layout and
// validation: crates/solidcore-ffi/src/units/probe_issue_1036.rs. Every read is bounds
// checked here as well (a malformed program is SC_ERR_ARG).
//
// `out`: per skin `status [faces edges vertices volume_status volume]` (bracketed on
// success), then for ops 1..3 the Boolean `status [has faces edges vertices volume_status
// volume]`.
#include "sc_shim.h"

#include "StdAfx.h"
#include <SmApiBrep.h>
#include <SmApiGeneral.h>
#include <SmApiPrimitives.h>
#include <SmApiQueries.h>
#include <SmBSplineCurve.h>
#include <SmBrep.h>
#include <SmContext.h>
#include <SmCurve.h>
#include <SmTArray.h>
#include <nurbs.h>

#include <stdexcept>
#include <vector>

#include "sc_sm.h"

using sc::flg;
using sc::st;

namespace
{
enum : int32_t
{
    OP_MEASURE = 0,
    OP_DIFFERENCE = 1,
    OP_INTERSECTION = 2,
    OP_UNION = 3,
};

struct Bad : std::runtime_error
{
    Bad() : std::runtime_error("bad program") {}
};

struct Rd
{
    const double* in;
    int64_t n;
    int64_t at;
    double num()
    {
        if (at >= n)
            throw Bad();
        return in[at++];
    }
    int64_t idx(int64_t lo, int64_t hi)
    {
        double d = num();
        if (!(d >= (double)lo && d <= (double)hi) || d != (double)(int64_t)d)
            throw Bad();
        return (int64_t)d;
    }
};

struct NurbSpec
{
    int64_t dim = 3, n = 0, p = 0, m = 0;
    std::vector<double> data;
};

struct SkinSpec
{
    std::vector<NurbSpec> curves;
    ULONG degree = 1;
    bool cap = true;
};

NurbSpec read_nurb(Rd& rd)
{
    NurbSpec c;
    c.dim = rd.idx(1, 3);
    c.n = rd.idx(1, 4096);
    c.p = rd.idx(1, 17);
    c.m = rd.idx(2, 4096 + 18);
    if (c.m != c.n + c.p + 1)
        throw Bad();
    for (int64_t i = 0; i < 4 * (c.n + 1) + c.m + 1; i++)
        c.data.push_back(rd.num());
    return c;
}

SkinSpec read_skin(Rd& rd)
{
    SkinSpec s;
    int64_t nc = rd.idx(2, 64);
    for (int64_t i = 0; i < nc; i++)
        s.curves.push_back(read_nurb(rd));
    s.degree = (ULONG)rd.idx(1, 17);
    s.cap = rd.idx(0, 1) != 0;
    return s;
}

SmBSplineCurve* build_nurb(const NurbSpec& c)
{
    const SmContext& ctx = *SmApiGetOrCreateContext();
    gw_CURVE* cur = sm_AllocateNurbCurve((gw_INDEX)c.n, (gw_DEGREE)c.p, (gw_INDEX)c.m);
    const double* d = c.data.data();
    int64_t at = 0;
    for (int64_t i = 0; i <= c.n; i++)
    {
        cur->pol->Pw[i].x = d[at++];
        cur->pol->Pw[i].y = d[at++];
        cur->pol->Pw[i].z = d[at++];
        cur->pol->Pw[i].w = d[at++];
    }
    for (int64_t i = 0; i <= c.m; i++)
        cur->knt->U[i] = d[at++];
    return new (ctx) SmBSplineCurve(cur, (ULONG)c.dim, FALSE);
}

// Builds the skin (NULL on failure) and appends `status [faces edges vertices
// volume_status volume]`. One profile curve per section, as the #1036 test passes.
SmBrep* build_skin(std::vector<double>& o, const SkinSpec& s)
{
    std::vector<SmBSplineCurve*> owned;
    SmTArray<ULONG> counts;
    SmTArray<const SmCurve*> curves;
    for (size_t i = 0; i < s.curves.size(); i++)
    {
        owned.push_back(build_nurb(s.curves[i]));
        counts.Add(1);
        curves.Add(owned.back());
    }
    SmBrep* res = NULL;
    SmApiStatus stat =
        SmApiCreateSkinPrimitive(counts, curves, s.degree, s.cap ? TRUE : FALSE, res, NULL, NULL, NULL);
    for (size_t i = 0; i < owned.size(); i++)
        delete owned[i];
    o.push_back(st(stat));
    if (stat != SM_SUCCESS || res == NULL)
    {
        if (stat == SM_SUCCESS)
            o.push_back(-1.0); // success without a B-rep
        delete res;
        return NULL;
    }
    o.push_back((double)res->GetNumFaces());
    o.push_back((double)res->GetNumEdges());
    o.push_back((double)res->GetNumVertices());
    double volume = 0.0;
    SmApiStatus vstat = SmApiBrepComputeVolume(res, 1e-6, volume);
    o.push_back(st(vstat));
    o.push_back(vstat == SM_SUCCESS ? volume : 0.0);
    return res;
}

void run(Rd& rd, std::vector<double>& o)
{
    int32_t op = (int32_t)rd.idx(0, 3);
    SkinSpec a = read_skin(rd);
    SkinSpec b = read_skin(rd);
    if (rd.at != rd.n)
        throw Bad();

    SmBrep* ba = build_skin(o, a);
    SmBrep* bb = build_skin(o, b);
    if (op == OP_MEASURE || ba == NULL || bb == NULL)
    {
        delete ba;
        delete bb;
        return;
    }
    SmBrep* res = NULL;
    SmApiStatus stat = SM_ERR;
    switch (op)
    {
    case OP_DIFFERENCE: stat = SmApiBooleanDifference(ba, bb, res); break;
    case OP_INTERSECTION: stat = SmApiBooleanIntersection(ba, bb, res); break;
    case OP_UNION: stat = SmApiBooleanUnion(ba, bb, res); break;
    default: throw Bad();
    }
    o.push_back(st(stat));
    if (stat != SM_SUCCESS)
    {
        // Both operands remain allocated on failure.
        if (res != NULL && res != ba && res != bb)
            delete res;
        delete ba;
        delete bb;
        return;
    }
    // On success pBrep1 becomes the result and pBrep2 is deleted (SmApiBrep.h).
    o.push_back(flg(res != NULL));
    if (res == NULL)
        return;
    o.push_back((double)res->GetNumFaces());
    o.push_back((double)res->GetNumEdges());
    o.push_back((double)res->GetNumVertices());
    double volume = 0.0;
    SmApiStatus vstat = SmApiBrepComputeVolume(res, 1e-6, volume);
    o.push_back(st(vstat));
    o.push_back(vstat == SM_SUCCESS ? volume : 0.0);
    delete res;
}
} // namespace

SC_EXPORT int32_t sc_probe_issue_1036_run(const double* in, int64_t n_in, double* out,
                                         int64_t cap, int64_t* out_len)
{
    if (in == NULL || n_in < 0 || out_len == NULL)
        return SC_ERR_ARG;
    SC_TRY
    SmApiGetOrCreateContext();
    std::vector<double> o;
    try
    {
        Rd rd{in, n_in, 0};
        run(rd, o);
    }
    catch (const Bad&)
    {
        return SC_ERR_ARG;
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
