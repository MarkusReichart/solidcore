// SPDX-License-Identifier: Apache-2.0
// The analytic B-rep fixtures of prog_test smTestBrep.cpp (v1.0.1), built with the same
// native constructors: B-rep blocks 10 (sm_MakeLoopOrientationCase), 11
// (sm_MakeTorusVertexLoopCase) and 12 (sm_make_xy_sheet) of sc_topology_fixture.h.
// Each builder returns the B-rep it made (nullptr with the status on failure); the
// caller owns it. The upstream test checks inside the builders become statuses.
#pragma once

#include "StdAfx.h"
#include <NL_Geometry.h>
#include <SmAxis2Placement.h>
#include <SmBSplineCurve.h>
#include <SmBrep.h>
#include <SmBrepData.h>
#include <SmContext.h>
#include <SmCrvOnSurf.h>
#include <SmCurveClass.h>
#include <SmEdgeuse.h>
#include <SmFace.h>
#include <SmLine.h>
#include <SmLoop.h>
#include <SmLoopuse.h>
#include <SmPrimitiveCreation.h>
#include <SmTArray.h>
#include <SmTorus.h>
#include <SmVertex.h>

namespace sc_test_brep
{
// Largest case numbers of blocks 10 and 11 (the upstream loops: 13 orientation cases,
// vertex-loop fixtures `lCase % 10`).
const int64_t MAX_LOOP_ORIENTATION_CASE = 12;
const int64_t MAX_TORUS_VERTEX_LOOP_CASE = 9;

// SER / NER of the upstream fixture: on failure delete the B-rep under construction and
// return its status.
#define SC_TB_SER(expr)                                                                            \
    do                                                                                             \
    {                                                                                              \
        SmStatus sc_tb_s = (expr);                                                                 \
        if (sc_tb_s != SM_SUCCESS)                                                                 \
        {                                                                                          \
            delete pBrep;                                                                          \
            status = sc_tb_s;                                                                      \
            return nullptr;                                                                        \
        }                                                                                          \
    } while (0)
#define SC_TB_NER(expr) SC_TB_SER((expr) ? SM_SUCCESS : SM_ERR)

// sm_make_xy_sheet: a dSize square sheet in the XY plane's positive quadrant.
inline SmBrep* xy_sheet(const SmContext& crContext, double dSize, SmStatus& status)
{
    SmAxis2Placement sRefFrame;
    sRefFrame.SetCanonical(SmPoint3d(0.0, 0.0, 0.0), SmVector3d(1.0, 0.0, 0.0), SmVector3d(0.0, 1.0, 0.0));
    SmBrep* pBrep = SmPrimitiveCreation::CreateRectangle(crContext, 1.0e-5, dSize, dSize, sRefFrame);
    status = pBrep ? SM_SUCCESS : SM_ERR;
    return pBrep;
}

// sm_MakeLoopOrientationCase(lCase): a torus face bounded by one closed SmCrvOnSurf edge
// (a biased-knot polygon, or a rational quadratic UV circle for lCase >= 6), reversed for
// cases 2, 4, 7, 8 and >= 10; cases 8 and 9 add a 4-edge hole, cases 3 and 4 give the edge
// use a stale reversed UV trim.
inline SmBrep* loop_orientation_case(const SmContext& crContext, ULONG lCase, SmStatus& status)
{
    SmBrepConstructor sBuilder;
    SmBrep* pBrep = sBuilder.StartBrep(crContext);
    status = SM_SUCCESS;
    if (pBrep == NULL)
    {
        status = SM_ERR;
        return nullptr;
    }
    pBrep->SetTolerance(1e-7, FALSE);
    SmShell *pOuter = NULL, *pInner = NULL;
    SC_TB_SER(sBuilder.StartShell(pOuter, pInner));
    sBuilder.StartFace(SM_OT_SAME);
    SmTorus* pTorus = NULL;
    SC_TB_SER(SmTorus::CreateCanonical(crContext, SmAxis2Placement(), 2.0, 0.5, pTorus));
    SC_TB_SER(sBuilder.SetFaceSurface(pTorus, pTorus->GetNaturalUVDomain()));
    SmLoop* pLoop = sBuilder.StartLoop(SM_OT_SAME);
    SmPoint2d sPoints[9] = {SmPoint2d(.4, .4), SmPoint2d(.4, .6), SmPoint2d(.6, .6),
                            SmPoint2d(.6, .4), SmPoint2d(.9, .4), SmPoint2d(.9, .9),
                            SmPoint2d(.1, .9), SmPoint2d(.1, .4), SmPoint2d(.4, .4)};
    double dTimes[9] = {0, .3, .6, .9, .911, .922, .933, .944, 1};
    const SmBoolean bReverse = lCase == 2 || lCase == 4 || lCase == 7 || lCase == 8 || lCase >= 10;
    if (lCase == 0)
        for (ULONG lPoint = 0; lPoint < 9; ++lPoint)
            dTimes[lPoint] = double(lPoint) / 8.0;
    if (lCase == 5)
        for (ULONG lPoint = 0; lPoint < 9; ++lPoint)
            sPoints[lPoint] = SmPoint2d(.5, .5) + (sPoints[lPoint] - SmPoint2d(.5, .5)) * 1e-5;
    gw_CPOINT sCP[9];
    double dKnots[12] = {0};
    for (ULONG lPoint = 0; lPoint < 9; ++lPoint)
    {
        ULONG lNext = bReverse ? 8 - lPoint : lPoint;
        sCP[lPoint] = {sPoints[lNext].x, sPoints[lNext].y, NL_NOZ, NL_NOW};
        dKnots[lPoint + 1] = bReverse ? 1.0 - dTimes[lNext] : dTimes[lNext];
    }
    dKnots[10] = 1.0;
    if (lCase >= 6)
    {
        // A rational quadratic UV circle exercises curved-span refinement.
        const SmPoint2d sCircle[9] = {SmPoint2d(.6, .5), SmPoint2d(.6, .6), SmPoint2d(.5, .6),
                                      SmPoint2d(.4, .6), SmPoint2d(.4, .5), SmPoint2d(.4, .4),
                                      SmPoint2d(.5, .4), SmPoint2d(.6, .4), SmPoint2d(.6, .5)};
        const double dCircleKnots[12] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
        for (ULONG lPoint = 0; lPoint < 9; ++lPoint)
        {
            ULONG lNext = bReverse ? 8 - lPoint : lPoint;
            double dWeight = (lPoint % 2) ? smos_Sqrt(.5) : 1.0;
            sCP[lPoint] = {dWeight * sCircle[lNext].x, dWeight * sCircle[lNext].y, NL_NOZ, dWeight};
        }
        for (ULONG lPoint = 0; lPoint < 12; ++lPoint)
            dKnots[lPoint] = dCircleKnots[lPoint];
    }
    SmBSplineCurve sUV(2, 8, lCase >= 6 ? 11 : 10, lCase >= 6 ? 2 : 1, sCP, dKnots, &crContext);
    SmCrvOnSurf* pCurve = new (crContext) SmCrvOnSurf(sUV, *pTorus, NULL, 2, &crContext);
    SmPoint3d sPoint;
    SmStatus sEval = pCurve->EvaluatePoint(0.0, sPoint);
    if (sEval != SM_SUCCESS)
        delete pCurve;
    SC_TB_SER(sEval);
    SmVertex* pVertex = sBuilder.StartVertexOfLoop(sPoint);
    if (!sBuilder.StartEdge(SM_OT_SAME, NULL, pVertex, pVertex))
    {
        delete pCurve;
        SC_TB_SER(SM_ERR);
    }
    SC_TB_SER(sBuilder.SetEdgeCurve(pCurve, FALSE));
    SC_TB_SER(sBuilder.SetLimits(SmVertexNode(sPoint), SmVertexNode(sPoint)));
    SC_TB_SER(sBuilder.EndEdge());
    SC_TB_SER(sBuilder.EndLoop());
    if (lCase == 8 || lCase == 9)
    {
        // Two disjoint CW holes require seams. An outer CCW loop with a
        // contained CW hole does not; checking only one loop is insufficient.
        double dCenter = lCase == 8 ? .2 : .5;
        SmPoint2d sHole[4] = {SmPoint2d(dCenter - .02, dCenter - .02), SmPoint2d(dCenter - .02, dCenter + .02),
                              SmPoint2d(dCenter + .02, dCenter + .02), SmPoint2d(dCenter + .02, dCenter - .02)};
        SC_TB_NER(sBuilder.StartLoop(SM_OT_SAME));
        SmVertex* pVertices[4];
        for (ULONG lPoint = 0; lPoint < 4; ++lPoint)
        {
            SmPoint3d sHolePoint;
            SC_TB_SER(pTorus->EvaluatePoint(sHole[lPoint], sHolePoint));
            pVertices[lPoint] = sBuilder.StartVertexOfLoop(sHolePoint);
        }
        for (ULONG lPoint = 0; lPoint < 4; ++lPoint)
        {
            ULONG lNext = (lPoint + 1) % 4;
            SC_TB_NER(sBuilder.StartEdge(SM_OT_SAME, NULL, pVertices[lPoint], pVertices[lNext]));
            SmLine sLine(sHole[lPoint], sHole[lNext], 2, &crContext);
            SC_TB_SER(sBuilder.SetEdgeCurve(new (crContext) SmCrvOnSurf(sLine, *pTorus, NULL, 2, &crContext), FALSE));
            SC_TB_SER(sBuilder.SetLimits(SmVertexNode(pVertices[lPoint]->GetPoint()),
                                         SmVertexNode(pVertices[lNext]->GetPoint())));
            SC_TB_SER(sBuilder.EndEdge());
        }
        SC_TB_SER(sBuilder.EndLoop());
    }
    SC_TB_SER(sBuilder.EndFace());
    SC_TB_SER(sBuilder.EndShell(FALSE));
    SC_TB_SER(sBuilder.EndBrep());

    SmTArray<SmEdgeuse*> sUses;
    pLoop->GetLoopuse()->GetEdgeuses(sUses);
    if (lCase == 3 || lCase == 4)
    {
        // A plausible closed cached trim with the opposite orientation.
        // TRUE must still measure the defining 3D curve, leaving this cache intact.
        gw_CPOINT sReverseCP[9];
        double dReverseKnots[11];
        for (ULONG lPoint = 0; lPoint < 9; ++lPoint)
            sReverseCP[lPoint] = sCP[8 - lPoint];
        for (ULONG lPoint = 0; lPoint < 11; ++lPoint)
            dReverseKnots[lPoint] = 1.0 - dKnots[10 - lPoint];
        sUses[0]->SetUVTrimCurve(new (crContext) SmBSplineCurve(2, 8, 10, 1, sReverseCP, dReverseKnots, &crContext),
                                 0.0, TRUE);
    }
    return pBrep;
}

// sm_MakeTorusVertexLoopCase(lCase): a torus face with a 4-edge SmCrvOnSurf UV-square loop,
// then 1, 2, 3 or 5 vertex loops at (.2, .3) and on the seams (fixture `lCase % 10`), made
// with editing enabled while a Boolean is in progress, then HasSeamProblem with the
// fixture's two checks (both seams missing, none crossed).
inline SmBrep* torus_vertex_loop_case(SmContext& rContext, ULONG lCase, SmStatus& status)
{
    SmBrepConstructor sBuilder;
    SmBrep* pBrep = sBuilder.StartBrep(rContext);
    status = SM_SUCCESS;
    if (pBrep == NULL)
    {
        status = SM_ERR;
        return nullptr;
    }
    SmShell *pOuter = NULL, *pInner = NULL;
    pBrep->SetTolerance(1e-7, FALSE);
    SC_TB_SER(sBuilder.StartShell(pOuter, pInner));
    SmFace* pFace = sBuilder.StartFace(SM_OT_SAME);
    SmTorus* pTorus = NULL;
    SC_TB_SER(SmTorus::CreateCanonical(rContext, SmAxis2Placement(), 2.0, 0.5, pTorus));
    SC_TB_SER(sBuilder.SetFaceSurface(pTorus, pTorus->GetNaturalUVDomain()));
    SC_TB_NER(sBuilder.StartLoop(SM_OT_SAME));
    SmPoint2d sUV[4] = {SmPoint2d(.4, .4), SmPoint2d(.4, .6), SmPoint2d(.6, .6), SmPoint2d(.6, .4)};
    SmPoint3d sPoints[4];
    SmVertex* pHoleVertices[4];
    for (ULONG lPoint = 0; lPoint < 4; ++lPoint)
    {
        SC_TB_SER(pTorus->EvaluatePoint(sUV[lPoint], sPoints[lPoint]));
        pHoleVertices[lPoint] = sBuilder.StartVertexOfLoop(sPoints[lPoint]);
    }
    for (ULONG lPoint = 0; lPoint < 4; ++lPoint)
    {
        ULONG lNext = (lPoint + 1) % 4;
        SmLine sUVLine(sUV[lPoint], sUV[lNext], 2, &rContext);
        SmCrvOnSurf* pCurve = new (rContext) SmCrvOnSurf(sUVLine, *pTorus, NULL, 2, &rContext);
        if (!sBuilder.StartEdge(SM_OT_SAME, NULL, pHoleVertices[lPoint], pHoleVertices[lNext]))
        {
            delete pCurve;
            SC_TB_SER(SM_ERR);
        }
        SC_TB_SER(sBuilder.SetEdgeCurve(pCurve, FALSE));
        SC_TB_SER(sBuilder.SetLimits(SmVertexNode(sPoints[lPoint]), SmVertexNode(sPoints[lNext])));
        SC_TB_SER(sBuilder.EndEdge());
    }
    SC_TB_SER(sBuilder.EndLoop());
    SC_TB_SER(sBuilder.EndFace());
    SC_TB_SER(sBuilder.EndShell(FALSE));
    SC_TB_SER(sBuilder.EndBrep());
    pBrep->m_bEditingEnabled = TRUE;
    SmTemporaryChangeValue<SmBoolean> sDoingBool(rContext.GetDoingBooleanRef(), TRUE);
    ULONG lFixture = lCase % 10;
    ULONG lKind = lFixture / 2;
    SmBoolean bFirstU = lFixture % 2 == 0;
    ULONG lNumPoints = lKind == 0 ? 1 : lKind == 3 ? 2 : lKind == 4 ? 5 : 3;
    for (ULONG lPoint = 0; lPoint < lNumPoints; ++lPoint)
    {
        SmPoint2d sPointUV(.2, .3);
        if (lPoint > 0)
        {
            SmBoolean bOnU = lKind == 4 ? (lPoint <= 2) == bFirstU : lKind == 1 ? bFirstU : !bFirstU;
            double dParameter = .25 + .5 * ((lPoint - 1) % 2);
            sPointUV = lKind == 3 ? SmPoint2d(0, 0) : bOnU ? SmPoint2d(0, dParameter) : SmPoint2d(dParameter, 0);
        }
        SmPoint3d sPoint;
        SmLoop* pVertexLoop = NULL;
        SmVertex* pVertex = NULL;
        SC_TB_SER(pTorus->EvaluatePoint(sPointUV, sPoint));
        SC_TB_SER(pBrep->MakeVertexLoop(pFace, sPoint, pVertexLoop, pVertex));
    }
    SmCurveClassification sU, sV;
    SmSurfParamType eCrossed = SM_SP_NEITHER, eMissing = SM_SP_NEITHER, eNear = SM_SP_NEITHER;
    pFace->HasSeamProblem(sU, sV, eCrossed, eMissing, eNear);
    // "Vertex loops hid a missing torus seam", "Vertex loops were counted as seam crossings"
    SC_TB_SER(eMissing == SM_SP_BOTH ? SM_SUCCESS : SM_ERR);
    SC_TB_SER(eCrossed == SM_SP_NEITHER ? SM_SUCCESS : SM_ERR);
    return pBrep;
}

#undef SC_TB_NER
#undef SC_TB_SER
} // namespace sc_test_brep
