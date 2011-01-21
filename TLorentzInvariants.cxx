#include "TLorentzInvariants.h"
#include "TRadCor.h"
#include "HapradErrors.h"
#include "TGlobalConfig.h"
#include "TKinematicalVariables.h"
#include "THadronKinematics.h"
#include "THapradException.h"
#include "haprad_constants.h"
#include "square_power.h"
#include <iostream>


TLorentzInvariants::TLorentzInvariants(const TRadCor* rc)
 : fS(0), fX(0), fSx(0), fSp(0), fQ2(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fHadKin = rc->GetHadronKinematics();
}



TLorentzInvariants::~TLorentzInvariants()
{
    // Do nothing
}



void TLorentzInvariants::Evaluate(void)
{
    using namespace TMath;

    Double_t M = kMassProton;
    Double_t m;

    switch(fConfig->PolarizationType()) {
        case 1:
            m = kMassElectron;
            break;
        case 2:
            m = kMassMuon;
            break;
        default:
            m = kMassElectron;
    }

    fS = 2. * M * fKin->E();

    Double_t y;
    if (fKin->Y() >= 0.) {
        y = - fKin->Y();
        fQ2 = fS * fKin->X() * fKin->Y();
    } else {
        fQ2 = - fKin->Y();
        y = fQ2 / (fS * fKin->X());
    }


    Double_t y_max = 1. / (1. + SQ(M) * fKin->X() / fS);
    Double_t y_min = (kMassC2 - SQ(M)) / (fS * (1. - fKin->X()));

    if (y > y_max || y < y_min || fKin->X() > 1. || fKin->X() < 0.)
        throw TKinematicException();

    fX  = fS * (1. - y);
    fSx = fS - fX;
    fSp = fS + fX;
    fW2 = fS - fX - fQ2 + SQ(M);

    fLambdaS = SQ(fS) - 4. * SQ(m) * SQ(M);
    fLambdaX = SQ(fX) - 4. * SQ(m) * SQ(M);
    fLambdaM = SQ(fQ2) + 4. * SQ(m) * fQ2;
    fLambdaQ = SQ(fSx) + 4. * SQ(M) * fQ2;


    if (fLambdaS < 0) HAPRAD_WARN_MSG("Invariants", "lambda_s < 0");
    if (fLambdaX < 0) HAPRAD_WARN_MSG("Invariants", "lambda_x < 0");
    if (fLambdaQ < 0) HAPRAD_WARN_MSG("Invariants", "lambda_q < 0");
    if (fLambdaM < 0) HAPRAD_WARN_MSG("Invariants", "lambda_m < 0");

    fSqrtLs = Sqrt(Max(0.,fLambdaS));
    fSqrtLx = Sqrt(Max(0.,fLambdaX));
    fSqrtLq = Sqrt(Max(0.,fLambdaQ));
}



void TLorentzInvariants::EvaluateV12(void)
{
    using namespace TMath;

    Double_t M = kMassProton;
    Double_t m;

    switch(fConfig->PolarizationType()) {
        case 1:
            m = kMassElectron;
            break;
        case 2:
            m = kMassMuon;
            break;
        default:
            m = kMassElectron;
    }

    Double_t costs, costx, sints, sintx;
    Double_t lambda;

    costs = (fS * (fS - fX) + 2. * SQ(M) * fQ2) / fSqrtLs / fSqrtLq;
    costx = (fX * (fS - fX) - 2. * SQ(M) * fQ2) / fSqrtLx / fSqrtLq;

    lambda = fS * fX * fQ2 - SQ(M) * SQ(fQ2) - SQ(m) * fLambdaQ;

    if (lambda > 0) {
        sints = 2. * M * Sqrt(lambda) / fSqrtLs / fSqrtLq;
        sintx = 2. * M * Sqrt(lambda) / fSqrtLx / fSqrtLq;
    } else {
        HAPRAD_WARN_MSG("sphi", "sints = NaN");
        HAPRAD_WARN_MSG("sphi", "sintx = NaN");
        sints = 0.;
        sintx = 0.;
    }

    Double_t v1, v2;
    v1 = costs * fHadKin->Pl() + sints * fHadKin->Pt() * Cos(fKin->PhiH());
    v2 = costx * fHadKin->Pl() + sintx * fHadKin->Pt() * Cos(fKin->PhiH());

    fV1 = (fS * fHadKin->Eh() - fSqrtLs * v1) / M;
    fV2 = (fX * fHadKin->Eh() - fSqrtLx * v2) / M;
}
