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
 : fS(0), fX(0), fSx(0), fSp(0), fQ2(0), fY(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fHadKin = rc->GetHadronKinematics();
}



TLorentzInvariants::~TLorentzInvariants()
{
    // Do nothing
}



void TLorentzInvariants::SetS()
{
    fS = 2 * kMassProton * fKin->E();
}



void TLorentzInvariants::SetX()
{
    fX  = fS * (1 - fY);
}



void TLorentzInvariants::SetSx()
{
    fSx = fS - fX;
}



void TLorentzInvariants::SetSp()
{
    fSp = fS + fX;
}



void TLorentzInvariants::SetQ2()
{
    if (fKin->Y() >= 0.) {
        fY = - fKin->Y();
        fQ2 = fS * fKin->X() * fY;
    } else {
        fQ2 = - fKin->Y();
        fY = fQ2 / (fS * fKin->X());
    }
}



void TLorentzInvariants::SetW2()
{
    fW2 = fS - fX - fQ2 + SQ(kMassProton);
}



void TLorentzInvariants::SetLambdas()
{
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

    fLambdaS = SQ(fS)  - 4 * SQ(m) * SQ(M);
    fLambdaX = SQ(fX)  - 4 * SQ(m) * SQ(M);
    fLambdaQ = SQ(fSx) + 4 * SQ(M) * fQ2;
    fLambdaM = SQ(fQ2) + 4 * SQ(m) * fQ2;

    if (fLambdaS < 0) HAPRAD_WARN_MSG("Invariants", "lambda_s < 0");
    if (fLambdaX < 0) HAPRAD_WARN_MSG("Invariants", "lambda_x < 0");
    if (fLambdaQ < 0) HAPRAD_WARN_MSG("Invariants", "lambda_q < 0");
    if (fLambdaM < 0) HAPRAD_WARN_MSG("Invariants", "lambda_m < 0");

    fSqrtLs = TMath::Sqrt(TMath::Max(0.,fLambdaS));
    fSqrtLq = TMath::Sqrt(TMath::Max(0.,fLambdaQ));
    fSqrtLx = TMath::Sqrt(TMath::Max(0.,fLambdaX));
}



void TLorentzInvariants::SetV12(void)
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

    costs = (fS * (fS - fX) + 2 * SQ(M) * fQ2) / fSqrtLs / fSqrtLq;
    costx = (fX * (fS - fX) - 2 * SQ(M) * fQ2) / fSqrtLx / fSqrtLq;

    lambda = fS * fX * fQ2 - SQ(M) * SQ(fQ2) - SQ(m) * fLambdaQ;

    if (lambda > 0) {
        sints = 2 * M * Sqrt(lambda) / fSqrtLs / fSqrtLq;
        sintx = 2 * M * Sqrt(lambda) / fSqrtLx / fSqrtLq;
    } else {
        HAPRAD_WARN_MSG("sphi", "sints = NaN");
        HAPRAD_WARN_MSG("sphi", "sintx = NaN");
        sints = 0;
        sintx = 0;
    }

    Double_t v1, v2;
    v1 = costs * fHadKin->Pl() + sints * fHadKin->Pt() * Cos(fKin->PhiH());
    v2 = costx * fHadKin->Pl() + sintx * fHadKin->Pt() * Cos(fKin->PhiH());

    fV1 = (fS * fHadKin->Eh() - fSqrtLs * v1) / M;
    fV2 = (fX * fHadKin->Eh() - fSqrtLx * v2) / M;
}
