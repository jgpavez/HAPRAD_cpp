#include "TDelta.h"
#include "TRadCor.h"
#include "TGlobalConfig.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "THapradUtils.h"
#include "haprad_constants.h"
#include "square_power.h"


TDelta::TDelta(const TRadCor* rc)
  : fVR(0), fInf(0), fVac(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();
}



TDelta::~TDelta()
{
    // Do nothing
}



void TDelta::Evaluate(void)
{
    using namespace TMath;

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

    Double_t S_ = fInv->S() - fInv->Q2() - fInv->V1();
    Double_t X_ = fInv->X() + fInv->Q2() - fInv->V2();

    Double_t l_m  = Log(fInv->Q2() / SQ(m));
    Double_t Li_2 = HapradUtils::fspen(1. - fHadKin->Px2() * fInv->Q2() / (S_ * X_));

    fVR  = 1.5 * l_m - 2. - 0.5 * SQ(Log(X_ / S_)) + Li_2 - SQ(kPi) / 6.;
    fInf = (l_m - 1.) * Log(SQ(fHadKin->Px2() - kMassC2) / S_ / X_);
    fVac = VacPol(fInv->Q2());

    fVR  *= kAlpha / kPi;
    fInf *= kAlpha / kPi;
    fVac *= kAlpha / kPi;
}



Double_t TDelta::VacPol(const Double_t Q2)
{
    Double_t leptonMass[3] = { 0.26110E-6,
                               0.111637E-1,
                               3.18301 };

    Double_t suml = 0;
    for (Int_t i = 0; i < 3; ++i) {
        Double_t a2    = 2. * leptonMass[i];
        Double_t sqlmi = TMath::Sqrt(SQ(Q2) + 2. * a2 * Q2);
        Double_t allmi = TMath::Log((sqlmi + Q2) / (sqlmi - Q2)) / sqlmi;

        suml = suml + 2. * (Q2 + a2) * allmi / 3. - 10. / 9. +
                    4. * a2 * (1. - a2 * allmi) / 3. / Q2;
    }

    Double_t a, b, c;

    if (Q2 < 1) {
        a = -1.345E-9;
        b = -2.302E-3;
        c = 4.091;
    } else if (Q2 < 64) {
        a = -1.512E-3;
        b = -2.822E-3;
        c =  1.218;
    } else {
        a = -1.1344E-3;
        b = -3.0680E-3;
        c =  9.9992E-1;
    }

    Double_t sumh;
    sumh = - (a + b * TMath::Log(1. + c * Q2)) * 2. * kPi / kAlpha;

    return suml + sumh;
}
