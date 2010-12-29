#include "TBorn.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "THadronKinematics.h"
#include "haprad_constants.h"

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif


TBorn::TBorn(const TRadCor* rc)
  : fH(rc)
{
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
    fHadKin = rc->GetHadronKinematics();

    const Double_t& M   = kMassProton;
    const Double_t& m_h = kMassDetectedHadron;

    fThetaB[0] = fInv->Q2();
    fThetaB[1] = (fInv->S() * fInv->X() - M * M * fInv->Q2()) / 2.;
    fThetaB[2] = (fInv->V1() * fInv->V2() - m_h * m_h * fInv->Q2()) / 2.;
    fThetaB[3] = (fInv->V2() * fInv->S() + fInv->V1() * fInv->X() -
                            fKin->Z() * fInv->Q2() * fInv->Sx()) / 2.;
}



Double_t TBorn::GetValue(Double_t N)
{
    fH.Evaluate(0.,0.,0.);

    Double_t sum = 0.;
    for (Int_t i = 0; i < 4; ++i) {
        sum = sum + fThetaB[i] * fH[i];
    }

    return sum * N / fInv->Q2() / fInv->Q2() * 2.;
}
