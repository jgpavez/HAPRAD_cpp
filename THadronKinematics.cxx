#include "THadronKinematics.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "THapradException.h"
#include "TLorentzInvariants.h"
#include "haprad_constants.h"
#include "square_power.h"

#include <iostream>
#ifdef DEBUG
#include <iomanip>
#endif


THadronKinematics::THadronKinematics(const TRadCor* rc)
 : fEh(0), fPl(0), fPt(0), fNu(0), fPx2(0), fPh(0)
{
    fConfig = rc->GetConfig();
    fKin    = rc->GetKinematicalVariables();
    fInv    = rc->GetLorentzInvariants();
}



THadronKinematics::~THadronKinematics()
{
    // Do nothing
}



void THadronKinematics::SetNu(void)
{
    fNu = fInv->Sx() / (2 * kMassProton);
}



void THadronKinematics::SetEh(void)
{
    fEh = fNu * fKin->Z();

    if (fEh < kMassDetectedHadron) {
        std::cout << "    E_h: " << fEh
                  << std::endl;
        throw TKinematicException();
    }
}



void THadronKinematics::SetSqNuQ(void)
{
    fSqNuQ = TMath::Sqrt(SQ(fNu) + fInv->Q2());
}



void THadronKinematics::SetMomentum(void)
{
    fPh = TMath::Sqrt(SQ(fEh) - SQ(kMassDetectedHadron));

    if (fKin->T() >= 0) {
        fPt = fKin->T();

        if (fPh < fPt) {
            std::cout << "    p_h: " << fPh
                      << "    p_t: " << fPt
                      << std::endl;
            throw TKinematicException();
        }

        fPl = TMath::Sqrt(SQ(fPh) - SQ(fPt));
    } else {
        fPl = (fKin->T() + fInv->Q2() - SQ(kMassDetectedHadron) +
                2 * fNu * fEh) / 2 / fSqNuQ;

        if (fPh < TMath::Abs(fPl)) {
            Double_t eps1, eps2, eps3, eps4, eps5, sum;

            eps1 = fKin->T() * kEpsMachine / fSqNuQ;
            eps2 = 2 * SQ(kMassDetectedHadron)  * kEpsMachine / fSqNuQ;
            eps3 = 2 * fNu * fEh * kEpsMachine / fSqNuQ;
            eps4 = fKin->T() + fInv->Q2() - SQ(kMassDetectedHadron) + 2 * fNu * fEh;
            eps5 = eps4 / fSqNuQ * kEpsMachine;

            sum = SQ(eps1) + SQ(eps2) + 2 * SQ(eps3) + SQ(eps5);

            Double_t epspl  = TMath::Sqrt(sum) / 2;
            Double_t mineps = fPh - TMath::Abs(fPl);

            if (TMath::Abs(mineps) > epspl) {
                std::cout << "    p_h: " << fPh
                          << "    p_l: " << fPl
                          << std::endl;
                throw TKinematicException();
            } else {
               std::cout << "    Zero p_t! " << fPl
                         << "\t"             << mineps
                         << "\t"             << epspl << std::endl;
               fPl = TMath::Sign(1., fPl) * fPh;
            }
        }

        fPt = TMath::Sqrt(SQ(fPh) - SQ(fPl));
    }
}



void THadronKinematics::SetPx2(void)
{
    Double_t px2_max = fInv->W2() - SQ(kMassDetectedHadron);

    fPx2 = SQ(kMassProton) + fInv->Sx() * (1 - fKin->Z()) + fKin->T();

    if (fPx2 < kMassC2 || fPx2 > px2_max) {
        std::cout << "    px2:    " << fPx2
                  << "    mc2:    " << kMassC2
                  << "    px2max: " << px2_max
                  << std::endl;
        throw TKinematicException();
    }
}



void THadronKinematics::Clear(void)
{
    fEh = 0;
    fPl = 0;
    fPt = 0;
    fNu = 0;
    fSqNuQ = 0;
    fPx2 = 0;
    fPh = 0;
}
