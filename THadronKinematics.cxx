#include "THadronKinematics.h"
#include "TRadCor.h"
#include "THapradConfig.h"
#include "TLorentzInvariants.h"
#include "TKinematicalVariables.h"
#include "THapradException.h"
#include "HapradErrors.h"
#include "haprad_constants.h"
#include "square_power.h"

#include <iostream>

THadronKinematics::THadronKinematics(const THapradConfig* config,
                                     const TKinematicalVariables* kin,
                                     const TLorentzInvariants* inv)
 : fConfig(config), fKin(kin), fInv(inv),
   fEh(0), fPl(0), fPt(0), fNu(0), fPx2(0), fPh(0)
{
    // Do nothing
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



void THadronKinematics::SetV12(void)
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

    // Just use short names to make formulas clear
    const Double_t S  = fInv->S();
    const Double_t X  = fInv->X();
    const Double_t Q2 = fInv->Q2();

    const Double_t SqrtLs = fInv->SqrtLs();
    const Double_t SqrtLx = fInv->SqrtLx();
    const Double_t SqrtLq = fInv->SqrtLq();


    costs = (S * (S - X) + 2 * SQ(M) * Q2) / SqrtLs / SqrtLq;
    costx = (X * (S - X) - 2 * SQ(M) * Q2) / SqrtLx / SqrtLq;

    lambda = S * X * Q2 - SQ(M) * SQ(Q2) - SQ(m) * fInv->LambdaQ();

    if (lambda > 0) {
        sints = 2 * M * Sqrt(lambda) / SqrtLs / SqrtLq;
        sintx = 2 * M * Sqrt(lambda) / SqrtLx / SqrtLq;
    } else {
        HAPRAD_WARN_MSG("sphi", "sints = NaN");
        HAPRAD_WARN_MSG("sphi", "sintx = NaN");
        sints = 0;
        sintx = 0;
    }

    Double_t v1, v2;
    v1 = costs * fPl + sints * fPt * Cos(fKin->PhiH());
    v2 = costx * fPl + sintx * fPt * Cos(fKin->PhiH());

    fV1 = (fInv->S() * fEh - SqrtLs * v1) / M;
    fV2 = (fInv->X() * fEh - SqrtLx * v2) / M;
}
