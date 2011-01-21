#include "TRadCor.h"
#include "THapradException.h"
#include "TDelta.h"
#include "TQQTPhi.h"
#include "TBorn.h"
#include "TRV2TR.h"
#include "haprad_constants.h"
#include "square_power.h"
#include "Math/GSLIntegrator.h"
#include "Math/GSLMCIntegrator.h"
#include <iostream>
#include <iomanip>


TRadCor::TRadCor()
: fInv(this), fHadKin(this),
  maxMx2(0.), Mx2(0.),
  sigma_born(0.), sig_obs(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron)
{
    // Default constructor
}



TRadCor::TRadCor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                 Double_t p_t, Double_t phi, Double_t maxMx2)
: fInv(this), fHadKin(this),
  maxMx2(0.), Mx2(0.),
  sigma_born(0.), sig_obs(0.), tail(0.),
  M(kMassProton), m(kMassElectron), m_h(kMassDetectedHadron)
{
    // Normal constructor for a radiative correction object
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount
    // of missing mass.

    fKin.SetAll(x,-Q2,z,p_t,phi/kRadianDeg,E);
    Setup();
}



TRadCor::~TRadCor()
{
    // Default destructor
}



void TRadCor::SetParameters(Double_t E, Double_t x, Double_t Q2, Double_t z,
                            Double_t p_t, Double_t phi, Double_t maxMx2)
{
    // Set the values for the variables used to calculate the radiative
    // correction.
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount
    // of missing mass.

    this->maxMx2 = maxMx2;

    fKin.SetAll(x,-Q2,z,p_t,phi/kRadianDeg,E);
    Setup();
}



void TRadCor::SetEbeam(Double_t E)
{
    fKin.SetE(E);
}



void TRadCor::SetX(Double_t x)
{
    fKin.SetX(x);
    Setup();
}



void TRadCor::SetQ2(Double_t Q2)
{
    fKin.SetY(-Q2);
    Setup();
}



void TRadCor::SetZ(Double_t z)
{
    fKin.SetZ(z);
    Setup();
}



void TRadCor::SetPt(Double_t p_t)
{
    fKin.SetT(p_t);
    Setup();
}



void TRadCor::SetPhi(Double_t phi)
{
    fKin.SetPhiH(phi/kRadianDeg);
}



void TRadCor::SetMaxMx(Double_t maxMx2)
{
    this->maxMx2 = maxMx2;
}



void TRadCor::RegisteredLepton(Int_t type)
{
    // Set the registere lepton: 1 -- electron; 2 -- muon.
    //
    // Default is 1.

    fConfig.SetLepton(type);
}



void TRadCor::IntegratePhiRad(Int_t type)
{
    // Set whether to integrate over phi_{rad} (1) or to approximate (0).
    //
    // Default is 0;

    fConfig.SetIntegrationPhiRad(type);
}



void TRadCor::IntegratePhiHad(Int_t type)
{
    // Set whether to integrate over phi_{had} (1) or not (0).
    //
    // Default is 0;

    fConfig.SetIntegrationPhiHad(type);
}



void TRadCor::SetPolarization(Int_t type)
{
    // Set the type of the polarization:
    //
    //     1 -- long; 2 -- tran; 0 -- unpol.
    //
    // Default is 0;

    fConfig.SetPolarization(type);
}



void TRadCor::Setup(void)
{
    // Calculate the missing mass squared

    Double_t S_x = - fKin.Y() / fKin.X();
    Mx2 = SQ(M) + S_x * (1.0 - fKin.Z()) + fKin.T();
}



void TRadCor::CalculateRCFactor(void)
{
    if (Mx2 > maxMx2)
        Haprad();
}



void TRadCor::CalculateRCFactor(Double_t E, Double_t x, Double_t Q2,
                                Double_t z, Double_t p_t, Double_t phi,
                                Double_t maxMx2)
{
    SetParameters(E,x,Q2,z,p_t,phi,maxMx2);
    CalculateRCFactor();
}



Double_t TRadCor::GetFactor1(void)
{
    Double_t sigma_obs_f1 = sig_obs + tai[0];

    return sigma_obs_f1 / sigma_born;
}



Double_t TRadCor::GetFactor2(void)
{
    Double_t sigma_obs_f2 = sig_obs + tai[1];

    return sigma_obs_f2 / sigma_born;
}



Double_t TRadCor::GetFactor3(void)
{
    Double_t sigma_obs_f3 = sig_obs + tai[0] + tai[1] / 2;

    return sigma_obs_f3 / sigma_born;
}



Double_t TRadCor::GetRCFactor(void)
{
    // Get the radiative correction factor. You must set the parameters before
    // using this method.

    Double_t tail;

    if (Mx2 > maxMx2) {
        Haprad();

        sig_obs += tai[0] + tai[1];
        rc = sig_obs / sigma_born;
        sig_obs *= 1E-3;
        sigma_born *= 1E-3;
        tail = tai[1] * 1E-3;

        std::cout << std::endl;
        std::cout.setf(std::ios::fixed);
        std::cout << "    sib    " << std::setw(10)
                                   << std::setprecision(7) << sigma_born
                  << "    sig    " << std::setw(10)
                                   << std::setprecision(7) << sig_obs
                  << "    tail   " << std::setw(10)
                                   << std::setprecision(7) << tail << std::endl;
    } else {
        rc = 0;
    }

    return rc;
}



Double_t TRadCor::GetRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                              Double_t p_t, Double_t phi, Double_t maxMx2)
{
    // Get the radiative correction factor for the given parameters.
    //
    // The E parameter is the energy of the beam, the x, Q2, z, p_t and phi
    // are the values of the kinematical variables who describe the cross
    // section of hadron electroproduction, and maxMx2 is the maximum amount of
    // missing mass.

    SetParameters(E,x,Q2,z,p_t,phi,maxMx2);
    return GetRCFactor();
}



void TRadCor::Haprad(void)
{
    pl = 0.;

    try {
        fInv.Evaluate();
        if (fKin.Y() < 0.) {
            Double_t y = fInv.Q2() / (fInv.S() * fKin.X());
            fKin.SetY(y);
        }
    } catch (TKinematicException& wrongKin) {
        std::cout << wrongKin.what() << std::endl;
        return;
    }

    try {
        fHadKin.Evaluate();
    } catch (TKinematicException& wrongKin) {
        std::cout << wrongKin.what() << std::endl;
        return;
    }


    N = kPi * SQ(kAlpha) * fKin.Y() * fInv.Sx() * M / 2. / fInv.SqrtLq() * kBarn;
#ifdef DEBUG
    std::cout.setf(std::ios::fixed);
    std::cout << "N      " << std::setw(20)
                           << std::setprecision(10) << N << std::endl;
#endif

    if (fKin.T() >= 0.) {
        if (fHadKin.Ph() > fHadKin.Pt()) {
            N = N * fInv.SqrtLq() / 2. / M / fHadKin.Pl();
        } else {
            N = 0.;
        }

        Double_t t = SQ(m_h) - fInv.Q2() + 2. * (fHadKin.SqNuQ() * fHadKin.Pl() - fHadKin.Nu() * fHadKin.Eh());
        fKin.SetT(t);

        std::cout << "p_l: " << fHadKin.Pl() << "\t" << t << "\t"
                  << fHadKin.Pl() - t + fInv.Q2() - SQ(m_h) + 2. * fHadKin.Nu() * fHadKin.Eh() / 2. / fHadKin.SqNuQ()
                  << std::endl;
    }
#ifdef DEBUG
    std::cout << "Pt     " << std::setw(20) << std::setprecision(10)
                           << fHadKin.Pt() << std::endl;
    std::cout << "Pl     " << std::setw(20) << std::setprecision(10)
                           << fHadKin.Pl() << std::endl;
    std::cout << "tdif   " << std::setw(20) << std::setprecision(10)
                           << fKin.T() << std::endl;
#endif

    fInv.EvaluateV12();
    fHadKin.EvaluatePx2();

    t_min = SQ(m_h) - fInv.Q2() + 2. * (fHadKin.SqNuQ() * fHadKin.Ph() - fHadKin.Nu() * fHadKin.Eh());

    Double_t tdmax;
    tdmax = SQ(m_h) - fInv.Q2() + 2. * (- fHadKin.SqNuQ() * fHadKin.Ph() - fHadKin.Nu() * fHadKin.Eh());

    try {
        if ((fKin.T() - t_min) > kEpsMachine || fKin.T() < tdmax) {
            throw TKinematicException();
        }
    } catch (TKinematicException& wrongKin) {
        std::cout << wrongKin.what()
                << " t     = " << fKin.T() << std::endl
                << " t_max = " << tdmax << std::endl
                << " t_min = " << t_min << " - " << fKin.T() - t_min << std::endl;
        return;
    }

    SPhiH();
}



void TRadCor::SPhiH(void)
{
    TDelta fDeltas(this);
    fDeltas.Evaluate();

    std::cout << "********** ita: " << 1 << " *********" << std::endl;
    TBorn fBornin(this);
    sigma_born = fBornin.GetValue(N);
    std::cout << "sib: " << sigma_born << std::endl;
    if (sigma_born == 0.0) {
        tai[0] = 0.;
    }

    qqt(tai); // fix this

    Double_t extai1 = TMath::Exp(fDeltas.Inf());
    sig_obs = sigma_born * extai1 * (1. + fDeltas.VR() + fDeltas.Vac());
}



void TRadCor::BorninTest(Double_t& sigma_born)
{

}



void TRadCor::qqt(Double_t tai[])
{
    TQQTPhi qphi(this);
    if (fConfig.IntegratePhiRad() == 1) {
        ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kNONADAPTIVE);
        ig.SetFunction(qphi);
        ig.SetRelTolerance(fConfig.EpsPhiR());
        ig.SetAbsTolerance(1E-18);
        tai[0] = ig.Integral(0, TMath::TwoPi());
        tai[0] = N * kAlpha * tai[0] / SQ(kPi) / 4. / fInv.SqrtLq();
    } else if (fConfig.IntegratePhiRad() == 0) {
        tai[0] = N * kAlpha / kPi * qphi(0.) / 2 / fInv.SqrtLq();
    }

    std::cout << "tai[" << 0 << "]\t"  << tai[0] << std::endl;

    Double_t tau_1, tau_2;
    Double_t tau[6];
    Double_t phi[4];

    tau_1 = - fInv.Q2() / fInv.S();
    tau_2 =   fInv.Q2() / fInv.X();

    phi[0] = 0.;
    phi[1] = 0.01 * kPi;
    phi[2] = 2. * kPi - 0.01 * kPi;
    phi[3] = 2. * kPi;

    double tau_max = (fInv.Sx() + fInv.SqrtLq()) / (2. * SQ(M));
    double tau_min = - fInv.Q2() / SQ(M) / tau_max;
    tau[0] = tau_min;
    tau[1] = tau_1 - 0.15 * (tau_1 - tau_min);
    tau[2] = tau_1 + 0.15 * (tau_2 - tau_1);
    tau[3] = tau_2 - 0.15 * (tau_2 - tau_1);
    tau[4] = tau_2 + 0.15 * (tau_max - tau_2);
    tau[5] = tau_max;

    std::cout << "********** ita: " << 2 << " *********" << std::endl;

    ROOT::Math::GSLMCIntegrator ig(ROOT::Math::IntegrationMultiDim::kMISER,
                                   1E-6,
                                   1E-3,
                                   20000);
    TRV2TR rv2tr(this);
    ig.SetFunction(rv2tr);

    Double_t rere = 0.;
    for (Int_t iph = 0; iph < 3; iph++) {
        for (Int_t ico = 0; ico < 5; ico++) {
            Double_t lower_lim[2], upper_lim[2];

            lower_lim[0] = tau[ico];
            upper_lim[0] = tau[ico+1];
            lower_lim[1] = phi[iph];
            upper_lim[1] = phi[iph+1];

            if (lower_lim[1] > upper_lim[1])
                std::cout << lower_lim[1] << " < " << upper_lim[1] << std::endl;
            if (lower_lim[0] > upper_lim[0])
                std::cout << lower_lim[0] << " < " << upper_lim[1] << std::endl;

            Double_t re = ig.Integral(lower_lim,upper_lim);

            std::cout.setf(std::ios::fixed);
            std::cout << " tai: "
                        << std::setw(4)  << ico + 1 << "  "
                        << std::setw(4)  << iph + 1 << "  "
                        << std::setw(10) << std::setprecision(4) << re << "  "
                        << std::endl;

            rere = rere + re;
        }
    }

    tai[1] = - kAlpha / (64. * TMath::Power(kPi,5.) * fInv.SqrtLq() * M) * N * rere;
    std::cout << "tai[" << 1 << "]\t"  << tai[1] << std::endl;
}
