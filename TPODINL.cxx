#include "TPODINL.h"
#include "TRadCor.h"
#include "TKinematicalVariables.h"
#include "TLorentzInvariants.h"
#include "TStructFunctionArray.h"
#include "TThetaMatrix.h"
#include "haprad_constants.h"
#include "square_power.h"
#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif


TPODINL::TPODINL(const TRadCor* rc, double tau, double mu,
                 const TStructFunctionArray& H0, const TThetaMatrix& theta)
  : fTau(tau), fMu(mu), fH0(H0), fTheta(theta)
{
    fRC     = rc;
    fInv    = rc->GetLorentzInvariants();
}



TPODINL::~TPODINL()
{

}



ROOT::Math::IBaseFunctionOneDim* TPODINL::Clone() const
{
    return 0;
}



double TPODINL::DoEval(double R) const
{
    TStructFunctionArray H(fRC);
    H.Evaluate(fTau, fMu, R);

    double pp = 0.;
    double pres = 0.;
    double podinl = 0.;
    for (int isf = 0 ; isf < 4; ++isf) {
        for (int irr = 0; irr < 3; ++irr) {
            pp = H[isf];
            if (irr == 0) pp = pp - fH0[isf] * SQ(1 + R * fTau / fInv->Q2());

            pres = pp * TMath::Power(R, (irr - 1)) / SQ(fInv->Q2() + R * fTau);
            podinl = podinl - fTheta[isf][irr] * pres;
        }
    }
    return podinl;
}
