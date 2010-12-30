#ifndef TRV2TR_H
#define TRV2TR_H

#include "Math/IFunction.h"

class TRadCor;
class TGlobalConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TRV2TR : public ROOT::Math::IBaseFunctionMultiDim {
public:
    TRV2TR(const TRadCor* rc);
    ~TRV2TR();

    virtual IBaseFunctionMultiDim* Clone() const;
    virtual unsigned int NDim() const;

private:
    virtual double DoEval(const double *x) const;

    const TRadCor*                  fRC;
    const TGlobalConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;

    double M;
    double M2;
    double mh2;
    double mu2;

    double tau_max;
    double tau_min;
};

#endif
