#ifndef TRADCOR_H
#define TRADCOR_H

#include "TROOT.h"

class THapradConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TRadCor {
public:
    TRadCor();
    ~TRadCor();

    void        CalculateRCFactor(Double_t E, Double_t x, Double_t Q2,
                                  Double_t z, Double_t p_t, Double_t phi, Double_t maxMx2);
    Double_t    GetRCFactor(Double_t E, Double_t x, Double_t Q2, Double_t z,
                            Double_t p_t, Double_t phi, Double_t maxMx2);


    Double_t    GetFactor1(void);
    Double_t    GetFactor2(void);
    Double_t    GetFactor3(void);

    void        RegisteredLepton(Int_t type = 1);
    void        IntegratePhiRad(Int_t type = 0);
    void        IntegratePhiHad(Int_t type = 0);
    void        SetPolarization(Int_t type = 0);

    const THapradConfig*            GetConfig(void) const { return fConfig; };
    const TKinematicalVariables*    GetKinematicalVariables(void) const { return fKin; };
    const TLorentzInvariants*       GetLorentzInvariants(void) const { return fInv; };
    const THadronKinematics*        GetHadronKinematics(void) const { return fHadKin; };

private:
    void        CreateVariables(Double_t E, Double_t x, Double_t Q2,
                                Double_t z, Double_t p_t, Double_t phi);
    void        DeleteVariables(void);

    void        Initialization(void);
    void        SPhiH(void);
    void        qqt(Double_t tai[]);

    THapradConfig*          fConfig;
    TKinematicalVariables*  fKin;
    TLorentzInvariants*     fInv;
    THadronKinematics*      fHadKin;

    Double_t     E;        // The beam energy

    //  Kinematic variables
    Double_t     t_min;
    Double_t     t_max;

    // Results
    Double_t     rc;
    Double_t     sigma_born;    // sigma_0
    Double_t     sig_obs;       // sigma_{obs}
    Double_t     tail;
    Double_t     tai[3];

    // Masses
    Double_t     M;
    Double_t     m;
    Double_t     m_h;

    // Integration
    Double_t     N;                // Normalization factor
    Double_t     pl;
    Int_t        ita;
    Double_t     phi_rad;
};

#endif
