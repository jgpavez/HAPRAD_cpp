#ifndef THADRONKINEMATICS_H
#define THADRONKINEMATICS_H

#include "TROOT.h"

class THapradConfig;
class TKinematicalVariables;
class TLorentzInvariants;


class THadronKinematics {
public:
    THadronKinematics(const THapradConfig* config,
                      const TKinematicalVariables* kin,
                      const TLorentzInvariants* inv);
    ~THadronKinematics();

    // Getters
    Double_t    Eh(void)  const { return fEh; };
    Double_t    Pl(void)  const { return fPl; };
    Double_t    Pt(void)  const { return fPt; };
    Double_t    Nu(void)  const { return fNu; };

    Double_t    SqNuQ(void)  const { return fNu; };

    Double_t    Px2(void) const { return fPx2; };
    Double_t    Ph(void)  const { return fPh; };
    Double_t    V1(void)  const { return fV1; };
    Double_t    V2(void)  const { return fV2; };

    // Setters
    void    SetEh(void);
    void    SetNu(void);
    void    SetSqNuQ(void);
    void    SetMomentum(void);
    void    SetPx2(void);
    void    SetV12(void);

private:
    const THapradConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;

    Double_t    fEh;
    Double_t    fPl;
    Double_t    fPt;
    Double_t    fNu;

    Double_t    fSqNuQ;

    Double_t    fPx2;
    Double_t    fPh;

    Double_t    fV1;
    Double_t    fV2;
};

#endif
