#ifndef THADRONKINEMATICS_H
#define THADRONKINEMATICS_H

#include "TROOT.h"

class TRadCor;
class TGlobalConfig;
class TKinematicalVariables;
class TLorentzInvariants;


class THadronKinematics {
public:
    THadronKinematics(const TRadCor* rc);
    ~THadronKinematics();

    // Getters
    Double_t    Eh(void)  const { return fEh; };
    Double_t    Pl(void)  const { return fPl; };
    Double_t    Pt(void)  const { return fPt; };
    Double_t    Nu(void)  const { return fNu; };

    Double_t    SqNuQ(void)  const { return fNu; };

    Double_t    Px2(void) const { return fPx2; };
    Double_t    Ph(void)  const { return fPh; };

    // Setters
    void    SetEh(void);
    void    SetNu(void);
    void    SetSqNuQ(void);
    void    SetMomentum(void);
    void    SetPx2(void);

private:
    const TGlobalConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;

    Double_t    fEh;
    Double_t    fPl;
    Double_t    fPt;
    Double_t    fNu;

    Double_t    fSqNuQ;

    Double_t    fPx2;
    Double_t    fPh;
};

#endif
