#ifndef TTHETAMATRIX_H
#define TTHETAMATRIX_H

#include "TMatrixD.h"

class TRadCor;
class TGlobalConfig;
class TKinematicalVariables;
class TLorentzInvariants;
class THadronKinematics;


class TThetaMatrix {
public:
    TThetaMatrix(const TRadCor* rc);
    TThetaMatrix(Int_t rows, Int_t cols, const TRadCor* rc);
    ~TThetaMatrix();

    void        Evaluate(Double_t tau, Double_t mu,
                         Int_t ita, Double_t phi_k);

    Double_t&   operator() (unsigned row, unsigned col)
                        { return fData[fCols * row + col]; };
    Double_t    operator() (unsigned row, unsigned col) const
                        { return fData[fCols * row + col]; };

private:
    const TGlobalConfig*            fConfig;
    const TKinematicalVariables*    fKin;
    const TLorentzInvariants*       fInv;
    const THadronKinematics*        fHadKin;

    Int_t       fRows;
    Int_t       fCols;
    Double_t*   fData;
};

#endif
