#include "TSemiInclusiveModel.h"
#include "TMath.h"
#include "Partons.h"
#include "ConfigFile.h"
#include "haprad_constants.h"
#include "square_power.h"
#include <iostream>

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "THn.h"

/*
 *This C extern declaration is used for calling the pdf libraries
 *form the Cern Libs directly using the previous fortran code.
 *
 */

extern "C" {
    void init_pdf_(int*, int*);

    void exec_structm_(double*, double*, double*, double*, double*, double*,
                       double*, double*, double*, double*, double*);

    void exec_pkhff_(int*, int*, double*, double*, double*, double*, double*,
                     double*, double*, double*);
}



namespace HapradUtils {

void SemiInclusiveModel(Double_t  q2,  Double_t  X,
                        Double_t  Y,   Double_t  Z,
                        Double_t  pt2, Double_t  mx2,
                        Double_t  pl,
                        Double_t& H1z, Double_t& H2z,
                        Double_t& H3z, Double_t& H4z)
{
    using namespace TMath;

    static Int_t nc = 0;
    nc++;

    static Double_t par0, par1, par2, par3, par4;
    static Double_t A1, B1, C1, D1, E1;
    //static TH3F *H3Hist;
    Int_t bins[4] = {6,5,10,5};
    Double_t xmins[4] = {1.,0.55,0.,0.};
    Double_t xmaxs[4] = {4.,0.1,1.,1.5};
    static THnD H3Hist("h3_hist", "H3 histo", 4, bins, xmins, xmaxs);
    static THnD H4Hist("h4_hist", "H4 histo", 4, bins, xmins, xmaxs);
    Double_t values[4];
    if (nc == 1) {
        ConfigFile parameters("parameters");

        std::cout << "   Reading parameters from configuration file..."
                    << std::endl;
        parameters.readInto(par0, "par0");
        parameters.readInto(par1, "par1");
        parameters.readInto(par2, "par2");
        parameters.readInto(par3, "par3");
        parameters.readInto(par4, "par4");

        parameters.readInto(A1, "A1");
	parameters.readInto(B1, "B1");
	parameters.readInto(C1, "C1");
	parameters.readInto(D1, "D1");
	parameters.readInto(E1, "E1");

	gSystem->Load("libTree");

	TFile file("newphihist.root");
	TNtuple *tuple = (TNtuple*)file.Get("AAcAcc_data");
	Float_t rXb, rQ2, rZh, rPt, rAc, rAcc;
	tuple->SetBranchAddress("Q2",&rQ2);
	tuple->SetBranchAddress("Xb",&rXb);
	tuple->SetBranchAddress("Zh",&rZh);
	tuple->SetBranchAddress("Ac",&rAc);
	tuple->SetBranchAddress("Pt",&rPt);
	tuple->SetBranchAddress("Acc",&rAcc);
	Int_t entries = tuple->GetEntries();
	for (int i = 0; i < tuple->GetEntries(); i++){
		tuple->GetEntry(i);
		values[0] = rQ2; values[1] = rXb;
		values[2] = rZh; values[3] = rPt;
		H3Hist.Fill(values,rAc);
		H4Hist.Fill(values,rAcc);

	}
    }

    Double_t GTMD, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, XD, GL;

    /*Check the kinematics*/
    if (X < 0 || X > 1) return;
    if (Z < 0 || Z > 1) return;


    if (q2 > 1) SCALE = Sqrt(q2);
    else        SCALE = 1.0000001;

    Int_t GPDF      =   5;
    Int_t SPDF      =   5;
    Int_t ISET      =   1;
    Int_t ICHARGE   =   1;

    if (nc == 1) init_pdf_(&GPDF, &SPDF);

    Double_t r  = Sqrt(1. + SQ(2 * kMassProton * X) / q2);
    Double_t xi = 2. * X / (1 + r);

    XD = xi;
    GL = 1.033;

    exec_structm_(&XD, &SCALE, &UPV, &DNV, &USEA, &DSEA, &STR, &CHM, &BOT, &TOP, &GL);

    Double_t    ZD;
    Double_t    Q2D;
    Int_t       IFINI;

    if (Z > 0.01)   ZD = Z;
    else            ZD = 0.01000001;

    if (q2 > 1)     Q2D = q2;
    else            Q2D = 1.0000001;

    if (nc == 1)    IFINI = 0;
    else            IFINI = 1;

    Double_t uff[2], dff[2], sff[2], cff[2], bff[2], gff[2];

    exec_pkhff_(&ISET, &ICHARGE, &ZD, &Q2D, uff, dff, sff, cff, bff, gff);

    Double_t meanpt = par0 + par1 * Z + par2 * Z * Z + par3 * X + par4 * X * X;
    Double_t sgmpt  = A1 + B1 * Z + C1 * Z * Z + D1 * X + E1 * X * X;

    if (sgmpt < 0.02) sgmpt = 0.02;
    if (sgmpt > 0.15) sgmpt = 0.15;

    GTMD = Exp(-SQ(Sqrt(pt2) - meanpt) / (2 * sgmpt)) / (2 * Pi() * sgmpt);

    if (mx2 < SQ(kMassProton + kMassPion)) return;

    Double_t uq, dq, sq, cq, bq, tq, gg;
    Double_t pi_thresh;
    Double_t H1, H2, H3m, H4m;
    Double_t xv, zv, q2v;
    Double_t rlt = 0.14;

    pi_thresh = Sqrt(1. - SQ(kMassProton + kMassPion) / mx2);
    uq = eu2 * ((UPV + USEA) * uff[0] + USEA * uff[1]);
    dq = ed2 * ((DNV + DSEA) * dff[0] + DSEA * dff[1]);
    sq = es2 * (STR * sff[0] + STR * sff[1]);
    cq = ec2 * (CHM * cff[0] + CHM * cff[1]);
    bq = eb2 * (BOT * bff[0] + BOT * bff[1]);
    tq = 0.0;
    gg = 0.0;
    H2 = (uq + dq + sq + cq + bq + tq + gg) * GTMD * pi_thresh;
    H1 = H2 / (2 * X * (1 + rlt)) * (1 + 4 * SQ(kMassProton) * SQ(X) / q2);
    xv = X;
    zv = Z;
    q2v = q2;

    if (X  < 0.1) xv  = 0.1;
    if (q2 < 1.0) q2v = 1.;
    if (Z  < 0.1) zv  = 0.1;

    Double_t Ebeam, rt, rtz, cterm, m_cos_phi, m_cos_2phi;
    Double_t bin[4];
    bin[0] = q2; bin[1] = X; bin[2] = Z; bin[3] = sqrt(pt2);

//    std::cout << q2 << " " << X << " " << X << " " << H3Hist.GetBinContent(H3Hist.GetBin(bin)) << std::endl;
//    std::cout << q2 << " " << X << " " << X << " " << H4Hist.GetBinContent(H4Hist.GetBin(bin))  << std::endl;

    H3m = H3Hist.GetBinContent(H3Hist.GetBin(bin)) * GTMD * pi_thresh;
    H4m = H4Hist.GetBinContent(H4Hist.GetBin(bin)) * GTMD * pi_thresh;
//    H3m = h3(X, q2, Z) * GTMD * pi_thresh;
//    H4m = h4(X, q2, Z) * GTMD * pi_thresh;

    // Check that cos(phi) and cos(2phi) are less than 1
    Ebeam   = q2 / (2 * kMassProton * X * Y);
    rt      = 1 - Y - kMassProton * X * Y / (2 * Ebeam);
    rtz     = Sqrt(rt / (1 + 2 * kMassProton * X / (Y * Ebeam)));
    cterm   = X * SQ(Y) * H1 + rt * H2;

    m_cos_phi = Sqrt(pt2 / q2) * (2 - Y) * rtz * H3m / (2 * cterm);

    if (Abs(4. * m_cos_phi) > 0.9) {
        H3m = 0.9 * Sign(1., m_cos_phi) * (2 * cterm) /
                (Sqrt(pt2 / q2) * (2 - Y) * rtz) / 4.;
    }

    m_cos_2phi = pt2 / q2 * SQ(rtz) * H4m / (2 * cterm);

    if (Abs(4 * m_cos_2phi) > 0.9)
        H4m = 0.9 * Sign(1., m_cos_2phi) * (2 * cterm) / (pt2 / q2 * SQ(rtz)) / 4.;

    H1z = H1;
    H2z = H2;
    H3z = H3m;
    H4z = H4m;

    return;
}

} // End namespace HapradUtils



Double_t h3(Double_t X, Double_t q2, Double_t Z)
{

    using namespace TMath;
    // Initialized data
    Double_t q0     =    1.0;
    Double_t lambda =    0.25;
    Double_t a      =   -3.6544E-4;
    Double_t a1     =   -2.1855;
    Double_t a2     =    3.4176;
    Double_t b1     =   -1.7567;
    Double_t b2     =    1.1272;
    Double_t bb     =    8.9985;

    if (q2 > q0) {
        return a * Power(X, a1) *
                   Power(1 - X, a2) *
                   Power(Z, b1) *
                   Power(1 - Z, b2) *
                   Power(Log(q2 / SQ(lambda)) / Log(q0 / SQ(lambda)), bb);
    } else {
        return a * Power(X, a1) *
                   Power(1 - X, a2) *
                   Power(Z, b1) *
                   Power(1 - Z, b2);
    }
}



Double_t h4(Double_t X, Double_t q2, Double_t Z)
{
    using namespace TMath;

    // Initialized data
    Double_t q0     =   1.0;
    Double_t lambda =   0.25;
    Double_t a      =   0.0010908;
    Double_t a1     =  -3.5265E-7;
    Double_t a2     =   3.0276E-8;
    Double_t b1     =  -0.66787;
    Double_t b2     =   3.5868;
    Double_t bb     =   6.8777;

    if (q2 > q0) {
        return a * Power(X, a1) *
                   Power(1 - X, a2) *
                   Power(Z, b1) *
                   Power(1 - Z, b2) *
                   Power(Log(q2 / SQ(lambda)) / Log(q0 / SQ(lambda)), bb / X);
    } else {
        return a * Power(X, a1) *
                   Power(1 - X, a2) *
                   Power(Z, b1) *
                   Power(1 - Z, b2);
    }
}
