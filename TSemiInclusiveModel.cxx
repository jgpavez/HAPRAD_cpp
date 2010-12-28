#include "TSemiInclusiveModel.h"
#include "TMath.h"
#include "Partons.h"
#include "haprad_constants.h"

/*
 *This C extern declaration is used for calling the pdf libraries
 *form the Cern Libs directly using the previous fortran code.
 *
 */

extern "C" {
    void init_pdf_(int&, int&);

    void exec_structm_(double&, double&, double&, double&, double&, double&,
                       double&, double&, double&, double&, double&);

    void exec_pkhff_(int&, int&, double&, double&, double*, double*, double*,
                     double*, double*, double*);
}

namespace HapradUtils {

    void SemiInclusiveModel(Double_t q2, Double_t X,
                            Double_t Y, Double_t Z,
                            Double_t pt2, Double_t mx2,
                            Double_t pl, Double_t& H1z,
                            Double_t& H2z, Double_t& H3z,
                            Double_t& H4z)
    {
        using namespace TMath;

        Double_t ac = 1.2025 * Power(10, -10);
        Double_t bc = -5.2703 * Power(10, -2);
        Double_t cc = 3.7467 * Power(10, -1);
        Double_t dc = 6.5397 * Power(10, -2);
        Double_t ec = -2.2136 * Power(10, -1);
        Double_t fc = -1.0621 * Power(10, -1);
        Double_t GTMD, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, XD, GL;
        Double_t r, xi;
        Int_t GPDF = 5;
        Int_t SPDF = 5;
        Int_t ISET = 1;
        Int_t ICHARGE = 1;
        Int_t nc = 0;
        /*Check the kinematics*/
        if (X < 0 || X > 1) return;

        if (Z < 0 || Z > 1) return;

        nc++;
        r = Sqrt(1. + Power(2 * kMassProton * X, 2) / q2);
        xi = 2.*X / (1 + r);

        if (q2 > 1) SCALE = Sqrt(q2);
        else SCALE = 1.0000001;

        if (nc == 1) init_pdf_(GPDF, SPDF);

        GL = 1.033;

        exec_structm_(XD, SCALE, UPV, DNV, USEA, DSEA, STR, CHM, BOT, TOP, GL);

        Double_t ZD;
        if (Z > 0.01) ZD = Z;
        else ZD = 0.01000001;

        Double_t Q2D;
        if (q2 > 1) Q2D = q2;
        else Q2D = 1.0000001;

        Int_t IFINI;
        if (nc == 1) IFINI = 0;
        else IFINI = 1;

        Double_t uff[2], dff[2], sff[2], cff[2], bff[2], gff[2];

        exec_pkhff_(ISET, ICHARGE, ZD, Q2D, uff, dff, sff, cff, bff, gff);

        Double_t sgmpt = ac + bc * X + cc * Z + dc * Power(Z, 2) +
                                                ec * Power(Z, 2) + fc * X * Z;

        if (sgmpt < 0.02) sgmpt = 0.02;
        if (sgmpt > 0.15) sgmpt = 0.15;

        if (pl > 0.15)
            GTMD  = Exp(-pt2 / (2.*sgmpt)) / (2.* kPi * sgmpt);
        else
            GTMD  = Exp(-(pt2 + 2.*Power(pl, 2)) / (2.*sgmpt)) / (2.* kPi * sgmpt);

        if (mx2 < Power((kMassProton + kMassPion), 2)) return;

        Double_t uq, dq, sq, cq, bq, tq, gg;
        Double_t pi_thresh;
        Double_t H1, H2, H3m, H4m;
        Double_t xv, zv, q2v;
        Double_t rlt = 0.14;

        pi_thresh = Sqrt(1. - Power((kMassProton + kMassPion), 2) / mx2);
        uq = eu2 * ((UPV + USEA) * uff[0] + USEA * uff[1]);
        dq = ed2 * ((DNV + DSEA) * dff[0] + DSEA * dff[1]);
        sq = es2 * (STR * sff[0] + STR * sff[1]);
        cq = ec2 * (CHM * cff[0] + CHM * cff[1]);
        bq = eb2 * (BOT * bff[0] + BOT * bff[1]);
        tq = 0.0;
        gg = 0.0;
        H2 = (uq + dq + sq + cq + bq + tq + gg) * GTMD * pi_thresh;
        H1 = H2 / (2. * X * (1. + rlt)) *
                  (1. + 4. * Power(kMassProton, 2) * Power(X, 2) / q2);
        xv = X;
        zv = Z;
        q2v = q2;

        if (X < 0.1) xv = 0.1;
        if (q2 < 1.) q2v = 1.;
        if (Z < 0.1) zv = 0.1;

        Double_t Ebeam, rt, rtz, cterm, m_cos_phi, m_cos_2phi;

        H3m = h3(X, q2, Z) * GTMD * pi_thresh;
        H4m = h4(X, q2, Z) * GTMD * pi_thresh;

// Check that <cos(phi)> and <cos(2phi)> are < 1

        Ebeam = q2 / (2. * kMassProton * X * Y);
        rt = 1. - Y - kMassProton * X * Y / (2. * Ebeam);
        rtz = Sqrt(rt / (1. + 2. * kMassProton * X / (Y * Ebeam)));
        cterm = X * Power(Y, 2) * H1 + rt * H2;
        m_cos_phi = Sqrt(pt2 / q2) * (2. - Y) * rtz * H3m / (2.*cterm);

        if (Abs(4. * m_cos_phi) > 0.9) {
            H3m = 0.9 * Sign(1., m_cos_phi) * (2.*cterm) /
                    (Sqrt(pt2 / q2) * (2. - Y) * rtz) / 4.;
        }

        m_cos_2phi = pt2 / q2 * Power(rtz, 2) * H4m / (2.*cterm);

        if (Abs(4.*m_cos_2phi) > 0.9)
            H4m = 0.9 * Sign(1., m_cos_2phi) * (2.*cterm) /
                    (pt2 / q2 * Power(rtz, 2)) / 4.;

        H1z = H1;
        H2z = H2;
        H3z = H3m;
        H4z = H4m;

        return;
    }
}



Double_t h3(Double_t X, Double_t q2, Double_t Z)
{
    using namespace TMath;
    /* InitialiZed data */
    Double_t q0 = 1.;
    Double_t lambda = .25;
    Double_t a = -3.6544e-4;
    Double_t a1 = -2.1855;
    Double_t a2 = 3.4176;
    Double_t b1 = -1.7567;
    Double_t b2 = 1.1272;
    Double_t bb = 8.9985;
    if (q2 > q0) {
        return a * Power(X, a1) * Power((1 - X), a2) * Power(Z, b1) *
                   Power((1 - Z), b2) *
                   Power((Log(q2 / Power(lambda, 2))) /
                                   Power(Log(q0 / Power(lambda, 2)), 2), bb);
    } else {
        return a * Power(X, a1) * Power((1 - X), a2) *
                   Power(Z, b1) * Power((1 - Z), b2);
    }

}



Double_t h4(Double_t X, Double_t q2, Double_t Z)
{
    using namespace TMath;

    /* InitialiZed data */
    Double_t q0 = 1.;
    Double_t lambda = .25;
    Double_t a = .0010908;
    Double_t a1 = -3.5265e-7;
    Double_t a2 = 3.0276e-8;
    Double_t b1 = -.66787;
    Double_t b2 = 3.5868;
    Double_t bb = 6.8777;

    if (q2 > q0) {
        return a * Power(X, a1) * Power((1. - X), a2) * Power(Z, b1) *
                   Power((1. - Z), b2) *
                   Power(Log(q2 / Power(lambda, 2)) /
                         Log(q0 / Power(lambda, 2)), (bb / X));
    } else {
        return a * Power(X, a1) * Power((1. - X), a2) *
                   Power(Z, b1) * Power((1. - Z), b2);
    }
}
