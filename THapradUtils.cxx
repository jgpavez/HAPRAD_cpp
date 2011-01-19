#include "THapradUtils.h"
#include "haprad_constants.h"
#include "square_power.h"
#include "Partons.h"
#include <iostream>
#include <string>

#include "TMath.h"
#include "TROOT.h"


namespace HapradUtils {

    using namespace TMath;

    Double_t fspen(const Double_t x)
    {
        const Double_t f1 = 1.644934;
        Double_t result;

        if (x < -1) {
            Double_t logprod;
            logprod = TMath::Log(1 - x) * TMath::Log(SQ(x) / (1 - x));
            result = - 0.5 * logprod - f1 + fspens(1 / (1 - x));
        } else if (x < 0) {
            Double_t log2;
            log2 = TMath::Log(1 - x) * TMath::Log(1 - x);
            result = - 0.5 * log2 - fspens(x / (x - 1));
        } else if (x <= 0.5) {
            result = fspens(x);
        } else if (x <= 1) {
            Double_t logprod;
            logprod = TMath::Log(x) * TMath::Log(1 - x + 1E-10);
            result = f1 - logprod - fspens(1 - x);
        } else if (x <= 2) {
            Double_t logprod;
            logprod = TMath::Log(x) * TMath::Log((x - 1) * (x - 1) / x);
            result = f1 - 0.5  * logprod + fspens(1 - 1 / x);
        } else {
            Double_t log2;
            log2 = TMath::Log(x) * TMath::Log(x);
            result = 2 * f1 - 0.5 * log2 - fspens(1 / x);
        }

        return result;
    }



    Double_t fspens(const Double_t x)
    {
        Double_t f   = 0;
        Double_t a   = 1;
        Double_t an  = 0;
        Double_t tch = 1E-16;
        Double_t b;

        do {
            an = an + 1;
            a = a * x;
            b = a / SQ(an);
            f = f + b;
        } while (b - tch > 0);

        return f;
    }



    Double_t dfint(Int_t narg, double *arg, Int_t *nent, Double_t *ent, Double_t *table)
    {
        Double_t kd = 0;
        Double_t m = 1;
        Int_t ja = 0;
        Int_t jb, k;
        Int_t ncomb[narg];
        Double_t d[narg];
        Int_t ient[narg];
        Double_t dfint;
        Int_t i = 0;
        while (i < narg) {
            ncomb[i] = 1;
            jb = ja + nent[i]; // VERIFICAR SI LOS LIMITES ESTAN CORRECTOS!!!!
            Int_t j = ja;
            while (j <= jb) {
                if (arg[i] <= ent[j]) goto exit_inner_loop;
                j++;
            }
            j = jb;
    exit_inner_loop:
            if (j == ja) ++j;
            d[i] = (ent[j] - arg[i]) / (ent[j] - ent[j - 1]);
            ient[i] = j - ja;
            kd = kd + ient[i] * m;
            m = m * nent[i];
            ja = jb + 1;
            i++;
        }

        dfint = 0;
    do_again:
        Double_t fac = 1.;
        Int_t iadr = kd;
        Int_t ifadr = 1;
        //cambiado a 0 dado el indice de los arreglos en fortran
        i = 0;
        while (i < narg) {
            if (ncomb[i] == 0) {
                fac = fac * d[i];
                iadr  = iadr - ifadr;
            } else {
                fac = fac * (1. -  d[i]);
            }
            ifadr = ifadr * nent[i];
            ++i;
        }

        dfint = dfint + fac * table[iadr];
        Int_t il = narg - 1;
        while (ncomb[il] == 0) {
            --il;
            if (il < 0) return dfint;
        }

        ncomb[il] = 0;
        if (il == narg) goto do_again;
        k = il + 1;
        while (k < narg) {
            ncomb[k] = 1;
            k++;
        }
        goto do_again;
    }
}
