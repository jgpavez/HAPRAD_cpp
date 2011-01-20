#include <stdlib.h>
#include "haprad_constants.h"
#include "TRadCor.h"
#include <iostream>


int main(int argc, char *argv[])
{

    TRadCor rc;
    Double_t m = TMath::Power((kMassNeutron + kMassPion),2);
    rc.SetParameters(6,0.32,2.5,0.1, -4.6002448552,140,m);
    Double_t factor = rc.GetRCFactor();

    std::cout << std::endl << "RC factor: " << factor << std::endl;
    return 0;
}
