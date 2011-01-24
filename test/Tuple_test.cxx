/*
 * ===========================================================================
 *
 *       Filename:  Tuple_test.cxx
 *
 *    Description:  Create tuples from HapRad
 *
 *        Version:  1.0
 *        Created:  21/01/11 16:16:09
 *  Last Modified:  24/01/11 15:47:24
 *       Compiler:  gcc
 *
 *         Author:  Ricardo Oyarzun
 *          Email:  royarzun@alumnos.inf.utfsm.cl
 *
 * ===========================================================================
 */

#include <iostream>
#include <exception>
#include "TFile.h"
#include "TRadCor.h"
#include "TNtuple.h"
#include "haprad_constants.h"

void CreateTuple()
{
    TFile *input = new TFile("tupla.root");
    TFile *output = new TFile("Haprad_Tuples.root", "RECREATE");

    TNtuple *t_input = (TNtuple *)input->Get("Tuple-0");
    TNtuple *t_output  = new TNtuple("Tuple", "Tuple from Haprad", "Q2:Xb:Pt:Zh:Phi:f1:f2:f3");

    Float_t x, Q2, pt, z, phi;
    Double_t f1, f2, f3;
    TRadCor rc;
    Double_t m = TMath::Power((kMassNeutron + kMassPion), 2);

    t_input->SetBranchAddress("Q2", &Q2);
    t_input->SetBranchAddress("x", &x);
    t_input->SetBranchAddress("pt", &pt);
    t_input->SetBranchAddress("z", &z);
    t_input->SetBranchAddress("phi", &phi);

    Int_t nEntries = t_input->GetEntries();
    //Int_t nEntries = 1000;
    for (Int_t i = 0; i < nEntries; ++i) {
        t_input->GetEntry(i);
        if (Q2 < 4 && Q2 > 1 && z < 0.7 && z > 0.4 &&
            pt < 4 && pt > 0 && phi > -180 && phi < 180){
            rc.CalculateRCFactor(5.015, x, Q2, z, pt, phi, m);
            f1 = rc.GetFactor1();
            f2 = rc.GetFactor2();
            f3 = rc.GetFactor3();
            t_output->Fill(Q2, x, pt, z, phi, f1, f2, f3);
        }
    }
    input->Close();
    t_output->Write();
    output->Close();
    return;

}

int main()
{

    CreateTuple();
    return 0;

}

