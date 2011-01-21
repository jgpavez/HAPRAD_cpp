/*
 * ===========================================================================
 *
 *       Filename:  Tuple_test.cxx
 *
 *    Description:  Create tuples from HapRad
 *
 *        Version:  1.0
 *        Created:  21/01/11 16:16:09
 *  Last Modified:  21/01/11 17:55:48
 *       Compiler:  gcc
 *
 *         Author:  Ricardo Oyarzun
 *          Email:  royarzun@alumnos.inf.utfsm.cl
 *
 * ===========================================================================
 */

#include <iostream>
#include "TFile.h"
#include "TRadCor.h"
#include "TNtuple.h"
#include "haprad_constants.h"

void CreateTuple(){
    TFile *input = new TFile("tupla.root");
    TFile *output = new TFile("Haprad_Tuples.root","RECREATE");

    TNtuple *t_input = (TNtuple *)input->Get("Tuple-0");
    TNtuple *t_output  = new TNtuple("Tuple","Tuple from Haprad","Q2:Xb:Pt:Zh:Phi:f1:f2:f3");
    

    Double_t Q2,Xb,Pt,Zh,Phi,f1,f2,f3;
    TRadCor rc;
    Double_t m = TMath::Power((kMassNeutron + kMassPion),2);

    t_input->SetBranchAddress("Q2", &Q2);
    t_input->SetBranchAddress("x", &Xb);
    t_input->SetBranchAddress("pt", &Pt);
    t_input->SetBranchAddress("z", &Zh);
    t_input->SetBranchAddress("phi", &Phi);

    Int_t nEntries = t_input->GetEntries();
    Int_t i;
    for(i = 0; i < nEntries;++i){
        t_input->GetEntry(i);    
        rc.CalculateRCFactor(6,Xb,Q2,Zh,Pt,Phi,m);
        f1 = rc.GetFactor1();
        f2 = rc.GetFactor2();
        f3 = rc.GetFactor3();
        t_output->Fill(Q2,Xb,Pt,Zh,Phi,f1,f2,f3);
        std::cout << Q2 <<" "<< Xb <<" "<< Pt <<" "<< Zh <<" "<< Phi <<" "<<f1<<" "<<f2<<" "<<f3<<" "<< std::endl;
    }
    input->Close();
    t_output->Write();
    output->Close();
    return;

}

int main(){

    CreateTuple();
    return 0;

}

