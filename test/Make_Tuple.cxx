/*
 * ===========================================================================
 *
 *       Filename:  Make_Tupla.cxx
 *
 *    Description:  Makes a TNtuple root file
 *
 *        Version:  1.0
 *        Created:  26/01/11 09:23:46
 *  Last Modified:  26/01/11 11:55:37
 *       Compiler:  gcc
 *
 *         Author:  Ricardo Oyarzun
 *          Email:  royarzun@alumnos.inf.utfsm.cl
 *
 * ===========================================================================
 */
#include<iostream>
#include "TNtuple.h"
#include "TFile.h"
#include "haprad_constants.h"
#include "TRadCor.h"

void MakeTuple(){
    TNtuple *t_output = new TNtuple("Tuple","Tuple from HapradCpp", "Q2:Xb:Pt:Zh:Phi:f1:f2:f3");
    Double_t q2Min,q2Max,Nu_min,Nu_max,pt2_min,pt2_max,z_min,z_max;
    Double_t phi_min,phi_max;
    Double_t M = TMath::Power((kMassNeutron + kMassPion),2);
    
    TRadCor rc;

    q2Min = 1.;
    q2Max = 4.;
    Nu_min = 2.2;
    Nu_max = 4.;
    pt2_min = 0.;
    pt2_max = 4.;
    z_min = 0.4;
    z_max = 0.7;
    phi_min = -3.1415926;
    phi_max = 3.1415926;

    Int_t Nu_n = 3;
    Int_t Zh_n = 3;
    Int_t Pt2_n = 20;
    Int_t Q2_n = 3;
    Int_t Phi_n = 12;

    Double_t Q2,Xb,Pt,Zh,Phi,f1,f2,f3,Nu;
    TFile *f = new TFile("RC_Tuples.root","RECREATE");
    Int_t Nu_i = 0;
    while(Nu_i < Nu_n){
        Nu = Nu_min + (double)((double)Nu_i - 0.5) * (Nu_max - Nu_min)/Nu_n;
        Int_t Q2_i = 0;
        while(Q2_i < Q2_n){
            Q2 = q2Min + (double)((double)Q2_i - 0.5)*(q2Max - q2Min)/Q2_n;
            Xb = Q2/2./Nu/kMassUndetectedHadron;
            Int_t Zh_i = 0;
            while(Zh_i < Zh_n){
                Zh = z_min + (double)((double)(Zh_i) - 0.5) * (z_max - z_min)/Zh_n;
                Int_t Pt2_i = 0;
                while(Pt2_i < Pt2_n){
                    Pt = pt2_min + (double)((double)(Pt2_i) - 0.5) * (pt2_max - pt2_min)/Pt2_n;
                    Int_t Phi_i = 0;
                    while(Phi_i < Phi_n){
                        Phi = phi_min + (double)((double)(Phi_i) - 0.5) * (phi_max - phi_min)/Phi_n;
                        rc.CalculateRCFactor(5.015,Xb,Q2,Zh,Pt,Phi,M);
                        f1 = rc.GetFactor1();
                        f2 = rc.GetFactor2();
                        f3 = rc.GetFactor3();
                        t_output->Fill(Q2,Xb,Pt,Zh,Phi,f1,f2,f3);
                        Phi_i++;
                    }
                    Pt2_i++;
                }
                Zh_i++;
            }
            Q2_i++;
        }
        Nu_i++;
    }
    t_output->Write();
    delete t_output;
    f->Close();

}


int main(){
    MakeTuple();
    return 0;
}


