//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 19 11:19:47 2012 by ROOT version 5.27/06b
// from TTree HiTree/
// found on file: ../Pythia30_HydjetDrum_mix01_HiForest2_v19.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Evts {
public :
   Evts(){};
   ~Evts(){};

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiZDC;
   Float_t         hiZDCplus;
   Float_t         hiZDCminus;
   Float_t         hiHFhit;
   Float_t         hiHFhitPlus;
   Float_t         hiHFhitMinus;
   Float_t         hiET;
   Float_t         hiEE;
   Float_t         hiEB;
   Float_t         hiEEplus;
   Float_t         hiEEminus;
   Int_t           hiNpix;
   Int_t           hiNpixelTracks;
   Int_t           hiNtracks;
   Int_t           hiNtracksPtCut;
   Int_t           hiNtracksEtaCut;
   Int_t           hiNtracksEtaPtCut;
   Int_t           hiNevtPlane;
   Float_t         hiEvtPlanes[76];   //[hiNevtPlane]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiZDC;   //!
   TBranch        *b_hiZDCplus;   //!
   TBranch        *b_hiZDCminus;   //!
   TBranch        *b_hiHFhit;   //!
   TBranch        *b_hiHFhitPlus;   //!
   TBranch        *b_hiHFhitMinus;   //!
   TBranch        *b_hiET;   //!
   TBranch        *b_hiEE;   //!
   TBranch        *b_hiEB;   //!
   TBranch        *b_hiEEplus;   //!
   TBranch        *b_hiEEminus;   //!
   TBranch        *b_hiNpix;   //!
   TBranch        *b_hiNpixelTracks;   //!
   TBranch        *b_hiNtracks;   //!
   TBranch        *b_hiNtracksPtCut;   //!
   TBranch        *b_hiNtracksEtaCut;   //!
   TBranch        *b_hiNtracksEtaPtCut;   //!
   TBranch        *b_hiNevtPlane;   //!
   TBranch        *b_hiEvtPlanes;   //!

};


void setupEvtTree(TTree *t,Evts &tEvts,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("run", &tEvts.run, &tEvts.b_run);
   t->SetBranchAddress("evt", &tEvts.evt, &tEvts.b_evt);
   t->SetBranchAddress("lumi", &tEvts.lumi, &tEvts.b_lumi);
   t->SetBranchAddress("vx", &tEvts.vx, &tEvts.b_vx);
   t->SetBranchAddress("vy", &tEvts.vy, &tEvts.b_vy);
   t->SetBranchAddress("vz", &tEvts.vz, &tEvts.b_vz);
   t->SetBranchAddress("hiBin", &tEvts.hiBin, &tEvts.b_hiBin);
   t->SetBranchAddress("hiHF", &tEvts.hiHF, &tEvts.b_hiHF);
   t->SetBranchAddress("hiHFplus", &tEvts.hiHFplus, &tEvts.b_hiHFplus);
   t->SetBranchAddress("hiHFminus", &tEvts.hiHFminus, &tEvts.b_hiHFminus);
   t->SetBranchAddress("hiZDC", &tEvts.hiZDC, &tEvts.b_hiZDC);
   t->SetBranchAddress("hiZDCplus", &tEvts.hiZDCplus, &tEvts.b_hiZDCplus);
   t->SetBranchAddress("hiZDCminus", &tEvts.hiZDCminus, &tEvts.b_hiZDCminus);
   t->SetBranchAddress("hiHFhit", &tEvts.hiHFhit, &tEvts.b_hiHFhit);
   t->SetBranchAddress("hiHFhitPlus", &tEvts.hiHFhitPlus, &tEvts.b_hiHFhitPlus);
   t->SetBranchAddress("hiHFhitMinus", &tEvts.hiHFhitMinus, &tEvts.b_hiHFhitMinus);
   t->SetBranchAddress("hiET", &tEvts.hiET, &tEvts.b_hiET);
   t->SetBranchAddress("hiEE", &tEvts.hiEE, &tEvts.b_hiEE);
   t->SetBranchAddress("hiEB", &tEvts.hiEB, &tEvts.b_hiEB);
   t->SetBranchAddress("hiEEplus", &tEvts.hiEEplus, &tEvts.b_hiEEplus);
   t->SetBranchAddress("hiEEminus", &tEvts.hiEEminus, &tEvts.b_hiEEminus);
   t->SetBranchAddress("hiNpix", &tEvts.hiNpix, &tEvts.b_hiNpix);
   t->SetBranchAddress("hiNpixelTracks", &tEvts.hiNpixelTracks, &tEvts.b_hiNpixelTracks);
   t->SetBranchAddress("hiNtracks", &tEvts.hiNtracks, &tEvts.b_hiNtracks);
   t->SetBranchAddress("hiNtracksPtCut", &tEvts.hiNtracksPtCut, &tEvts.b_hiNtracksPtCut);
   t->SetBranchAddress("hiNtracksEtaCut", &tEvts.hiNtracksEtaCut, &tEvts.b_hiNtracksEtaCut);
   t->SetBranchAddress("hiNtracksEtaPtCut", &tEvts.hiNtracksEtaPtCut, &tEvts.b_hiNtracksEtaPtCut);
   t->SetBranchAddress("hiNevtPlane", &tEvts.hiNevtPlane, &tEvts.b_hiNevtPlane);
   t->SetBranchAddress("hiEvtPlanes", tEvts.hiEvtPlanes, &tEvts.b_hiEvtPlanes);
   if (doCheck) {
      if (t->GetMaximum("hiNevtPlane")>76) cout <<"FATAL ERROR: Arrary size of hiNevtPlane too small!!!  "<<t->GetMaximum("hiNevtPlane")<<endl;
   }
}

