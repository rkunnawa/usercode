//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 19 11:35:34 2012 by ROOT version 5.27/06b
// from TTree HltTree/
// found on file: ../HiForest-promptskim-hihighpt-hltjet80-pt90-v2_v3_part2.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Skims {
public :
   Skims(){};
   ~Skims(){};

   // Declaration of leaf types
   Int_t           reco_extra;
   Int_t           reco_extra_jet;
   Int_t           pat_step;
   Int_t           ana_step;
   Int_t           pcollisionEventSelection;
   Int_t           pPAcollisionEventSelectionPA;
   Int_t           pHBHENoiseFilter;
   Int_t           hltAna;
   Int_t           rechitAna;

   // List of branches
   TBranch        *b_reco_extra;   //!
   TBranch        *b_reco_extra_jet;   //!
   TBranch        *b_pat_step;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pPAcollisionEventSelectionPA;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_hltAna;   //!
   TBranch        *b_rechitAna;   //!

};


void setupSkimTree(TTree *t,Skims &tSkims,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("reco_extra", &tSkims.reco_extra, &tSkims.b_reco_extra);
   t->SetBranchAddress("reco_extra_jet", &tSkims.reco_extra_jet, &tSkims.b_reco_extra_jet);
   t->SetBranchAddress("pat_step", &tSkims.pat_step, &tSkims.b_pat_step);
   t->SetBranchAddress("ana_step", &tSkims.ana_step, &tSkims.b_ana_step);
   t->SetBranchAddress("pcollisionEventSelection", &tSkims.pcollisionEventSelection, &tSkims.b_pcollisionEventSelection);
   t->SetBranchAddress("pPAcollisionEventSelectionPA", &tSkims.pPAcollisionEventSelectionPA, &tSkims.b_pPAcollisionEventSelectionPA);
   t->SetBranchAddress("pHBHENoiseFilter", &tSkims.pHBHENoiseFilter, &tSkims.b_pHBHENoiseFilter);
   t->SetBranchAddress("hltAna", &tSkims.hltAna, &tSkims.b_hltAna);
   t->SetBranchAddress("rechitAna", &tSkims.rechitAna, &tSkims.b_rechitAna);
   if (doCheck) {
   }
}

