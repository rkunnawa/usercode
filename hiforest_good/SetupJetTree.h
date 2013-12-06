//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 19 12:25:44 2012 by ROOT version 5.27/06b
// from TTree t/akPu3PFpatJets Jet Analysis Tree
// found on file: ../Pythia30_HydjetDrum_mix01_HiForest2_v19.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Jets {
public :
   Jets(){};
   ~Jets(){};

   // Declaration of leaf types
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[450];   //[nref]
   Float_t         jtpt[450];   //[nref]
   Float_t         jteta[450];   //[nref]
   Float_t         jty[450];   //[nref]
   Float_t         jtphi[450];   //[nref]
   Float_t         jtpu[450];   //[nref]
   Float_t         jtm[450];   //[nref]
   Float_t         discr_fr01[450];   //[nref]
   Float_t         trackMax[450];   //[nref]
   Float_t         trackSum[450];   //[nref]
   Int_t           trackN[450];   //[nref]
   Float_t         trackHardSum[450];   //[nref]
   Int_t           trackHardN[450];   //[nref]
   Float_t         chargedMax[450];   //[nref]
   Float_t         chargedSum[450];   //[nref]
   Int_t           chargedN[450];   //[nref]
   Float_t         chargedHardSum[450];   //[nref]
   Int_t           chargedHardN[450];   //[nref]
   Float_t         photonMax[450];   //[nref]
   Float_t         photonSum[450];   //[nref]
   Int_t           photonN[450];   //[nref]
   Float_t         photonHardSum[450];   //[nref]
   Int_t           photonHardN[450];   //[nref]
   Float_t         neutralMax[450];   //[nref]
   Float_t         neutralSum[450];   //[nref]
   Int_t           neutralN[450];   //[nref]
   Float_t         eMax[450];   //[nref]
   Float_t         eSum[450];   //[nref]
   Int_t           eN[450];   //[nref]
   Float_t         muMax[450];   //[nref]
   Float_t         muSum[450];   //[nref]
   Int_t           muN[450];   //[nref]
   Float_t         matchedPt[450];   //[nref]
   Float_t         matchedR[450];   //[nref]
   Int_t           beamId1;
   Int_t           beamId2;
   Float_t         pthat;
   Float_t         refpt[450];   //[nref]
   Float_t         refeta[450];   //[nref]
   Float_t         refy[450];   //[nref]
   Float_t         refphi[450];   //[nref]
   Float_t         refdphijt[450];   //[nref]
   Float_t         refdrjt[450];   //[nref]
   Float_t         refparton_pt[450];   //[nref]
   Int_t           refparton_flavor[450];   //[nref]
   Int_t           refparton_flavorForB[450];   //[nref]
   Int_t           ngen;
   Int_t           genmatchindex[26];   //[ngen]
   Float_t         genpt[26];   //[ngen]
   Float_t         geneta[26];   //[ngen]
   Float_t         geny[26];   //[ngen]
   Float_t         genphi[26];   //[ngen]
   Float_t         gendphijt[26];   //[ngen]
   Float_t         gendrjt[26];   //[ngen]
   Int_t           gensubid[26];   //[ngen]
   Float_t         smpt[450];   //[nref]

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_discr_fr01;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_trackSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_trackHardSum;   //!
   TBranch        *b_trackHardN;   //!
   TBranch        *b_chargedMax;   //!
   TBranch        *b_chargedSum;   //!
   TBranch        *b_chargedN;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_photonMax;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_photonN;   //!
   TBranch        *b_photonHardSum;   //!
   TBranch        *b_photonHardN;   //!
   TBranch        *b_neutralMax;   //!
   TBranch        *b_neutralSum;   //!
   TBranch        *b_neutralN;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_eN;   //!
   TBranch        *b_muMax;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_muN;   //!
   TBranch        *b_matchedPt;   //!
   TBranch        *b_matchedR;   //!
   TBranch        *b_beamId1;   //!
   TBranch        *b_beamId2;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!
   TBranch        *b_gensubid;   //!
   TBranch        *b_smpt;   //!

};


void setupJetTree(TTree *t,Jets &tJets,bool doCheck = 0)
{
   // Set branch addresses and branch pointers
   t->SetBranchAddress("evt", &tJets.evt, &tJets.b_evt);
   t->SetBranchAddress("b", &tJets.b, &tJets.b_b);
   t->SetBranchAddress("nref", &tJets.nref, &tJets.b_nref);
   t->SetBranchAddress("rawpt", tJets.rawpt, &tJets.b_rawpt);
   t->SetBranchAddress("jtpt", tJets.jtpt, &tJets.b_jtpt);
   t->SetBranchAddress("jteta", tJets.jteta, &tJets.b_jteta);
   t->SetBranchAddress("jty", tJets.jty, &tJets.b_jty);
   t->SetBranchAddress("jtphi", tJets.jtphi, &tJets.b_jtphi);
   t->SetBranchAddress("jtpu", tJets.jtpu, &tJets.b_jtpu);
   t->SetBranchAddress("jtm", tJets.jtm, &tJets.b_jtm);
   t->SetBranchAddress("discr_fr01", tJets.discr_fr01, &tJets.b_discr_fr01);
   t->SetBranchAddress("trackMax", tJets.trackMax, &tJets.b_trackMax);
   t->SetBranchAddress("trackSum", tJets.trackSum, &tJets.b_trackSum);
   t->SetBranchAddress("trackN", tJets.trackN, &tJets.b_trackN);
   t->SetBranchAddress("trackHardSum", tJets.trackHardSum, &tJets.b_trackHardSum);
   t->SetBranchAddress("trackHardN", tJets.trackHardN, &tJets.b_trackHardN);
   t->SetBranchAddress("chargedMax", tJets.chargedMax, &tJets.b_chargedMax);
   t->SetBranchAddress("chargedSum", tJets.chargedSum, &tJets.b_chargedSum);
   t->SetBranchAddress("chargedN", tJets.chargedN, &tJets.b_chargedN);
   t->SetBranchAddress("chargedHardSum", tJets.chargedHardSum, &tJets.b_chargedHardSum);
   t->SetBranchAddress("chargedHardN", tJets.chargedHardN, &tJets.b_chargedHardN);
   t->SetBranchAddress("photonMax", tJets.photonMax, &tJets.b_photonMax);
   t->SetBranchAddress("photonSum", tJets.photonSum, &tJets.b_photonSum);
   t->SetBranchAddress("photonN", tJets.photonN, &tJets.b_photonN);
   t->SetBranchAddress("photonHardSum", tJets.photonHardSum, &tJets.b_photonHardSum);
   t->SetBranchAddress("photonHardN", tJets.photonHardN, &tJets.b_photonHardN);
   t->SetBranchAddress("neutralMax", tJets.neutralMax, &tJets.b_neutralMax);
   t->SetBranchAddress("neutralSum", tJets.neutralSum, &tJets.b_neutralSum);
   t->SetBranchAddress("neutralN", tJets.neutralN, &tJets.b_neutralN);
   t->SetBranchAddress("eMax", tJets.eMax, &tJets.b_eMax);
   t->SetBranchAddress("eSum", tJets.eSum, &tJets.b_eSum);
   t->SetBranchAddress("eN", tJets.eN, &tJets.b_eN);
   t->SetBranchAddress("muMax", tJets.muMax, &tJets.b_muMax);
   t->SetBranchAddress("muSum", tJets.muSum, &tJets.b_muSum);
   t->SetBranchAddress("muN", tJets.muN, &tJets.b_muN);
   t->SetBranchAddress("matchedPt", tJets.matchedPt, &tJets.b_matchedPt);
   t->SetBranchAddress("matchedR", tJets.matchedR, &tJets.b_matchedR);
   t->SetBranchAddress("beamId1", &tJets.beamId1, &tJets.b_beamId1);
   t->SetBranchAddress("beamId2", &tJets.beamId2, &tJets.b_beamId2);
   t->SetBranchAddress("pthat", &tJets.pthat, &tJets.b_pthat);
   t->SetBranchAddress("refpt", tJets.refpt, &tJets.b_refpt);
   t->SetBranchAddress("refeta", tJets.refeta, &tJets.b_refeta);
   t->SetBranchAddress("refy", tJets.refy, &tJets.b_refy);
   t->SetBranchAddress("refphi", tJets.refphi, &tJets.b_refphi);
   t->SetBranchAddress("refdphijt", tJets.refdphijt, &tJets.b_refdphijt);
   t->SetBranchAddress("refdrjt", tJets.refdrjt, &tJets.b_refdrjt);
   t->SetBranchAddress("refparton_pt", tJets.refparton_pt, &tJets.b_refparton_pt);
   t->SetBranchAddress("refparton_flavor", tJets.refparton_flavor, &tJets.b_refparton_flavor);
   t->SetBranchAddress("refparton_flavorForB", tJets.refparton_flavorForB, &tJets.b_refparton_flavorForB);
   t->SetBranchAddress("ngen", &tJets.ngen, &tJets.b_ngen);
   t->SetBranchAddress("genmatchindex", tJets.genmatchindex, &tJets.b_genmatchindex);
   t->SetBranchAddress("genpt", tJets.genpt, &tJets.b_genpt);
   t->SetBranchAddress("geneta", tJets.geneta, &tJets.b_geneta);
   t->SetBranchAddress("geny", tJets.geny, &tJets.b_geny);
   t->SetBranchAddress("genphi", tJets.genphi, &tJets.b_genphi);
   t->SetBranchAddress("gendphijt", tJets.gendphijt, &tJets.b_gendphijt);
   t->SetBranchAddress("gendrjt", tJets.gendrjt, &tJets.b_gendrjt);
   t->SetBranchAddress("gensubid", tJets.gensubid, &tJets.b_gensubid);
   t->SetBranchAddress("smpt", tJets.smpt, &tJets.b_smpt);
   if (doCheck) {
      if (t->GetMaximum("nref")>450) cout <<"FATAL ERROR: Arrary size of nref too small!!!  "<<t->GetMaximum("nref")<<endl;
      if (t->GetMaximum("ngen")>26) cout <<"FATAL ERROR: Arrary size of ngen too small!!!  "<<t->GetMaximum("ngen")<<endl;
   }
}

