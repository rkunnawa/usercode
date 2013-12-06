//
// writetree_pp.C 
//
// Raghav Kunnawalkam Elayavalli - 12/3/2013 
// Macro to read a hiforest (pp) and remove the duplicate events and create a simpler tree with the same structure
// as the hiforests.  
//
//


// Header files
#include "hiforest_good/hiForest.h"
#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <cstdlib>
#include <cmath>
#include <iostream>

// namespace
using namespace std;

//global variables
const float ketacut=2.0;
const float kptrecocut=30.;
const float cPP = 1;

//function definitions 
void LoadLib();
void ShutoffBranches(HiForest */*hi*/);
//void FindLeadSubLeadJets(Jets */*mJets*/, int */*ljet*/);
TStopwatch timer;

void LoadLib()
{
    gSystem->Load("hiforest_good/hiForest_h.so");
}

void ShutoffBranches(HiForest *hi)
{
    
    //! added by pawan
    //! Select only the branches you want to use for the analysis
    //! This increases the speed for running
    
    //! For Hlt
    hi->hltTree->SetBranchStatus("*",0,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet40_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet60_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet80_NoJetID_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_PAJet100_NoJetID_v1",1,0);
    
    //! for Skim Tree
    hi->skimTree->SetBranchStatus("*",0,0);
    hi->skimTree->SetBranchStatus("pcollisionEventSelection",1,0);
    hi->skimTree->SetBranchStatus("pPAcollisionEventSelectionPA",1,0);
    hi->skimTree->SetBranchStatus("pHBHENoiseFilter",1,0);
    
    //! Evt tree
    hi->evtTree->SetBranchStatus("*",0,0);
    hi->evtTree->SetBranchStatus("run",1,0);
    hi->evtTree->SetBranchStatus("evt",1,0);
    hi->evtTree->SetBranchStatus("vx",1,0);
    hi->evtTree->SetBranchStatus("vy",1,0);
    hi->evtTree->SetBranchStatus("vz",1,0);
    hi->evtTree->SetBranchStatus("hiNtracks",1,0);
    /*
    hi->ak2jetTree->SetBranchStatus("*",0,0);
    hi->ak2jetTree->SetBranchStatus("nref",1,0);
    hi->ak2jetTree->SetBranchStatus("rawpt",1,0);
    hi->ak2jetTree->SetBranchStatus("jtpt",1,0);
    hi->ak2jetTree->SetBranchStatus("jteta",1,0);
    hi->ak2jetTree->SetBranchStatus("jtphi",1,0);
    hi->ak2jetTree->SetBranchStatus("chargedMax",1,0);
    hi->ak2jetTree->SetBranchStatus("trackMax",1,0);
    hi->ak2jetTree->SetBranchStatus("chargedSum",1,0);
    hi->ak2jetTree->SetBranchStatus("photonSum",1,0);
    hi->ak2jetTree->SetBranchStatus("neutralSum",1,0);
    hi->ak2jetTree->SetBranchStatus("trackSum",1,0);
    hi->ak2jetTree->SetBranchStatus("photonMax",1,0);
    hi->ak2jetTree->SetBranchStatus("neutralMax",1,0);
    */
    hi->ak3jetTree->SetBranchStatus("*",0,0);
    hi->ak3jetTree->SetBranchStatus("nref",1,0);
    hi->ak3jetTree->SetBranchStatus("rawpt",1,0);
    hi->ak3jetTree->SetBranchStatus("jtpt",1,0);
    hi->ak3jetTree->SetBranchStatus("jteta",1,0);
    hi->ak3jetTree->SetBranchStatus("jtphi",1,0);
    hi->ak3jetTree->SetBranchStatus("chargedMax",1,0);
    hi->ak3jetTree->SetBranchStatus("trackMax",1,0);
    hi->ak3jetTree->SetBranchStatus("chargedSum",1,0);
    hi->ak3jetTree->SetBranchStatus("photonSum",1,0);
    hi->ak3jetTree->SetBranchStatus("neutralSum",1,0);
    hi->ak3jetTree->SetBranchStatus("trackSum",1,0);
    hi->ak3jetTree->SetBranchStatus("photonMax",1,0);
    hi->ak3jetTree->SetBranchStatus("neutralMax",1,0);
    
    hi->ak4jetTree->SetBranchStatus("*",0,0);
    hi->ak4jetTree->SetBranchStatus("nref",1,0);
    hi->ak4jetTree->SetBranchStatus("rawpt",1,0);
    hi->ak4jetTree->SetBranchStatus("jtpt",1,0);
    hi->ak4jetTree->SetBranchStatus("jteta",1,0);
    hi->ak4jetTree->SetBranchStatus("jtphi",1,0);
    hi->ak4jetTree->SetBranchStatus("chargedMax",1,0);
    hi->ak4jetTree->SetBranchStatus("trackMax",1,0);
    hi->ak4jetTree->SetBranchStatus("chargedSum",1,0);
    hi->ak4jetTree->SetBranchStatus("photonSum",1,0);
    hi->ak4jetTree->SetBranchStatus("neutralSum",1,0);
    hi->ak4jetTree->SetBranchStatus("trackSum",1,0);
    hi->ak4jetTree->SetBranchStatus("photonMax",1,0);
    hi->ak4jetTree->SetBranchStatus("neutralMax",1,0);
    
    
    //! Track tree
    //hi->trackTree->SetBranchStatus("*",0,0);
    //hi->trackTree->SetBranchStatus("nTrk",1,0);
    // hi->trackTree->SetBranchStatus("trkPt",1,0);
    // hi->trackTree->SetBranchStatus("trkEta",1,0);
    // hi->trackTree->SetBranchStatus("trkPhi",1,0);
    // hi->trackTree->SetBranchStatus("highPurity",1,0);
    // hi->trackTree->SetBranchStatus("trkDz1",1,0);
    // hi->trackTree->SetBranchStatus("trkDzError1",1,0);
    // hi->trackTree->SetBranchStatus("trkDxy1",1,0);
    // hi->trackTree->SetBranchStatus("trkDxyError1",1,0);
    
    /*
     //! PF jet tree
     hi->ak2PFJetTree->SetBranchStatus("*",0,0);
     hi->ak2PFJetTree->SetBranchStatus("nref",1,0);
     hi->ak2PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak2PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak2PFJetTree->SetBranchStatus("jteta",1,0);
     hi->ak2PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak2PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    /*
     hi->ak3PFJetTree->SetBranchStatus("*",0,0);
     hi->ak3PFJetTree->SetBranchStatus("nref",1,0);
     hi->ak3PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak3PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak3PFJetTree->SetBranchStatus("jteta",1,0);
     hi->ak3PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak3PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    /*
     hi->ak4PFJetTree->SetBranchStatus("*",0,0);
     hi->ak4PFJetTree->SetBranchStatus("nref",1,0);
     hi->ak4PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak4PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak4PFJetTree->SetBranchStatus("jteta",1,0);
     hi->ak4PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak4PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    /*
     hi->ak5PFJetTree->SetBranchStatus("*",0,0);
     hi->ak5PFJetTree->SetBranchStatus("nref",1,0);
     hi->ak5PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak5PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak5PFJetTree->SetBranchStatus("jteta",1,0);
     hi->ak5PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak5PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    
    //
    /*
     hi->akPu2PFJetTree->SetBranchStatus("*",0,0);
     hi->akPu2PFJetTree->SetBranchStatus("nref",1,0);
     hi->akPu2PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu2PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu2PFJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu2PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu2PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    /*
     hi->akPu3PFJetTree->SetBranchStatus("*",0,0);
     hi->akPu3PFJetTree->SetBranchStatus("nref",1,0);
     hi->akPu3PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu3PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu3PFJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu3PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu3PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    /*
     hi->akPu4PFJetTree->SetBranchStatus("*",0,0);
     hi->akPu4PFJetTree->SetBranchStatus("nref",1,0);
     hi->akPu4PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu4PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu4PFJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu4PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu4PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    /*
     hi->akPu5PFJetTree->SetBranchStatus("*",0,0);
     hi->akPu5PFJetTree->SetBranchStatus("nref",1,0);
     hi->akPu5PFJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu5PFJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu5PFJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu5PFJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu5PFJetTree->SetBranchStatus("trackMax",1,0);
     */
    //
    
    
    //! Calo jet trees
    /*
     hi->ak2CaloJetTree->SetBranchStatus("*",0,0);
     hi->ak2CaloJetTree->SetBranchStatus("nref",1,0);
     hi->ak2CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak2CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak2CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->ak2CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak2CaloJetTree->SetBranchStatus("trackMax",1,0);
     
     hi->ak3CaloJetTree->SetBranchStatus("*",0,0);
     hi->ak3CaloJetTree->SetBranchStatus("nref",1,0);
     hi->ak3CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak3CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak3CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->ak3CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak3CaloJetTree->SetBranchStatus("trackMax",1,0);
     
     hi->ak4CaloJetTree->SetBranchStatus("*",0,0);
     hi->ak4CaloJetTree->SetBranchStatus("nref",1,0);
     hi->ak4CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak4CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak4CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->ak4CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak4CaloJetTree->SetBranchStatus("trackMax",1,0);
     
     hi->ak5CaloJetTree->SetBranchStatus("*",0,0);
     hi->ak5CaloJetTree->SetBranchStatus("nref",1,0);
     hi->ak5CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->ak5CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->ak5CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->ak5CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->ak5CaloJetTree->SetBranchStatus("trackMax",1,0);
     */
    
    ////
    /*
     hi->akPu2CaloJetTree->SetBranchStatus("*",0,0);
     hi->akPu2CaloJetTree->SetBranchStatus("nref",1,0);
     hi->akPu2CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu2CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu2CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu2CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu2CaloJetTree->SetBranchStatus("trackMax",1,0);
     
     hi->akPu3CaloJetTree->SetBranchStatus("*",0,0);
     hi->akPu3CaloJetTree->SetBranchStatus("nref",1,0);
     hi->akPu3CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu3CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu3CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu3CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu3CaloJetTree->SetBranchStatus("trackMax",1,0);
     
     hi->akPu4CaloJetTree->SetBranchStatus("*",0,0);
     hi->akPu4CaloJetTree->SetBranchStatus("nref",1,0);
     hi->akPu4CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu4CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu4CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu4CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu4CaloJetTree->SetBranchStatus("trackMax",1,0);
     
     hi->akPu5CaloJetTree->SetBranchStatus("*",0,0);
     hi->akPu5CaloJetTree->SetBranchStatus("nref",1,0);
     hi->akPu5CaloJetTree->SetBranchStatus("rawpt",1,0);
     hi->akPu5CaloJetTree->SetBranchStatus("jtpt",1,0);
     hi->akPu5CaloJetTree->SetBranchStatus("jteta",1,0);
     hi->akPu5CaloJetTree->SetBranchStatus("jtphi",1,0);
     hi->akPu5CaloJetTree->SetBranchStatus("trackMax",1,0);
     */
    ///
    std::cout<<"Loaded all tree variables "<<std::endl;
    
}




#ifdef _MAKECINT_
#pragma link C++ class HiForest;
#pragma link C++ class Hlts;
#pragma link C++ class Skims;
#pragma link C++ class Evts;
#pragma link C++ class Jets;
#endif

int writetree_pp(char *ksp = "ppJet40"){

  timer.Start();
  LoadLib();

 
  TString inname="";
  if(strcmp(ksp,"ppJet40")==0)inname = "/hadoop/store/user/belt/hiForest2/PP2013_HiForest_PromptReco_JSon_Jet40Jet60_ppTrack_forestv84.root";
  else if(strcmp(ksp,"ppJet80")==0)inname = "/hadoop/store/user/belt/hiForest2/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root";
    
  //! Load Lib
  //gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest_h.so");
    
  //! Define the input file and HiForest
  //! CMSSW_5_3_3
  HiForest *c = new HiForest(inname,Form("Forest%s",ksp),cPP);
  cout<<"Loaded the hiforest tree : "<<c->GetName()<<endl;
  ShutoffBranches(c);
    

  TFile *fout = new TFile(Form("ntuple_2013_%s_v2.root",ksp),"RECREATE");
    
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"My hiForest Tree : " <<c->GetName()<<"\t Entries "<<c->GetEntries()<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
    
  //! shut off jet trees
  //c->hasAk2CaloJetTree=0;
  //c->hasAk4CaloJetTree=0;
  //c->hasAk3CaloJetTree=0;
  //c->hasAk5CaloJetTree=0;
    
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  //c->hasAkPu5CaloJetTree=0;
    
  //c->hasAk2PFJetTree=0;
  //c->hasAk4PFJetTree=0;
  //c->hasAk5PFJetTree=0;
    
  //c->hasAkPu2PFJetTree=0;
  //c->hasAkPu4PFJetTree=0;
  //c->hasAkPu5PFJetTree=0;
    
  c->hasTrackTree=0;
    
  //! For jets
  //Jets *mJets2 = 0;
  Jets *mJets3 = 0;
  Jets *mJets4 = 0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;

  //TTree *jetR2Tree = new TTree("jetR2","ak2PF");
  TTree *jetR3Tree = new TTree("jetR3","ak3PF");
  TTree *jetR4Tree = new TTree("jetR4","ak4PF");
  TTree *evtTree = new TTree("evt","evt");

  // declare the event variables.
  int evt;
  int run;
  float vx;
  float vy;
  float vz;
  int jet40;
  int jet60;
  int jet80;
  int jet100;
  int ntrk;

  // declare the jet variables
  /*
  int nrefe2;
  float pt2[1000];
  float raw2[1000];
  float eta2[1000];
  float phi2[1000];
  float chMax2[1000];
  float trkMax2[1000];
  float chSum2[1000];
  float phSum2[1000];
  float neSum2[1000];
  float trkSum2[1000];
  float phMax2[1000];
  float neMax2[1000];
  */

  int nrefe3;
  float pt3[1000];
  float raw3[1000];
  float eta3[1000];
  float phi3[1000];
  float chMax3[1000];
  float trkMax3[1000];
  float chSum3[1000];
  float phSum3[1000];
  float neSum3[1000];
  float trkSum3[1000];
  float phMax3[1000];
  float neMax3[1000];

  int nrefe4;
  float pt4[1000];
  float raw4[1000];
  float eta4[1000];
  float phi4[1000];
  float chMax4[1000];
  float trkMax4[1000];
  float chSum4[1000];
  float phSum4[1000];
  float neSum4[1000];
  float trkSum4[1000];
  float phMax4[1000];
  float neMax4[1000];
    
  //set the branches in the trees. 
  evtTree->Branch("evt",&evt,"evt/I");
  evtTree->Branch("run",&run,"run/I");
  evtTree->Branch("vx",&vx,"vx/F");
  evtTree->Branch("vy",&vy,"vy/F");
  evtTree->Branch("vz",&vz,"vz/F");
  evtTree->Branch("jet40",&jet40,"jet40/I");
  evtTree->Branch("jet60",&jet60,"jet60/I");
  evtTree->Branch("jet80",&jet80,"jet80/I");
  evtTree->Branch("jet100",&jet100,"jet100/I");
  evtTree->Branch("ntrk",&ntrk,"ntrk/I");
  /*
  jetR2Tree->Branch("nrefe",&nrefe2,"nrefe/I");
  jetR2Tree->Branch("pt",&pt2,"pt[nrefe]/F");
  jetR2Tree->Branch("raw",&raw2,"raw[nrefe]/F");
  jetR2Tree->Branch("eta",&eta2,"eta[nrefe]/F");
  jetR2Tree->Branch("phi",&phi2,"phi[nrefe]/F");
  jetR2Tree->Branch("chMax",&chMax2,"chMax[nrefe]/F");
  jetR2Tree->Branch("trkMax",&trkMax2,"trkMax[nrefe]/F");
  jetR2Tree->Branch("phMax",&phMax2,"phMax[nrefe]/F");
  jetR2Tree->Branch("neMax",&neMax2,"neMax[nrefe]/F");
  jetR2Tree->Branch("chSum",&chSum2,"chSum[nrefe]/F");
  jetR2Tree->Branch("phSum",&phSum2,"phSum[nrefe]/F");
  jetR2Tree->Branch("neSum",&neSum2,"neSum[nrefe]/F");
  jetR2Tree->Branch("trkSum",&trkSum2,"trkSum[nrefe]/F");
  */
  jetR3Tree->Branch("nrefe",&nrefe3,"nrefe/I");
  jetR3Tree->Branch("pt",&pt3,"pt[nrefe]/F");
  jetR3Tree->Branch("raw",&raw3,"raw[nrefe]/F");
  jetR3Tree->Branch("eta",&eta3,"eta[nrefe]/F");
  jetR3Tree->Branch("phi",&phi3,"phi[nrefe]/F");
  jetR3Tree->Branch("chMax",&chMax3,"chMax[nrefe]/F");
  jetR3Tree->Branch("trkMax",&trkMax3,"trkMax[nrefe]/F");
  jetR3Tree->Branch("phMax",&phMax3,"phMax[nrefe]/F");
  jetR3Tree->Branch("neMax",&neMax3,"neMax[nrefe]/F");
  jetR3Tree->Branch("chSum",&chSum3,"chSum[nrefe]/F");
  jetR3Tree->Branch("phSum",&phSum3,"phSum[nrefe]/F");
  jetR3Tree->Branch("neSum",&neSum3,"neSum[nrefe]/F");
  jetR3Tree->Branch("trkSum",&trkSum3,"trkSum[nrefe]/F");

  jetR4Tree->Branch("nrefe",&nrefe4,"nrefe/I");
  jetR4Tree->Branch("pt",&pt4,"pt[nrefe]/F");
  jetR4Tree->Branch("raw",&raw4,"raw[nrefe]/F");
  jetR4Tree->Branch("eta",&eta4,"eta[nrefe]/F");
  jetR4Tree->Branch("phi",&phi4,"phi[nrefe]/F");
  jetR4Tree->Branch("chMax",&chMax4,"chMax[nrefe]/F");
  jetR4Tree->Branch("trkMax",&trkMax4,"trkMax[nrefe]/F");
  jetR4Tree->Branch("phMax",&phMax4,"phMax[nrefe]/F");
  jetR4Tree->Branch("neMax",&neMax4,"neMax[nrefe]/F");
  jetR4Tree->Branch("chSum",&chSum4,"chSum[nrefe]/F");
  jetR4Tree->Branch("phSum",&phSum4,"phSum[nrefe]/F");
  jetR4Tree->Branch("neSum",&neSum4,"neSum[nrefe]/F");
  jetR4Tree->Branch("trkSum",&trkSum4,"trkSum[nrefe]/F");

     
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
  //for (Long64_t ievt=0; ievt<100;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);

    bool evSel = false;
    evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pPAcollisionEventSelectionPA && (c->hlt.HLT_PAJet40_NoJetID_v1 || c->hlt.HLT_PAJet60_NoJetID_v1 || c->hlt.HLT_PAJet80_NoJetID_v1 || c->hlt.HLT_PAJet100_NoJetID_v1);
	
    if(!evSel)continue;
        
   
    run = c->evt.run;
    evt = c->evt.evt;
    vx = c->evt.vx;
    vy = c->evt.vy;
    vz = c->evt.vz;
    jet40  = c->hlt.HLT_PAJet40_NoJetID_v1;
    jet60  = c->hlt.HLT_PAJet60_NoJetID_v1;
    jet80  = c->hlt.HLT_PAJet80_NoJetID_v1;
    jet100 = c->hlt.HLT_PAJet100_NoJetID_v1;
    ntrk = c->evt.hiNtracks;

    if(ievt%10000 == 0) cout<<" ******** Event # "<< ievt <<"\t Run " <<run<<endl;
    /*
    mJets2 = &(c->ak2PF);
    nrefe2 = mJets2->nref;

    for (int i = 0; i<nrefe2; i++) {
      pt2[i]     = mJets2->jtpt[i];
      eta2[i]    = mJets2->jteta[i];
      phi2[i]    = mJets2->jtphi[i];
      raw2[i]    = mJets2->rawpt[i];
      chMax2[i]  = mJets2->chargedMax[i];
      trkMax2[i]  = mJets2->trackMax[i];
      chSum2[i]  = mJets2->chargedSum[i];
      phSum2[i]  = mJets2->photonSum[i];
      neSum2[i]  = mJets2->neutralSum[i];
      trkSum2[i] = mJets2->trackSum[i];
      phSum2[i]  = mJets2->photonMax[i];
      neMax2[i]  = mJets2->neutralMax[i];
    }
    */

    mJets3 = &(c->ak3PF);
    nrefe3 = mJets3->nref;

    for (int i = 0; i<nrefe3; i++) {
      pt3[i]     = mJets3->jtpt[i];
      eta3[i]    = mJets3->jteta[i];
      phi3[i]    = mJets3->jtphi[i];
      raw3[i]    = mJets3->rawpt[i];
      chMax3[i]  = mJets3->chargedMax[i];
      trkMax3[i]  = mJets3->trackMax[i];
      chSum3[i]  = mJets3->chargedSum[i];
      phSum3[i]  = mJets3->photonSum[i];
      neSum3[i]  = mJets3->neutralSum[i];
      trkSum3[i] = mJets3->trackSum[i];
      phSum3[i]  = mJets3->photonMax[i];
      neMax3[i]  = mJets3->neutralMax[i];
    }

    mJets4 = &(c->ak4PF);
    nrefe4 = mJets4->nref;

    for (int i = 0; i<nrefe4; i++) {
      pt4[i]     = mJets4->jtpt[i];
      eta4[i]    = mJets4->jteta[i];
      phi4[i]    = mJets4->jtphi[i];
      raw4[i]    = mJets4->rawpt[i];
      chMax4[i]  = mJets4->chargedMax[i];
      trkMax4[i]  = mJets4->trackMax[i];
      chSum4[i]  = mJets4->chargedSum[i];
      phSum4[i]  = mJets4->photonSum[i];
      neSum4[i]  = mJets4->neutralSum[i];
      trkSum4[i] = mJets4->trackSum[i];
      phSum4[i]  = mJets4->photonMax[i];
      neMax4[i]  = mJets4->neutralMax[i];
    }

    evtTree->Fill();
    //jetR2Tree->Fill();
    jetR3Tree->Fill();
    jetR4Tree->Fill();

  }  

  fout->cd();
  fout->Write();
  fout->Close();
    
  // time out:
  timer.Stop();
  float rtime  = timer.RealTime();
  float ctime  = timer.CpuTime();
    
  cout<<"\t"<<endl;
  cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<endl;
  cout<<"\t"<<endl;
  cout<<"Good bye : " <<"\t"<<endl;
    
  return 0;

  
  




}

