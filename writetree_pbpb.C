//
// writetree_pbpb.C 
//
// Raghav Kunnawalkam Elayavalli - 12/3/2013 
// Macro to read a hiforest (pbpb) and remove the duplicate events and create a simpler tree with the same structure
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
const float cPP = 0;

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
    //! PbPb
    hi->hltTree->SetBranchStatus("HLT_HIJet55_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_HIJet65_v1",1,0);
    hi->hltTree->SetBranchStatus("HLT_HIJet80_v1",1,0);
    
    //! pp 2013
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
    hi->evtTree->SetBranchStatus("hiBin",1,0);
    hi->evtTree->SetBranchStatus("hiNtracks",1,0);
    
    hi->akPu3jetTree->SetBranchStatus("*",0,0);
    hi->akPu3jetTree->SetBranchStatus("nref",1,0);
    hi->akPu3jetTree->SetBranchStatus("rawpt",1,0);
    hi->akPu3jetTree->SetBranchStatus("jtpt",1,0);
    hi->akPu3jetTree->SetBranchStatus("jteta",1,0);
    hi->akPu3jetTree->SetBranchStatus("jtphi",1,0);
    hi->akPu3jetTree->SetBranchStatus("trackMax",1,0);
    hi->akPu3jetTree->SetBranchStatus("chargedMax",1,0);
    hi->akPu3jetTree->SetBranchStatus("chargedSum",1,0);
    hi->akPu3jetTree->SetBranchStatus("photonSum",1,0);
    hi->akPu3jetTree->SetBranchStatus("neutralSum",1,0);
    hi->akPu3jetTree->SetBranchStatus("trackSum",1,0);
    hi->akPu3jetTree->SetBranchStatus("photonMax",1,0);
    hi->akPu3jetTree->SetBranchStatus("neutralMax",1,0);
    
    
    
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

class DuplicateEvents {
public:
    DuplicateEvents(TString infname){
        inf = TFile::Open(infname);
        t = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
    };
    ~DuplicateEvents(){
        delete inf;
    }
    
    void MakeList(){
        cout<<"starting make list to check for duplicate events"<<endl;
        evts.clear();
        occurence.clear();
        int run,evt;
        t->SetBranchAddress("run",&run);
        t->SetBranchAddress("evt",&evt);
        for(int i = 0;i<t->GetEntries();i++){
            t->GetEntry(i);
            if(i%100000==0) cout<<i<<" / "<<t->GetEntries()<<" run: "<<run<<" evt: "<<evt<<endl;
            int occur = (int)FindOccurences(run,evt);
            if(occur==0) occurence.push_back(1);
            else occurence.push_back(2);
            evts.push_back(std::make_pair(run,evt));
        }
    }
    int FindOccurences(int run, int evt){
        int noccur = count(evts.begin(),evts.end(),std::make_pair(run,evt));
        return noccur;
    }
    TFile* inf;
    TTree* t;
    vector <pair<int,int> > evts;
    vector <int> occurence;
};


int writetree_pbpb(char *ksp="pbpbJet55"){



  timer.Start();
    
  LoadLib();
    
  TString inname="";
  if(strcmp(ksp,"pbpbJet80")==0) inname = "/hadoop/store/user/belt/HiForest2_v21_UMD/HiForest-promptskim-hiForest2_v21_HLTFix.root";//filename
  else if(strcmp(ksp,"pbpbJet65")==0) inname = "/hadoop/store/user/belt/hiForest2/HiForest_JetRAA_HIJet65_v24.root";
  else if(strcmp(ksp,"pbpbJet55")==0) inname = "/hadoop/store/user/belt/hiForest2/HiForest_JetRAA_HIJet55_v25.root";
  
  // Load Lib
  //gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest_h.so");
  
  // Define the input file and HiForest
  // CMSSW_5_3_3
  HiForest *c = new HiForest(inname,Form("Forest%s",ksp),cPP);
  cout<<"Loaded the hiforest tree : "<<c->GetName()<<endl;
  ShutoffBranches(c);
  
  DuplicateEvents dupEvt(inname);
  dupEvt.MakeList();
  
    
  // Output file
  TFile *fout = new TFile(Form("ntuple_2011_%s_v2.root",ksp),"RECREATE");
  TTree *jetTree = new TTree("jet","jet");
  TTree *evtTree = new TTree("evt","evt");

  // declare the event variables.
  int evt;
  int run;
  float vx;
  float vy;
  float vz;
  int jet55;
  int jet65;
  int jet80;
  int nrefe;
  int ntrk;
  float bin;

  // declare the jet variables
  float pt[1000];
  float raw[1000];
  float eta[1000];
  float phi[1000];
  float chMax[1000];
  float trkMax[1000];
  float chSum[1000];
  float phSum[1000];
  float neSum[1000];
  float trkSum[1000];
  float phMax[1000];
  float neMax[1000];
    
  //set the branches in the trees. 
  evtTree->Branch("evt",&evt,"evt/I");
  evtTree->Branch("run",&run,"run/I");
  evtTree->Branch("bin",&bin,"bin/F");
  evtTree->Branch("vx",&vx,"vx/F");
  evtTree->Branch("vy",&vy,"vy/F");
  evtTree->Branch("vz",&vz,"vz/F");
  evtTree->Branch("jet55",&jet55,"jet55/I");
  evtTree->Branch("jet65",&jet65,"jet65/I");
  evtTree->Branch("jet80",&jet80,"jet80/I");
  evtTree->Branch("nrefe",&nrefe,"nrefe/I");
  evtTree->Branch("ntrk",&ntrk,"ntrk/I");

  jetTree->Branch("nrefe",&nrefe,"nrefe/I");
  jetTree->Branch("pt",&pt,"pt[nrefe]/F");
  jetTree->Branch("raw",&raw,"raw[nrefe]/F");
  jetTree->Branch("eta",&eta,"eta[nrefe]/F");
  jetTree->Branch("phi",&phi,"phi[nrefe]/F");
  jetTree->Branch("chMax",&chMax,"chMax[nrefe]/F");
  jetTree->Branch("trkMax",&trkMax,"trkMax[nrefe]/F");
  jetTree->Branch("phMax",&phMax,"phMax[nrefe]/F");
  jetTree->Branch("neMax",&neMax,"neMax[nrefe]/F");
  jetTree->Branch("chSum",&chSum,"chSum[nrefe]/F");
  jetTree->Branch("phSum",&phSum,"phSum[nrefe]/F");
  jetTree->Branch("neSum",&neSum,"neSum[nrefe]/F");
  jetTree->Branch("trkSum",&trkSum,"trkSum[nrefe]/F");

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
  // c->hasAk4CaloJetTree=0;
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
  Jets *mJets=0;
  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  
  for(Long64_t ievt = 0;ievt<nentries;ievt++){
  //for(Long64_t ievt =0;ievt<100;ievt++){
    
    //load the hiforest event 
    c->GetEntry(ievt);

    if(dupEvt.occurence[ievt]==2)continue;

    bool evSel = false;

    evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pcollisionEventSelection && (c->hlt.HLT_HIJet55_v1 || c->hlt.HLT_HIJet65_v1 || c->hlt.HLT_HIJet80_v1);
    if(!evSel)continue;
    
    run = c->evt.run;
    evt = c->evt.evt;
    vx = c->evt.vx;
    vy = c->evt.vy;
    vz = c->evt.vz;
    jet55 = c->hlt.HLT_HIJet55_v1;
    jet65 = c->hlt.HLT_HIJet65_v1;
    jet80 = c->hlt.HLT_HIJet80_v1;
    ntrk = c->evt.hiNtracks;
    bin = c->evt.hiBin;

    if(ievt%10000 == 0) cout<<" ******** Event # "<< ievt <<"\t Run " <<run<<endl;
    mJets = &(c->akPu3PF);
    nrefe = mJets->nref;

    for (int i = 0; i<mJets->nref; i++) {
      pt[i]     = mJets->jtpt[i];
      eta[i]    = mJets->jteta[i];
      phi[i]    = mJets->jtphi[i];
      raw[i]    = mJets->rawpt[i];
      chMax[i]  = mJets->chargedMax[i];
      trkMax[i]  = mJets->trackMax[i];
      chSum[i]  = mJets->chargedSum[i];
      phSum[i]  = mJets->photonSum[i];
      neSum[i]  = mJets->neutralSum[i];
      trkSum[i] = mJets->trackSum[i];
      phSum[i]  = mJets->photonMax[i];
      neMax[i]  = mJets->neutralMax[i];

    }
    evtTree->Fill();
    jetTree->Fill();

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


