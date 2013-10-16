
//
//  writentuple_pbpb.C
//  
//  Pawan - github.com/pawannetrakanti/pawan/blob/master/Data/writentuple_pbpb.C
//  Raghav Kunnawalkam Elayavalli on 9/18/13.
//  
//

// From the prompt reco
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
using namespace std;

const float ketacut=2.0;
const float kptrecocut=30.;
const float cPP = 0;

void LoadLib();
void ShutoffBranches(HiForest */*hi*/);
//void FindLeadSubLeadJets(Jets */*mJets*/, int */*ljet*/);



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




TStopwatch timer;
int writentuple_pbpb(char *ksp="pbpbJet55")
{
    
    timer.Start();
    
    LoadLib();
    
    TString inname="";
    if(strcmp(ksp,"pbpbJet80")==0) inname = "/hadoop/store/user/belt/HiForest2_v21_UMD/HiForest-promptskim-hiForest2_v21_HLTFix.root";//filename
    else if(strcmp(ksp,"pbpbJet65")==0) inname = "/hadoop/store/user/belt/hiForest2/HiForest_JetRAA_HIJet65_v24.root";
    else if(strcmp(ksp,"pbpbJet55")==0) inname = "/hadoop/store/user/belt/hiForest2/HiForest_JetRAA_HIJet55_v25.root";
    
    //! Load Lib
    //gSystem->Load("/afs/cern.ch/user/p/pawan/scratch0/CMSSW_6_2_0/src/work/pPb/HiForest/V3/hiForest_h.so");
    
    //! Define the input file and HiForest
    //! CMSSW_5_3_3
    HiForest *c = new HiForest(inname,Form("Forest%s",ksp),cPP);
    cout<<"Loaded the hiforest tree : "<<c->GetName()<<endl;
    ShutoffBranches(c);
    
    DuplicateEvents dupEvt(inname);
    dupEvt.MakeList();
    
    
    
    //! Output file
    //! HIHighPt
    TFile *fout = new TFile(Form("ntuple_2011_%s.root",ksp),"RECREATE");
    
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
    
    string jetVars = "";
    jetVars += "evt:run:bin:vz:trig:jet55:jet65:jet80:ntrk:nrefe:pt:raw:eta:phi:chMax:trkMax:chSum:phSum:neSum";
    TNtuple *ntjet=0;
    ntjet = new TNtuple("ntjet","",jetVars.data());
    
    for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
        //for (Long64_t ievt=0; ievt<100;ievt++) {//! event loop
        //! load the hiForest event
        c->GetEntry(ievt);
        
        if(dupEvt.occurence[ievt] == 2)continue;
        
        //! events with Single vertex
        bool evSel = false;
        float trig=-9;
        // if(strcmp(ksp,"pbpbJet55")==0){
        //   evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pcollisionEventSelection && c->hlt.HLT_HIJet55_v1;
        //   trig=1;
        // }
        // else if(strcmp(ksp,"pbpbJet65")==0){
        //   evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter && c->skim.pcollisionEventSelection && c->hlt.HLT_HIJet65_v1;
        //   trig=2;
        // }
        // else if(strcmp(ksp,"pbpbJet80")==0){
        //   evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter && c->skim.pcollisionEventSelection && c->hlt.HLT_HIJet80_v1;
        //   trig=3;
        // }
        
        
        evSel = fabs(c->evt.vz)<15. && c->skim.pHBHENoiseFilter  && c->skim.pcollisionEventSelection && (c->hlt.HLT_HIJet55_v1 || c->hlt.HLT_HIJet65_v1 || c->hlt.HLT_HIJet80_v1);
        if(!evSel)continue;
        
        float pt = -9,raw    = -9,eta    = -9,phi    = -9,chMax  = -9,trkMax = -9,chSum  = -9,phSum  = -9,neSum  = -9,nrefe = 0;
        
        
        float run   = c->evt.run;
        float evt   = c->evt.evt;
        float vz    = c->evt.vz;
        
        float jet55  = c->hlt.HLT_HIJet55_v1;
        float jet65  = c->hlt.HLT_HIJet65_v1;
        float jet80  = c->hlt.HLT_HIJet80_v1;
        
        float ntrk    = c->evt.hiNtracks;
        
        float bin = c->evt.hiBin;
        
        if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<run<<std::endl;
        //std::cout<<" ********** Event # " <<ievt<<"\t Run : "<<run<<std::endl;
        mJets = &(c->akPu3PF);
	nrefe = mJets->nref;
        
        //int *ljet = new int[3];
        //FindLeadSubLeadJets(mJets,ljet);
        
        //int jtLead = -1, jtSubLead = -1, jtThird = -1;
        
        //if(ljet[0] >=0 ) jtLead    = ljet[0];
        //if(ljet[1] >=0 ) jtSubLead = ljet[1];
        //if(ljet[2] >=0 ) jtThird   = ljet[2];
        
        //if(jtLead<0)continue;
        /*
        if(jtLead > -1){
            pt1     = mJets->jtpt[jtLead];
            eta1    = mJets->jteta[jtLead];
            phi1    = mJets->jtphi[jtLead];
            raw1    = mJets->rawpt[jtLead];
            //chMax1  = mJets->chargedMax[jtLead];
            chMax1  = mJets->trackMax[jtLead];
            chSum1  = mJets->chargedSum[jtLead];
            phSum1  = mJets->photonSum[jtLead];
            neSum1  = mJets->neutralSum[jtLead];
        }
        
        if(jtSubLead > -1){
            pt2     = mJets->jtpt[jtSubLead];
            eta2    = mJets->jteta[jtSubLead];
            phi2    = mJets->jtphi[jtSubLead];
            raw2    = mJets->rawpt[jtSubLead];
            //chMax2  = mJets->chargedMax[jtSubLead];
            chMax2  = mJets->trackMax[jtSubLead];
            chSum2  = mJets->chargedSum[jtSubLead];
            phSum2  = mJets->photonSum [jtSubLead];
            neSum2  = mJets->neutralSum[jtSubLead];
        }
        
        if(jtThird > -1){
            pt3     = mJets->jtpt[jtThird];
            eta3    = mJets->jteta[jtThird];
            phi3    = mJets->jtphi[jtThird];
            raw3    = mJets->rawpt[jtThird];
            //chMax3  = mJets->chargedMax[jtThird];
            chMax3  = mJets->trackMax[jtThird];
            chSum3  = mJets->chargedSum[jtThird];
            phSum3  = mJets->photonSum [jtThird];
            neSum3  = mJets->neutralSum[jtThird];
        }
        */
        for (int i = 0; i<mJets->nref; i++) {
            pt     = mJets->jtpt[i];
            eta    = mJets->jteta[i];
            phi    = mJets->jtphi[i];
            raw    = mJets->rawpt[i];
            chMax  = mJets->chargedMax[i];
            trkMax  = mJets->trackMax[i];
            chSum  = mJets->chargedSum[i];
            phSum  = mJets->photonSum[i];
            neSum  = mJets->neutralSum[i];
            
            float jentry[] = {evt,run,vz,bin,trig,jet55,jet65,jet80,ntrk,nrefe,
			      pt,raw,eta,phi,chMax,trkMax,chSum,phSum,neSum
            };
            
            ntjet->Fill(jentry);
        }

        //delete [] ljet;
    }//! event loop ends
    
    
    //! Write to output file
    fout->cd();
    fout->Write();
    fout->Close();
    
    //! Check
    timer.Stop();
    float rtime  = timer.RealTime();
    float ctime  = timer.CpuTime();
    
    std::cout<<"\t"<<std::endl;
    std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
    std::cout<<"\t"<<std::endl;
    std::cout<<"Good bye : " <<"\t"<<std::endl;
    
    return 0;
}
/*
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
    ljet[0]=-1; ljet[1]=-2; ljet[2]=-3;
    
    float tempt=-9;
    //! Get the leading jet
    for(int ij=0; ij<jetc->nref; ij++){
        if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]<kptrecocut || jetc->rawpt[ij]<15)continue;
        float jetpt = jetc->jtpt[ij];
        if(jetpt > tempt){
            tempt = jetpt;
            ljet[0] = ij;
        }
    }
    if(ljet[0]>=0){
        // Subleading
        tempt=-9;
        for(int ij=0; ij<jetc->nref; ij++){
            if(ij==ljet[0])continue;
            if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
            float jetpt = jetc->jtpt[ij];
            if (jetpt > tempt){
                tempt = jetpt;
                ljet[1] = ij;
            }
        }
        
        if(ljet[1]>=0){
            // third jet
            tempt=-9;
            for(int ij=0; ij<jetc->nref; ij++){
                if(ij==ljet[0] || ij==ljet[1])continue;
                if(fabs(jetc->jteta[ij])>ketacut || jetc->jtpt[ij]>jetc->jtpt[ljet[0]] || jetc->rawpt[ij]<15)continue;
                float jetpt = jetc->jtpt[ij];
                if (jetpt > tempt){
                    tempt = jetpt;
                    ljet[2] = ij;
                }
            }
        }//! third jet
    }
}*/
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
