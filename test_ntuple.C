#include <iostream>

#include <stdio.h>

#include <TRandom.h>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"


void test_ntuple(){

  TH1::SetDefaultSumw2();
  
  TFile *fpbpb1 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet80.root");
  TFile *fpbpb2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet80_v2.root");

  TTree *jet1 = (TTree*)fpbpb1->Get("ntjet");
  TTree *jet2 = (TTree*)fpbpb2->Get("ntjet");
  TTree *evt2 = (TTree*)fpbpb2->Get("ntevt");

  TH1F *hHiBin_2 = new TH1F("hHiBin_2","HiBin histogram from the v2 ntuple",40,0,40);
  TH1F *hHiBin_1 = new TH1F("hHiBin_1","HiBin histogram from the v1 ntuple",40,0,40);
  
  TH1F *hp_T_2 = new TH1F("hp_T_2","p_T distribution from v2 ntuple",50,0,500);
  TH1F *hp_T_1 = new TH1F("hp_T_1","p_T distribution from v1 ntuple",50,0,500);

  jet2->SetMakeClass(1);
  evt2->SetMakeClass(1);
  float hiBin2 = 0;
  float pt2 = 0;
  float jet80_2 = 0;
  float nrefe2 = 0;
  evt2->SetBranchAddress("bin",&hiBin2);
  jet2->SetBranchAddress("pt",&pt2);
  evt2->SetBranchAddress("jet80",&jet80_2);
  evt2->SetBranchAddress("nrefe",&nrefe2);

  jet1->SetMakeClass(1);
  float hiBin1 = 0;
  float pt1 = 0;
  float jet80_1 = 0;
  float nrefe1 = 0;
  jet1->SetBranchAddress("bin",&hiBin1);
  jet1->SetBranchAddress("pt",&pt1);
  jet1->SetBranchAddress("jet80",&jet80_1);
  jet1->SetBranchAddress("nrefe",&nrefe1);

  Long64_t njet1 = jet1->GetEntries();
  Long64_t njet2 = jet2->GetEntries();
  Long64_t nevt2 = evt2->GetEntries();
  
  cout<<"number of entries in jet tree ntuple v1 = "<<njet1<<endl;
  cout<<"number of entries in jet tree ntuple v2 = "<<njet2<<endl;
  cout<<"number of entries in evt tree ntuple v2 = "<<nevt2<<endl;

  cout<<endl<<endl<<endl<<"starting loop 1 - for the event variables"<<endl<<endl<<endl;
  //read the histograms from ntuple type 1;
  for(Long64_t i = 0;i<njet1;){ 
    // loop for jet1 starts - note only one loop for this ntuple 
    // which means that we can go either one of two ways, do i++ in the beginning and get only the jet variables like 
    // pt, eta, etc... or we can do i = i+ nrefe1 and get the event variables like vz, bin etc...  
    // there might be another way to get both of them in the same loop but i cant think of that now. 
    jet1->GetEntry(i);
    //if(i%10000 == 0)cout<<"entry number = "<<i<<endl;
    hHiBin_1->Fill(hiBin1);
    i = i+ nrefe1;//we increment i here to get to the next event 
  }

  //check the histograms
  hHiBin_1->Print("base");


  for(Long64_t i = 0;i<njet1;i++){ //again this is going for the jet variables. so only one loop
    //if(i%10000 == 0)cout<<"entry number = "<<i<<endl;
    jet1->GetEntry(i);
    if(jet80_1 == 1)
      hp_T_1->Fill(pt1);
  }
  hp_T_1->Print("base");

  // now entering the ntuple 2 loop to get the histograms. this is going to much simpler than the first one due to the 
  // structure of the loops. 

  Long64_t counter = 0; //counter for our jet loop inside here. 

  for(Long64_t i = 0;i<nevt2;i++){
    // get the information from the event loop here 
    //if(i%10000 == 0)cout<<"entry number = "<<i<<endl;
    evt2->GetEntry(i);
    hHiBin_2->Fill(hiBin2);
    //cout<<"counter = "<<counter<<endl;
    for(int j = counter;j<counter+nrefe2;j++){
      jet2->GetEntry(counter);
      if(jet80_2 == 1)
	hp_T_2->Fill(pt2);
    }
    counter = counter + nrefe2;
  }
  hHiBin_2->Print("base");
  hp_T_2->Print("base");

  // now lets draw the histograms and compare the results
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,2);
  c1->cd(1);
  hHiBin_1->Draw();
  c1->cd(2);
  hHiBin_2->Draw();
  c1->cd(3);
  c1->cd(3)->SetLogy();
  hp_T_1->Draw();
  c1->cd(4);
  c1->cd(4)->SetLogy();
  hp_T_2->Draw();

  c1->SaveAs("test_ntuple.gif","RECREATE");
  

}
