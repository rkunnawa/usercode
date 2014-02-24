// macro to merge different HLT for ppb 2013 


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

TH1F *functionHist(TF1 *f, TH1F* h,char *fHistname)
{
	TH1F *hF = (TH1F*)h->Clone(fHistname);
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
		hF->SetBinContent(i,var);
		hF->SetBinError(i,0);
	}
	return hF;
}


// divide by bin width
void divideBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val/=h->GetBinWidth(i);
		valErr/=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}

void multiplyBinWidth(TH1 *h)
{
	h->Sumw2();
	for (int i=0;i<=h->GetNbinsX();i++)
	{
		Float_t val = h->GetBinContent(i);
		Float_t valErr = h->GetBinError(i);
		val*=h->GetBinWidth(i);
		valErr*=h->GetBinWidth(i);
		h->SetBinContent(i,val);
		h->SetBinError(i,valErr);
	}
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
}



void drawText(const char *text, float xp, float yp, int size){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(63);
	tex->SetTextSize(size);
	tex->SetTextColor(kBlack);
	tex->SetLineWidth(1);
	//tex->SetTextFont(42);
	tex->SetNDC();
	tex->Draw();
}


void putCMSPrel(double x, double y, double size){
	TLatex *tex=0;
	tex = new TLatex(x,y,"CMS Preliminary");
	tex->SetTextSize(size);
	tex->SetLineWidth(2);
	tex->SetNDC();
	tex->Draw();
}


TLegend *myLegend(double x1,double y1,double x2, double y2)
{
	TLegend *leg = new TLegend(x1,y1,x2,y2);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	return leg; 
	
}

// Remove bins with error > central value
void cleanup(TH1F *h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		double val1 = h->GetBinContent(i);
		double valErr1 = h->GetBinError(i);
		if (valErr1>=val1) {
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}   
	
}

static const int nbins_rec = 100;
static const double boundaries_rec[nbins_rec+1] = {
        0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170,
	180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310,
	320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450,
        460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590,
        600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730,
        740, 750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870,
        880, 890, 900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000
	};

Double_t boundaries_yaxian[] = {0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
Int_t nbins_yaxian = 21;

//Double_t boundaries_yaxian[] = {3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 507, 592, 1000};
//Int_t nbins_yaxian = 31;

// rebin the spectra
TH1F *rebin(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),1000,0,1000);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      double val=h->GetBinContent(i);
      double valErr=h->GetBinError(i);
      int binNum = hRebin->FindBin(h->GetBinCenter(i));
      double val1 = hRebin->GetBinContent(binNum);
      double valErr1 = hRebin->GetBinError(binNum);
      hRebin->SetBinContent(binNum,val+val1);
      hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
    }
  cleanup(hRebin);
  hRebin->SetName(histName);
  return hRebin;
}


// rebin the spectra
TH1F *rebin_yaxian(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_yaxian,boundaries_yaxian);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      double val=h->GetBinContent(i);
      double valErr=h->GetBinError(i);
      int binNum = hRebin->FindBin(h->GetBinCenter(i));
      double val1 = hRebin->GetBinContent(binNum);
      double valErr1 = hRebin->GetBinError(binNum);
      hRebin->SetBinContent(binNum,val+val1);
      hRebin->SetBinError(binNum,sqrt(valErr1*valErr1+valErr*valErr));
    }
  cleanup(hRebin);
  hRebin->SetName(histName);
  return hRebin;
}

void merge_ppb_HLT(){


  //for(int i = 0; i<=1000; i++){
  //  cout<<i<<", ";
  //}

  TH1::SetDefaultSumw2();
  
  TFile *fppb0_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppb_ppbJetMB_v2.root");
  //TFile *fppb0_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppb_pbp_ppbJetMB_v2.root");
  TFile *fppb1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppbJet40_v2.root");
  TFile *fppb2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppb_ppbJet80_v2.root");
  //TFile *fppb2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppb_pbp_ppbJet80_v2.root");

  

  //Number convention Files  - 0 - MB, 1- 40_60, 2- 80_100 

  TTree *jetppb0_v2 = (TTree*)fppb0_v2->Get("jetR3");
  TTree *jetppb1_v2 = (TTree*)fppb1_v2->Get("jetR3");
  TTree *jetppb2_v2 = (TTree*)fppb2_v2->Get("jetR3");

  TTree *evtppb0_v2 = (TTree*)fppb0_v2->Get("evt");
  TTree *evtppb1_v2 = (TTree*)fppb1_v2->Get("evt");
  TTree *evtppb2_v2 = (TTree*)fppb2_v2->Get("evt");

  jetppb0_v2->AddFriend(evtppb0_v2);
  jetppb1_v2->AddFriend(evtppb1_v2);
  jetppb2_v2->AddFriend(evtppb2_v2);


  //If we are running on HiForests then include these trees:
  /*
    TTree* jetppb0_v2 = (TTree*)fppb0_v2->Get("akPu3PFJetAnalyzer/t");
    TTree* evtppb0_v2 = (TTree*)fppb0_v2->Get("hiEvtAnalyzer/HiTree");
    TTree* hltppb0_v2 = (TTree*)fppb0_v2->Get("hltanalysis/HltTree");
    TTree* skimppb0_v2 = (TTree*)fppb0_v2->Get("skimanalysis/HltTree");
    
    jetppb0_v2->AddFriend(evtppb0_v2);
    jetppb0_v2->AddFriend(hltppb0_v2);
    jetppb0_v2->AddFriend(skimppb0_v2);

  */

  //number convention is going to change here to better understand the merging. 
  // 0 - MB
  // 1 - 20
  // 2 - 40
  // 3 - 60
  // 4 - 80
  // 5 - 100
  // 6 - 120

  

  //these histograms are only for merging. 
  /*
  TH1F *hppb0 = new TH1F("hppb0","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppb1 = new TH1F("hppb1","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppb2 = new TH1F("hppb2","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppb3 = new TH1F("hppb3","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppb4 = new TH1F("hppb4","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppb5 = new TH1F("hppb5","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppb6 = new TH1F("hppb6","",nbins_yaxian,boundaries_yaxian);
  TH1F *hppbComb = new TH1F("hppbComb","",nbins_yaxian,boundaries_yaxian);
  TH1F* hMB = new TH1F("hMB","",nbins_yaxian,boundaries_yaxian);
  */
  TH1F *hppb0 = new TH1F("hppb0","",1000,0,1000);
  TH1F *hppb1 = new TH1F("hppb1","",1000,0,1000);
  TH1F *hppb2 = new TH1F("hppb2","",1000,0,1000);
  TH1F *hppb3 = new TH1F("hppb3","",1000,0,1000);
  TH1F *hppb4 = new TH1F("hppb4","",1000,0,1000);
  TH1F *hppb5 = new TH1F("hppb5","",1000,0,1000);
  TH1F *hppb6 = new TH1F("hppb6","",1000,0,1000);
  TH1F *hppbComb = new TH1F("hppbComb","",1000,0,1000);
  TH1F* hMB = new TH1F("hMB","",1000,0,1000);
  

  //define histograms which are used for trigger turn on curves. 
  TH1F* hHLT_120 = new TH1F("hHLT_120","",1000,0,200);
  TH1F* hHLT_100 = new TH1F("hHLT_100","",1000,0,200);
  TH1F* hHLT_80 = new TH1F("hHLT_80","",1000,0,200);
  TH1F* hHLT_60 = new TH1F("hHLT_60","",1000,0,200);
  TH1F* hHLT_40 = new TH1F("hHLT_40","",1000,0,200);
  TH1F* hHLT_20 = new TH1F("hHLT_20","",1000,0,200);
  
  TH1F* hMB_0 = new TH1F("hMB_0","",1000,0,200);
  TH1F* hMB_1 = new TH1F("hMB_1","",1000,0,200);
  TH1F* hMB_2 = new TH1F("hMB_2","",1000,0,200);
  TH1F* hMB_lead = new TH1F("hMB_lead","",1000,0,200);
  TH1F* hMB_1_lead = new TH1F("hMB_1_lead","",1000,0,200);
  TH1F* hMB_2_lead = new TH1F("hMB_2_lead","",1000,0,200);
  TH1F* hBase_100 = new TH1F("hBase_100","",1000,0,200);
  TH1F* hBase_40 = new TH1F("hBase_40","",1000,0,200);

  
  Float_t MB_count = (Float_t)jetppb0_v2->GetEntries();
  cout<<"MB_count = "<<MB_count<<endl; //2.89687e07
    
  //the following is how the prescl factor was calculated in 12-003 
  Float_t presclppb0 = (Float_t)jetppb2_v2->GetEntries("jet100&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jetMB&&jet100&&raw>20&&run>210658");
  Float_t presclppb1 = (Float_t)jetppb2_v2->GetEntries("jet100&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet20&&jet100&&raw>20&&run>210658");
  Float_t presclppb2 = (Float_t)jetppb2_v2->GetEntries("jet100&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet40&&jet100&&raw>20&&run>210658");
  Float_t presclppb3 = (Float_t)jetppb2_v2->GetEntries("jet100&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet60&&jet100&&raw>20&&run>210658");
  Float_t presclppb4 = (Float_t)jetppb2_v2->GetEntries("jet100&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet80&&jet100&&raw>20&&run>210658");
  //Float_t presclppb5 = (Float_t)jetppb1_v2->GetEntries("jet100")/jetppb1_v2->GetEntries("jet100&&jet100");
  
  //using a different way to calculate the prescl factor using the prescl - this i understand now is wrong. 
  /*
  Float_t Presclppb0 = (Float_t)jetppb2_v2->GetEntries("jetMB_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jetMB&&raw>20&&run>210658");
  Float_t Presclppb1 = (Float_t)jetppb2_v2->GetEntries("jet20_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet20&&raw>20&&run>210658");
  Float_t Presclppb2 = (Float_t)jetppb2_v2->GetEntries("jet40_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet40&&raw>20&&run>210658");
  Float_t Presclppb3 = (Float_t)jetppb2_v2->GetEntries("jet60_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet60&&raw>20&&run>210658");
  Float_t Presclppb4 = (Float_t)jetppb2_v2->GetEntries("jet80_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet80&&raw>20&&run>210658");
  Float_t Presclppb5 = (Float_t)jetppb2_v2->GetEntries("jet100_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet100&&raw>20&&run>210658");
  Float_t Presclppb6 = (Float_t)jetppb2_v2->GetEntries("jet120_p&&raw>20&&run>210658")/jetppb2_v2->GetEntries("jet120&&raw>20&&run>210658");
  */
  cout<<"ppb prescl0 = "<<presclppb0<<endl;//2202.28
  cout<<"ppb prescl1 = "<<presclppb1<<endl;//1163.58
  cout<<"ppb prescl2 = "<<presclppb2<<endl;//40.9373
  cout<<"ppb prescl3 = "<<presclppb3<<endl;//4.46062
  cout<<"ppb prescl4 = "<<presclppb4<<endl;//1.41429

  /*
  cout<<"ppb Prescl0 = "<<Presclppb0<<endl;//
  cout<<"ppb Prescl1 = "<<Presclppb1<<endl;//
  cout<<"ppb Prescl2 = "<<Presclppb2<<endl;//
  cout<<"ppb Prescl3 = "<<Presclppb3<<endl;//
  cout<<"ppb Prescl4 = "<<Presclppb4<<endl;//
  cout<<"ppb Prescl5 = "<<Presclppb5<<endl;//
  cout<<"ppb Prescl6 = "<<Presclppb6<<endl;//
  */
  /*
    When we have the file with ppb and pbp
    prescl0 = 2202.28
    prescl1 = 1163.58
    prescl2 = 40.9373
    prescl3 = 4.46062
    prescl4 = 1.41429

    when we remove the pbp from the above files we are left with different numbers 
    prescl0 = 3064.5
    prescl1 = 2097.69
    prescl2 = 38.3544
    prescl3 = 4.24566
    prescl4 = 1.31221

    Please change them and the Num of equivalent MB events. 

    after adding Yaxian's extra selection cuts the presl value was changed for the pPb files
    MB_count = 125552
    ppb prescl0 = 3065.45 , the one found with trial and error when comparing yaxian's spectra is 2346.71 - closest. 
    ppb prescl1 = 2107.11
    ppb prescl2 = 38.0767
    ppb prescl3 = 4.21271
    ppb prescl4 = 1.3074
    ppb Prescl0 = 2846.71
    ppb Prescl1 = 1922.5
    ppb Prescl2 = 33.9329
    ppb Prescl3 = 3.63636
    ppb Prescl4 = 1.1013
    ppb Prescl5 = 2.55614
    ppb Prescl6 = 6.47722

    after including Yaxian's selection cuts and the good run cut (should i include it???)
    ppb prescl0 = 3407.08
    ppb prescl1 = 2197.18
    ppb prescl2 = 40.2
    ppb prescl3 = 3.82987
    ppb prescl4 = 1.38968
    ppb Prescl0 = 3285.08
    ppb Prescl1 = 2001.8
    ppb Prescl2 = 35.063
    ppb Prescl3 = 3.19709
    ppb Prescl4 = 1.12841
    ppb Prescl5 = 2.46404
    ppb Prescl6 = 6.2406

  */

  //the way to estimate the equivalent number of min bias events 
  // Look at the the total number of HLT_100 events in the HLT file and the corresponding no in the min bias file. Then simply scale that number up accordingly. simple direct proportion. 
  
  //Float_t N_MB_ppb = 6.53104e10;
  //Float_t N_MB_ppb = 3.135e10;//this one is for the ppb file only //3.265e10 

  //the one from the minbias counter which is located in cern /afs/cern.ch/user/r/rkunnawa/WORK/MinBiasCounter/CMSSW_5_3_8_HI_patch2/src/mbCounter/MinBiasCounter/test 
  // and this one only uses run number from 210676 till 211631 (which is confusing because it uses both forward and backward runs. )
  /*
  Run numbers: pPb run 2013 
  210498 <= run < 210676 old runs
  210676 <= run <= 211256 new good runs
  211313 <= run <= 211631 reverse new good runs 
  */

  //Float_t N_MB_ppb = 2.659675e10; //from the merged_MinBiasCentrality_Histo_NoPileup.root
  Float_t N_MB_ppb = 2.6026e10; //from the merged_MinBiasCentrality_Histo.root


  //read the required histograms including the prescale factor 
  //removing the chMax/jtpt>0.01 cut for comparision with Yaxian. 
  // also have to add the required run>210658 cut because those are old runs - look at https://twiki.cern.ch/twiki/bin/viewauth/CMS/PAData2013
  // that run cut removes about 30% of the data. 29.12% for the HLT_100 file and 33.42% for the MB file. 

  TCut ppb0 = "abs(eta_CM)<1&&jetMB&&!jet80&&!jet100&&raw>20&&run>210658";
  TCut ppb1 = "abs(eta_CM)<1&&jet20&&!jet80&&!jet100&&raw>20&&run>210658";
  TCut ppb2 = "abs(eta_CM)<1&&jet40&&!jet80&&!jet100&&raw>20&&run>210658";
  TCut ppb3 = "abs(eta_CM)<1&&jet60&&!jet80&&!jet100&&raw>20&&run>210658";
  TCut ppb4 = "abs(eta_CM)<1&&jet80&&!jet100&&raw>20&&run>210658";//try with or without 60 first and then try without 60
  TCut ppb5 = "abs(eta_CM)<1&&jet100&&raw>20&&run>210658";
  //my ntuples already had the event selection cuts applied in them. 
  //add the required event selection cuts here. 
  TCut eventcut = "fabs(vz)<15&&pHBHENoiseFilter&&pprimaryvertexFilter&&pPAcollisionEventSelectionPA&&pVertexFilterCutGplus";
  
  //the value multiplying the ppb0 here is hard coded from the prescl value got from above. We can also do it with ""*ppb0
  
  jetppb0_v2->Project("hppb0","pt","3407.08"*ppb0)
    //jetppb0_v2->Project("hppb0","pt","3407.08"*(ppb0 && eventcut)); //if you are running on ntuples produced by my code then we dont need this evetcut. keep it for running on hiForest. 
  hppb0->Print("base");

  //jetppb1_v2->Project("hppb1","pt","1163.58"*ppb1);
  //hppb1->Print("base");

  //jetppb1_v2->Project("hppb2","pt","40.9373"*ppb2);
  //hppb2->Print("base");

  //jetppb1_v2->Project("hppb3","pt","4.46062"*ppb3);
  //hppb3->Print("base");

  jetppb2_v2->Project("hppb4","pt",ppb4);
  //jetppb2_v2->Project("hppb4","pt",ppb4&&eventcut);
  hppb4->Print("base");
  //divideBinWidth(hppb4);
  //hppb4->Scale(1./5.25585e10);
  //hppb4->Scale(1./4);

  jetppb2_v2->Project("hppb5","pt",ppb5);
  //jetppb2_v2->Project("hppb5","pt",ppb5&&eventcut);
  hppb5->Print("base");
  //divideBinWidth(hppb5);
  //hppb5->Scale(1./6.01453e10); //total no of min bias events for 80
  //hppb5->Scale(1./4);

  // Add the histograms
  hppbComb->Add(hppb0,1);
  //hppbComb->Add(hppb1,1);
  //hppbComb->Add(hppb2,1);
  //hppbComb->Add(hppb3,1);
  hppbComb->Add(hppb4,1);
  hppbComb->Add(hppb5,1);
  //hppbComb->Add(hppb6,1);

  hppbComb->Print("base");

  

  //plot the trigger turn-on curves 
  //TCut turnon_MB = "abs(eta_CM)<1&&jetMB&&chMax/pt>0.01";
  //TCut turnon_20 = "abs(eta_CM)<1&&jet20&&chMax/pt>0.01";
  //TCut turnon_40 = "abs(eta_CM)<1&&jet40&&chMax/pt>0.01";
  //TCut turnon_60 = "abs(eta_CM)<1&&jet60&&chMax/pt>0.01";
  //TCut turnon_80 = "abs(eta_CM)<1&&jet80&&chMax/pt>0.01";
  //TCut turnon_100 = "abs(eta_CM)<1&&jet100&&chMax/pt>0.01";
  //TCut turnon_120 = "abs(eta_CM)<1&&jet120&&chMax/pt>0.01";
  
  TCut turnon_MB = "jetMB";
  TCut turnon_20 = "jet20";
  TCut turnon_40 = "jet40";
  TCut turnon_60 = "jet60";
  TCut turnon_80 = "jet80";
  TCut turnon_100 = "jet100";
  //TCut turnon_120 = "jet120";
  
  jetppb0_v2->Project("hMB","pt","3064.5"*turnon_MB);
  hMB->Print("base");
  //jetppb0_v2->Project("hMB_2","pt","abs(eta_CM)<1&&chMax/pt>0.01");
  //hMB_2->Print("base");
  /*
  jetppb1_v2->Project("hHLT_20","pt","1163.58"*turnon_20);
  hHLT_20->Print("base");
  jetppb1_v2->Project("hHLT_40","pt","40.9373"*turnon_40);
  hHLT_40->Print("base");
  jetppb1_v2->Project("hHLT_60","pt","4.46062"*turnon_60);
  hHLT_60->Print("base");
  jetppb2_v2->Project("hHLT_80","pt","1.41429"*turnon_80);
  //jetppb2_v2->Project("hHLT_80","pt",turnon_80);
  hHLT_80->Print("base");
  jetppb2_v2->Project("hHLT_100","pt","1"*turnon_100);
  hHLT_100->Print("base");
  jetppb2_v2->Project("hBase_100","pt");
  hBase_100->Print("base");
  */  

  /*
  jetppb1_v2->Project("hBase_40","pt");
  hBase_40->Print("base");
  //jetppb2_v2->Project("h120","pt",""*turnon_120);
  

  jetppb0_v2->Project("hMB_lead","ptlead","jetMB");
  jetppb1_v2->Project("hMB_1_lead","ptlead","jetMB");
  jetppb2_v2->Project("hMB_2_lead","ptlead","jetMB");
  jetppb0_v2->Project("hMB","pt","jetMB");
  jetppb1_v2->Project("hMB_1","pt","jetMB");
  jetppb2_v2->Project("hMB_2","pt","jetMB");
  jetppb0_v2->Project("hHLT_20","ptlead","jetMB&&jet20");
  jetppb0_v2->Project("hHLT_40","ptlead","jetMB&&jet40");
  jetppb1_v2->Project("hHLT_60","ptlead","jetMB&&jet60");
  jetppb1_v2->Project("hHLT_80","ptlead","jetMB&&jet80");
  //jetppb2_v2->Project("hHLT_100","ptlead","jetMB&&jet100");
  
  jetppb2_v2->Project("hHLT_100","ptlead","1"*turnon_100);
  hHLT_100->Print("base");
  jetppb2_v2->Project("hBase_100","ptlead");
  hBase_100->Print("base");

  TH1F* hTurnon_20 = (TH1F*)hHLT_20->Clone("hTurnon_20");
  TH1F* hTurnon_40 = (TH1F*)hHLT_40->Clone("hTurnon_40");
  TH1F* hTurnon_60 = (TH1F*)hHLT_60->Clone("hTurnon_60");
  TH1F* hTurnon_80 = (TH1F*)hHLT_80->Clone("hTurnon_80");
  TH1F* hTurnon_100 = (TH1F*)hHLT_100->Clone("hTurnon_100");

  hTurnon_20->Divide(hMB_lead);
  hTurnon_40->Divide(hMB_lead);
  hTurnon_60->Divide(hMB_1_lead);
  hTurnon_80->Divide(hMB_1_lead);
  hTurnon_100->Divide(hBase_100);
  
  TFile fout("ppb_turnon_HLT.root","RECREATE");
  hTurnon_20->Write();
  hTurnon_40->Write();
  hTurnon_60->Write();
  hTurnon_80->Write();
  hTurnon_100->Write();
  hMB_1->Write();
  hMB_2->Write();
  hMB->Write();
  hMB_1_lead->Write();
  hMB_2_lead->Write();
  hMB_lead->Write();
  hHLT_20->Write();
  hHLT_40->Write();
  hHLT_60->Write();
  hHLT_80->Write();
  hHLT_100->Write();
  fout.Close();

  TCanvas *c1 = new TCanvas("c1","",800,600);
  hTurnon_40->SetAxisRange(0,2,"Y");
  hTurnon_40->SetAxisRange(0,140,"X");
  hTurnon_40->SetTitle("HLT Turn on curves");
  hTurnon_40->SetXTitle("p_{T} GeV/c");
  hTurnon_40->SetMarkerStyle(24);
  hTurnon_40->SetMarkerColor(2);
  //hTurnon_20->SetMarkerStyle(25);
  //hTurnon_20->SetMarkerColor(3);
  hTurnon_60->SetMarkerStyle(26);
  hTurnon_60->SetMarkerColor(4);
  hTurnon_80->SetMarkerStyle(27);
  hTurnon_80->SetMarkerColor(1);
  hTurnon_100->SetMarkerStyle(28);
  hTurnon_100->SetMarkerColor(6);
  hTurnon_40->Draw();
  //hTurnon_20->Draw("same");
  hTurnon_60->Draw("same");
  hTurnon_80->Draw("same");
  hTurnon_100->Draw("same");
  
  TLegend *legend = myLegend(0.34,0.65,0.65,0.9);
  //legend->AddEntry(hTurnon_20,"HLT_20","pl");
  legend->AddEntry(hTurnon_40,"HLT_40","pl");
  legend->AddEntry(hTurnon_60,"HLT_60","pl");
  legend->AddEntry(hTurnon_80,"HLT_80","pl");
  legend->AddEntry(hTurnon_100,"HLT_100","pl");
  legend->SetTextSize(0.04);
  legend->Draw();
  
  drawText("PPb 2013 Data Prompt Reco",0.15,0.8,20);
  drawText("Anti-k_{T} PU PF Jets R = 0.3,|vz|<15",0.15,0.7,20);

  c1->SaveAs("ppb_2013_HLT_turn_on_curve.jpg","RECREATE");
  */
  //
  TCanvas* c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();
  TH1F* hPPbComb = rebin(hppbComb,"hPPbComb");
  TH1F* hPPbComb_2 = (TH1F*)hppbComb->Clone("hPPbComb_2");
  TH1F* hPPbComb_3 = (TH1F*)hppbComb->Clone("hPPbComb_3");
  //hPPbComb->Scale(1./3.083e11);
  hPPbComb->Print("base");
  hPPbComb_2->Print("base");
  hPPbComb_3->Print("base");

  //hPPbComb->Scale(1./4);
  divideBinWidth(hPPbComb);

  hppbComb->SetMarkerStyle(20);
  hppbComb->SetMarkerColor(kBlack);
  //hPPbComb->SetYTitle("#frac{dN}{N_{MB} dp_{T} d#eta}");
  hppbComb->SetYTitle("counts");
  hppbComb->SetXTitle("Jet p_{T} GeV/c");
  hppbComb->Draw();
  hppb0->SetMarkerStyle(27);
  hppb0->SetMarkerColor(28);
  hppb0->Draw("same");
  //hppb3->SetMarkerStyle(25);
  //hppb3->SetMarkerColor(kBlue);
  //hppb3->Draw("same");
  hppb4->SetMarkerStyle(26);
  hppb4->SetMarkerColor(kRed);
  hppb4->Draw("same");
  hppb5->SetMarkerStyle(24);
  hppb5->SetMarkerColor(kGreen);
  hppb5->Draw("same");

  TLegend *title = myLegend(0.34,0.65,0.65,0.9);
  title->AddEntry(hppbComb,"PPb Merged","pl");
  title->AddEntry(hppb0,"PPb MB","pl");
  //title->AddEntry(hppb3,"PPb HLT_60","pl");
  title->AddEntry(hppb4,"PPb HLT_80","pl");
  title->AddEntry(hppb5,"PPb HLT_100","pl");
  title->SetTextSize(0.04);
  title->Draw();
  drawText("PPb 2013 Data PromptReco",0.3,0.65,20);  
  drawText("Anti-k_{T} PU PF Jets R = 0.3, |#eta_{CM}|<1, |vz|<15",0.3,0.56,20);

  //c2->SaveAs("ppb_2013_pt_merged_60_lowest.gif","RECREATE");
  c2->SaveAs("ppb_2013_pt_merged_MB_lowest_eta_CM_1_lowest.jpg","RECREATE");

  //comparison with the unfolded data: 
  TFile* fPPb_Unfo = TFile::Open("result-2013-ppb-akPu3PF-cent-1/ppb_merge_MB_eta_CM_1_lowest_pp_mc_Unfo_akPu3PF_cent_1.root");
  TH1F* hPPb_Unfo = (TH1F*)fPPb_Unfo->Get("Unfolded_cent0");
  //TH1F* hPPb_Unfo = rebin(hPPb_Unfo_test,"hPPb_Unfo");
  hPPb_Unfo->SetAxisRange(20,500,"X");
  TH1F* hPPb_Unfo_2 = (TH1F*)hPPb_Unfo->Clone("hPPb_Unfo_2");
  TH1F* hPPb_Unfo_3 = (TH1F*)hPPb_Unfo->Clone("hPPb_Unfo_3");
  TH1F* hPPb_Unfo_BinByBin = (TH1F*)fPPb_Unfo->Get("UnfoldedBinByBin_cent0");
  //TH1F* hPPb_Unfo_BinByBin = rebin(hPPb_Unfo_BinByBin_test,"hPPb_Unfo_BinByBin");
  hPPb_Unfo_BinByBin->SetAxisRange(30,500,"X");
  hPPb_Unfo_BinByBin->Print("base");
  hPPb_Unfo->Print("base");
  hPPb_Unfo_2->Print("base");
  hPPb_Unfo_3->Print("base");
  divideBinWidth(hPPb_Unfo);
  divideBinWidth(hPPb_Unfo_BinByBin);
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetLogy();
  hPPbComb->SetAxisRange(20,500,"X");
  hPPbComb->SetYTitle("counts");
  hPPbComb->SetXTitle("Jet p_{T} GeV/c");
  hPPbComb->SetTitle("pPb 2013 Merged");
  hPPbComb->SetMarkerStyle(7);
  hPPbComb->SetMarkerColor(kBlack);
  
  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]/pow(x,[1])");
  hPPbComb->Fit("fPowerLaw","","",30,500);
  hPPbComb->Fit("fPowerLaw","","",30,500);
  hPPbComb->Fit("fPowerLaw","","",30,500);
  hPPbComb->Fit("fPowerLaw","","",30,500);
  hPPbComb->Fit("fPowerLaw","","",30,500);
  
  hPPbComb->Draw();
  hPPb_Unfo->SetMarkerStyle(24);
  hPPb_Unfo->SetMarkerColor(kBlack);
  hPPb_Unfo->Draw("same");
  hPPb_Unfo_BinByBin->SetMarkerStyle(26);
  hPPb_Unfo_BinByBin->SetMarkerColor(kRed);
  hPPb_Unfo_BinByBin->Draw("same");

  TLegend *title1 = myLegend(0.34,0.65,0.65,0.9);
  title1->AddEntry(hPPbComb,"Measured","pl");
  title1->AddEntry(hPPb_Unfo,"Bayesian","pl");
  title1->AddEntry(hPPb_Unfo_BinByBin,"Bin by Bin","pl");
  title1->SetTextSize(0.04);
  title1->Draw();
  //c3->SaveAs("ppb_merged_60_lowest_measured_vs_unfolded.gif","RECREATE");
  c3->SaveAs("ppb_merged_MB_eta_cm_1_lowest_measured_vs_unfolded_correct_prescl.jpg","RECREATE");
  //c3->SaveAs("ppb_merged_40_eta_cm_1_lowest_measured_vs_unfolded.jpg","RECREATE");

  //create the spectra with event fraction. //preliminary changing from ppb number 6.53104e10 to half that value
  hPPbComb_2->Scale(1./N_MB_ppb);
  hPPbComb_2->Scale(1./2);
  hPPb_Unfo_2->Scale(1./N_MB_ppb);
  hPPb_Unfo_2->Scale(1./2);
  hMB->Scale(1./3.265e10);// i divided 2.89687e07 by 2 to get 1.448e07, now trying to use that other number from above
  //hMB->Scale(1./2.89687e07)
  hMB->Scale(1./2);
  divideBinWidth(hMB);
  divideBinWidth(hPPbComb_2);
  divideBinWidth(hPPb_Unfo_2);


 

  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetLogy();
  hPPbComb_2->SetAxisRange(30,500,"X");
  hPPbComb_2->SetYTitle("#frac{dN}{N_{MB} dp_{T} d#eta}");
  hPPbComb_2->SetXTitle("Jet p_{T} GeV/c");
  hPPbComb_2->SetMarkerStyle(7);
  hPPbComb_2->SetMarkerColor(kBlack);
  
  
  hMB->SetMarkerStyle(26);
  hMB->SetMarkerColor(kRed);
  hMB->Draw();
  hPPbComb_2->Draw("same");
  hPPb_Unfo_2->SetMarkerStyle(24);
  hPPb_Unfo_2->SetMarkerColor(kBlack);
  hPPb_Unfo_2->Draw("same");

  TLegend *title2 = myLegend(0.34,0.65,0.65,0.9);
  title2->AddEntry(hPPbComb_2,"PPb(merged) 2013 Measured","pl");
  title2->AddEntry(hMB,"PPb MB","pl");
  title2->AddEntry(hPPb_Unfo_2,"PPb(merged) 2013 Bayesian","pl");
  //title2->AddEntry(h_Yaxian,"Trigger Comb Spectra Yaxian","pl")
  title2->SetTextSize(0.04);
  title2->Draw();
  //c4->SaveAs("ppb_merged_60_lowest_evtfrac_measured_vs_unfolded.gif","RECREATE");  
  c4->SaveAs("ppb_merged_MB_eta_CM_1_lowest_evtfrac_measured_vs_unfolded.jpg","RECREATE");  
  //c4->SaveAs("ppb_merged_40_eta_CM_1_lowest_evtfrac_measured_vs_unfolded.jpg","RECREATE");  
  /*
  //create the spectra with cross section: 
  hPPbComb_3->Scale(1./30.9);
  hPPbComb_3->Scale(1./2);
  hPPb_Unfo_3->Scale(1./30.9);
  hPPb_Unfo_3->Scale(1./2);
  divideBinWidth(hPPbComb_3);
  divideBinWidth(hPPb_Unfo_3);

  TCanvas *c5 = new TCanvas("c5","",800,600);
  c5->SetLogy();
  hPPbComb_3->SetAxisRange(30,500,"X");
  hPPbComb_3->SetYTitle("#sigma (nb)");
  hPPbComb_3->SetXTitle("Jet p_{T} GeV/c");
  hPPbComb_3->SetMarkerStyle(7);
  hPPbComb_3->SetMarkerColor(kBlack);
  hPPbComb_3->Draw();
  hPPb_Unfo_3->SetMarkerStyle(24);
  hPPb_Unfo_3->SetMarkerColor(kBlack);
  hPPb_Unfo_3->Draw("same");

  TLegend *title3 = myLegend(0.34,0.65,0.65,0.9);
  title3->AddEntry(hPPbComb_2,"PPb(merged) 2013 Measured","pl");
  title3->AddEntry(hPPb_Unfo_2,"PPb(merged) 2013 Bayesian","pl");
  title3->SetTextSize(0.04);
  title3->Draw();
  //c5->SaveAs("ppb_merged_60_lowest_crosssection_measured_vs_unfolded.gif","RECREATE");  
  c5->SaveAs("ppb_merged_MB_eta_CM_1_lowest_crosssection_measured_vs_unfolded.jpg","RECREATE");  
  //c5->SaveAs("ppb_merged_40_eta_CM_1_lowest_crosssection_measured_vs_unfolded.jpg","RECREATE");  
  */

  
  //compare the MC reco and Gen for pPb used to create the unfolding matrix. 
  TH1F* hppbMC_Gen = (TH1F*)fPPb_Unfo->Get("hGen_cent0");
  TH1F* hppbMC_Reco = (TH1F*)fPPb_Unfo->Get("hRecoMC_cent0");
  hppbMC_Gen->Print("base");
  hppbMC_Reco->Print("base");
  
  TCanvas *c6 = new TCanvas("c6","",800,600);
  c6->SetLogy();
  hppbMC_Gen->SetTitle("");
  hppbMC_Gen->SetYTitle("arbitrary units");
  hppbMC_Gen->SetXTitle("Jet p_{T} GeV/c");
  hppbMC_Gen->SetMarkerStyle(22);
  hppbMC_Gen->SetMarkerColor(kBlack);
  hppbMC_Gen->Draw();
  hppbMC_Reco->SetMarkerStyle(23);
  hppbMC_Reco->SetMarkerColor(kRed);
  hppbMC_Reco->Draw("same");
  
  TLegend *title4 = myLegend(0.34,0.65,0.65,0.9);
  title4->SetTextSize(0.04);
  title4->AddEntry(hppbMC_Gen,"MC - Ref p_{T}","pl");
  title4->AddEntry(hppbMC_Reco,"MC - Jet p_{T}","pl");
  title4->Draw();

  drawText("PPb 2013 MC",0.3,0.65,20);
  drawText("|eta_CM|<1  akPu3PF",0.3,0.75,20);
  c6->SaveAs("ppb_2013_mc_reco_gen_eta_CM_1.jpg","RECREATE");

  //compare with Yaxian's epectra
  TFile *fYaxian = TFile::Open("PPbJetTrigHFsumEta4Bin1PYTHIAAkPu3PFJetSpectra.root");
  TH1F* h_Yaxian = fYaxian->Get("DataJetInEtaBin-10_10;1");
  //TH1F* h_Yaxian_evenbin = rebin(h_Yaxian,"h_Yaxian_evenbin");
  h_Yaxian->Print("base");
  divideBinWidth(h_Yaxian);
  h_Yaxian->Scale(1./2);
  //h_Yaxian_evenbin->Print("base");
  //divideBinWidth(h_Yaxian_evenbin);
  //h_Yaxian_evenbin->Scale(1./2);

  TCanvas *c7 = new TCanvas("c7","",800,600);
  c7->Divide(2,1);
  c7->cd(1)->SetLogy();
  h_Yaxian->SetMarkerStyle(25);
  h_Yaxian->SetMarkerColor(kBlue);
  h_Yaxian->SetYTitle("#frac{dN}{N_{MB} dp_{T} d#eta}");
  h_Yaxian->SetXTitle("Jet p_{T} GeV/c");
  h_Yaxian->Draw();
  //h_Yaxian->Draw();
  TH1F* hPPbComb_y_test = hppbComb->Clone("hPPbComb_y_test");
  TH1F* hPPbComb_y = rebin_yaxian(hPPbComb_y_test,"hPPbComb_y");
  divideBinWidth(hPPbComb_y);
  hPPbComb_y->Scale(1./2);
  hPPbComb_y->Scale(1./N_MB_ppb);
  hPPbComb_y->SetMarkerStyle(23);
  hPPbComb_y->SetMarkerColor(kRed);
  hPPbComb_y->Draw("same");

  //divideBinWidth(hPPbComb_y_test);
  //hPPbComb_y_test->Scale(1./2);
  //hPPbComb_y_test->Scale(1./N_MB_ppb);
  //hPPbComb_y_test->SetMarkerStyle(24);
  //hPPbComb_y_test->SetMarkerColor(kGreen);
  //hPPbComb_y_test->Draw("same");

  TLegend *title5 = myLegend(0.34,0.65,0.65,0.9);
  title5->SetTextSize(0.04);
  title5->AddEntry(h_Yaxian,"Yaxian","pl");
  title5->AddEntry(hPPbComb_y,"Raghav","pl");
  //title5->AddEntry(hPPbComb_y_test,"Raghav even bins","pl");
  title5->Draw();
  drawText("|eta_CM|<1, akPu3PF",0.3,0.75,20);

  c7->cd(2);
  TH1F* hComparison = (TH1F*)h_Yaxian->Clone("hComparison");
  hComparison->Print("base");
  hPPbComb_y->Print("base");
  hComparison->Divide(hPPbComb_y);
  hComparison->SetYTitle(" ");
  hComparison->Draw();

  c7->SaveAs("yaxian_ppb_spectra_comparison_correct_prescl.gif","RECREATE");

  //calculate RpA here as well along with Yaxian's 

 //get yaxian's mc spectra 
  TH1F* hPP_Refe_Yaxian = (TH1F*)fYaxian->Get("PYTHIAJetInEtaBin-10_10_Cen0-100%;1");
  //TH1F* hPP_Refe_Yaxian_evenbin = rebin(hPP_Refe_Yaxian,"hPP_Refe_Yaxian_evenbin");

  TH1F* hPP_Refe = (TH1F*)fPPb_Unfo->Get("hRecoMC_cent1");//change it from hRecoMC_cent1 and hGen_cent
  TH1F* hPP_Refe_y = rebin_yaxian(hPP_Refe,"hPP_Refe_y");
  cout<<"hPP_Refe"<<endl;
  hPP_Refe->Print("base");
  hPP_Refe->Scale(1./70);
  hPP_Refe->Scale(6.9);
  divideBinWidth(hPP_Refe);//p_T width
  hPP_Refe->Scale(1./2);//RAPIDITY

  cout<<"hPP_Refe_y"<<endl;
  hPP_Refe_y->Print("base");
  hPP_Refe_y->Scale(1./70);
  hPP_Refe_y->Scale(6.9);
  divideBinWidth(hPP_Refe_y);//p_T width
  hPP_Refe_y->Scale(1./2);//RAPIDITY

 

  hPP_Refe_Yaxian->Print("base");
  hPP_Refe_Yaxian->Scale(1./70);
  hPP_Refe_Yaxian->Scale(1./2);
  hPP_Refe_Yaxian->Scale(6.9);
  divideBinWidth(hPP_Refe_Yaxian);

  //hPP_Refe_Yaxian_evenbin->Print("base");
  //hPP_Refe_Yaxian_evenbin->Scale(1./70);
  //hPP_Refe_Yaxian_evenbin->Scale(1./2);
  //hPP_Refe_Yaxian_evenbin->Scale(6.9);
  //divideBinWidth(hPP_Refe_Yaxian_evenbin);
  
  //compare the pp reference mine vs yaxian's
  TCanvas *c9 = new TCanvas("c9","",800,600);
  c9->Divide(2,1);
  c9->cd(1);
  c9->cd(1)->SetLogy();

  hPP_Refe_Yaxian->SetMarkerStyle(25);
  hPP_Refe_Yaxian->SetMarkerColor(kBlue);
  hPP_Refe_Yaxian->SetTitle("PP Reference comparison");
  hPP_Refe_Yaxian->Draw();
  hPP_Refe_y->SetMarkerStyle(23);
  hPP_Refe_y->SetMarkerColor(kRed);
  hPP_Refe_y->Draw("same");

  TLegend *title7 = myLegend(0.5,0.35,0.65,0.5);
  title7->SetTextSize(0.04);
  title7->AddEntry(hPP_Refe_Yaxian,"Yaxian","pl");
  title7->AddEntry(hPP_Refe_y,"Raghav","pl");
  title7->Draw();
  drawText("|eta_CM|<1, PP MC ak3PF",0.5,0.5,20);
  
  c9->cd(2);
  TH1F* hPP_comparison = (TH1F*)hPP_Refe_Yaxian->Clone("hPP_comparison");
  hPP_comparison->Divide(hPP_Refe_y);
  hPP_comparison->Draw();

  c9->SaveAs("yaxian_pp_refe_comparison.gif","RECREATE");



  /*
  TCanvas *c10 = new TCanvas("c10","",800,600);
  c10->Divide(2,1);
  c10->cd(1)->SetLogy();
  h_Yaxian_evenbin->Draw();
  h_Yaxian->Draw("same");

  c10->cd(2)->SetLogy();
  hPP_Refe_Yaxian_evenbin->Draw();
  hPP_Refe_Yaxian->Draw("same");
  
  c10->SaveAs("yaxian_diff_bin_comparision.gif","RECREATE");
  */

  
  //TH1F* hRpA_Meas_r_even = (TH1F*)hPPbComb_y_test->Clone("hRpA_Meas_r_even");
  //cout<<"hRpA_Meas_r_even"<<endl;
  //hRpA_Meas_r_even->Print("base");
  //TH1F* hRpA_Unfo_r_even = (TH1F*)hPPb_Unfo_2->Clone("hRpA_Unfo_r_even");
  //cout<<"hRpA_Unfo_r_even"<<endl;
  //hRpA_Unfo_r_even->Print("base");
  //multiplyBinWidth(hRpA_Unfo_r_even);
  TH1F* hRpA_Meas_r = (TH1F*)hPPbComb_y->Clone("hRpA_Meas_r");
  cout<<"hRpA_Meas_r"<<endl;
  hRpA_Meas_r->Print("base");
  TH1F* hRpA_Meas_y = (TH1F*)h_Yaxian->Clone("hRpA_Meas_y");
  cout<<"hRpA_Meas_y"<<endl;
  hRpA_Meas_y->Print("base");
  //TH1F* hRpA_Unfo_r = rebin_yaxian(hRpA_Unfo_r_even,"hRpA_Unfo_r");
  //cout<<"hRpA_Unfo_r"<<endl;
  //hRpA_Unfo_r->Print("base");
  
  //hRpA_Meas_r_even->Divide(hPP_Refe);
  //hRpA_Unfo_r_even->Divide(hPP_Refe);
  hRpA_Meas_r->Divide(hPP_Refe_y);
  hRpA_Meas_y->Divide(hPP_Refe_Yaxian);
  //hRpA_Unfo_r->Divide(hPP_Refe_y);

  TCanvas *c8 = new TCanvas("c8","",800,600);
  //hRpA_Meas_r_even->Draw();
  //hRpA_Unfo_r_even->Draw("same");
  hRpA_Meas_r->SetTitle("Measured RpA - pPb");
  hRpA_Meas_r->Draw();
  hRpA_Meas_y->Draw("same");
  //hRpA_Unfo_r->Draw("same");

  TLegend *title6 = myLegend(0.5,0.35,0.65,0.5);
  title6->SetTextSize(0.04);
  title6->AddEntry(hRpA_Meas_y,"Yaxian","pl");
  title6->AddEntry(hRpA_Meas_r,"Raghav","pl");
  title6->Draw();
  drawText("|eta_CM|<1, akPu3PF",0.5,0.5,20);
  
  c8->SaveAs("RpA_multiple_ppb_pbp.gif");

  


  //TFile f("merge_ppb_60_lowest_HLT_V2.root","RECREATE");
  //TFile f("merge_ppb_40_eta_CM_1_lowest_HLT_V2.root","RECREATE");
  TFile f("merge_ppb_MB_eta_CM_1_lowest_HLT_V2.root","RECREATE");
  hppb1->Write();
  hppb2->Write();
  hppb3->Write();
  hppb4->Write();
  hppbComb->Write();
  hPPbComb->Write();
  hPPb_Unfo->Write();
  hPPbComb_2->Write();
  hPPb_Unfo_2->Write();
  hPPbComb_3->Write();
  hPPb_Unfo_3->Write();
  hPPbComb_y->Write();
  hPP_Refe->Write();
  hPP_Refe_y->Write();
  //hPP_Refe_Yaxian_evenbin->Write();
  hPP_Refe_Yaxian->Write();
  h_Yaxian->Write();
  //h_Yaxian_evenbin->Write();
  //hRpA_Meas_r_even->Write();
  //hRpA_Unfo_r_even->Write();
  hRpA_Meas_r->Write();
  //hRpA_Unfo_r->Write();
  hRpA_Meas_y->Write();
  f.Close();
  
  



}
