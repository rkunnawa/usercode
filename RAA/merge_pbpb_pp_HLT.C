// Raghav Kunnawalkam Elayavalli
// created: March 12th 2014

// macro to read in latest PbPb data files and merge those three triggers 
// effective minbias events is calculated outside(interactively) and then it is hardcoded in the macro

// make the change to centrality  = 6. just increase it.  

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

//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

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

const Double_t boundaries_jetPtBin[]={0,4,8,14,24,34,44,54,64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 507, 592, 1000};
const int nbins_jetPtBin = 25;

static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 790, 967};

// rebin the spectra
TH1F *rebin(TH1F *h, char *histName)
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

void merge_pbpb_pp_HLT(int radius = 3, char *algo = "Vs"){
  
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  // number convension:
  // 0 - MB
  // 1 - 55 or 65
  // 2 - 80 or 95 
  // 80 is the unprescaled trigger - yes
  // 
  
  //data files - PbPb 
  TFile *fpbpb0_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/HIMinBias2011_GR_R_53_LV6_CMSSW_5_3_16_Forest_Track8_Jet21/0.root");
  TFile *fpbpb1_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet55or65_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root");
  TFile *fpbpb1 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet55or65_GR_R_53_LV6_25Mar2014_0200CET_Track8_Jet26/0.root");

  TFile *fpbpb2_old = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_12Mar2014_0000CET_Track8_Jet21/0.root");
  TFile *fpbpb2 = TFile::Open("/mnt/hadoop/cms/store/user/velicanu/hiForest_Jet80or95_GR_R_53_LV6_25Mar2014_0200CET_Track8_Jet26/0.root");


  // data files - pp 
  TFile *fpp1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet80_v2.root");
  TFile *fpp2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet40_v2.root");

  /*
  //loading the files for doing the JEC's on the fly
  string fJECL2AK3PF = "HI_PythiaZ2_2760GeV_5316_v14_L2Relative_AK3PF_offline.txt";
  string fJECL3AK3PF = "HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_AK3PF_offline.txt";
  string fJECL2AKPu3PF = "HI_PythiaZ2_2760GeV_5316_v14_L2Relative_AKPu3PF_offline.txt";
  string fJECL3AKPu3PF = "HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_AKPu3PF_offline.txt";
  string fJECL2AKVs3PF = "HI_PythiaZ2_2760GeV_5316_v14_L2Relative_AKVs3PF_offline.txt";
  string fJECL3AKVs3PF = "HI_PythiaZ2_2760GeV_5316_v14_L3Absolute_AKVs3PF_offline.txt";

  // grab the JEC's
  vector <JetCorrectorParameters> vpar_HI53x;
  FactorizedJetCorrector *JEC_HI53;

  JetCorrectorParameters *parHI53x_l2 = new JetCorrectorParameters(fJECL2AKVs3PF.c_str());
  JetCorrectorParameters *parHI53x_l3 = new JetCorrectorParameters(fJECL3AKVs3PF.c_str());

  vpar_HI53x.push_back(*parHI53x_l2);
  vpar_HI53x.push_back(*parHI53x_l3);
  JEC_HI53x = new FactorizedJetCorrector(vpar_HI53x);
  */

  // grab the trees from the data files. 
  TTree *jetpbpb0_old = (TTree*)fpbpb0_old->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb1 = (TTree*)fpbpb1->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb1_old = (TTree*)fpbpb1_old->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb2 = (TTree*)fpbpb2->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));
  TTree *jetpbpb2_old = (TTree*)fpbpb2_old->Get(Form("ak%s%dPFJetAnalyzer/t",algo,radius));

  TTree *evtpbpb0_old = (TTree*)fpbpb0_old->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpbpb1 = (TTree*)fpbpb1->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpbpb1_old = (TTree*)fpbpb1_old->Get("hiEvtAnalyzer/HiTree");
  TTree *evtpbpb2 = (TTree*)fpbpb2->Get("hiEvtAnalyzer/HiTree");  
  TTree *evtpbpb2_old = (TTree*)fpbpb2_old->Get("hiEvtAnalyzer/HiTree");

  TTree* hltpbpb0_old = (TTree*)fpbpb0_old->Get("hltanalysis/HltTree");
  TTree* hltpbpb1 = (TTree*)fpbpb1->Get("hltanalysis/HltTree");
  TTree* hltpbpb1_old = (TTree*)fpbpb1_old->Get("hltanalysis/HltTree");
  TTree* hltpbpb2 = (TTree*)fpbpb2->Get("hltanalysis/HltTree");
  TTree* hltpbpb2_old = (TTree*)fpbpb2_old->Get("hltanalysis/HltTree");

  TTree* skmpbpb0_old = (TTree*)fpbpb0_old->Get("skimanalysis/HltTree");
  TTree* skmpbpb1 = (TTree*)fpbpb1->Get("skimanalysis/HltTree");
  TTree* skmpbpb1_old = (TTree*)fpbpb1_old->Get("skimanalysis/HltTree");
  TTree* skmpbpb2 = (TTree*)fpbpb2->Get("skimanalysis/HltTree");
  TTree* skmpbpb2_old = (TTree*)fpbpb2_old->Get("skimanalysis/HltTree");

  jetpbpb0_old->AddFriend(evtpbpb0_old);
  jetpbpb1->AddFriend(evtpbpb1);
  jetpbpb1_old->AddFriend(evtpbpb1_old);
  jetpbpb2->AddFriend(evtpbpb2);
  jetpbpb2_old->AddFriend(evtpbpb2_old);

  jetpbpb0_old->AddFriend(hltpbpb0_old);
  jetpbpb1->AddFriend(hltpbpb1);
  jetpbpb1_old->AddFriend(hltpbpb1_old);
  jetpbpb2->AddFriend(hltpbpb2);
  jetpbpb2_old->AddFriend(hltpbpb2_old);

  jetpbpb0_old->AddFriend(skmpbpb0_old);
  jetpbpb1->AddFriend(skmpbpb1);
  jetpbpb1_old->AddFriend(skmpbpb1_old);
  jetpbpb2->AddFriend(skmpbpb2);
  jetpbpb2_old->AddFriend(skmpbpb2_old);


  //do it for the pp - need to check up on this. 
  TTree *jetpp1_v2 = (TTree*)fpp1_v2->Get(Form("jetR%d",radius));
  TTree *jetpp2_v2 = (TTree*)fpp2_v2->Get(Form("jetR%d",radius));

  TTree *evtpp1_v2 = (TTree*)fpp1_v2->Get("evt");
  TTree *evtpp2_v2 = (TTree*)fpp2_v2->Get("evt");

  jetpp1_v2->AddFriend(evtpp1_v2);
  jetpp2_v2->AddFriend(evtpp2_v2);

  //get all the pp spectra here: 
  TCut pp3 = "abs(eta)<2&&jet40&&!jet60&&!jet80&&chMax/pt>0.01";
  
  TH1F *hpp1 = new TH1F("hpp1","",1000,0,1000);
  TH1F *hpp2 = new TH1F("hpp2","",1000,0,1000);
  TH1F *hpp3 = new TH1F("hpp3","",1000,0,1000);
  TH1F *hppComb = new TH1F("hppComb","",1000,0,1000);
  
  //get the prescl factor information. 
  //Float_t presclpbpb3 = (Float_t)jetpbpb1_v2->GetEntries("jet80")/jetpbpb1_v2->GetEntries("jet55&&jet80");
  //cout<<"pbpb prescl3 = "<<presclpbpb3<<endl;//1.99871
  //Float_t presclpp3 = (Float_t)jetpp1_v2->GetEntries("jet80")/jetpp1_v2->GetEntries("jet40&&jet80");
  //cout<<"pp prescl3 = "<<presclpp3<<endl; //9.24968

  //root [9] (Float_t)jet->GetEntries("HLT_HIJet80_v1")/jet->GetEntries("HLT_HIJet80_v1&&HLT_HIJet55_v1")
  //(double)2.34995051108819153e+00
  //ive commented this below - to just check for the pbpb histograms to load. 
  
  jetpp1_v2->Project("hpp1","pt","abs(eta)<2&&jet80&&(chMax/pt)>0.01");
  hpp1->Print("base");
 
  jetpp2_v2->Project("hpp2","pt","abs(eta)<2&&jet60&&!jet80&&(chMax/pt)>0.01");
  hpp2->Print("base");

  jetpp2_v2->Project("hpp3","pt","9.25038"*pp3);
  // 9.25038 was the value. 
  //jetpp2_v2->Project("hpp3","pt","jet40_p"*pp3);
  hpp3->Print("base");
 
  
  hpp1->Scale(1./5300e6);//pp lumi
  hpp2->Scale(1./5300e6);
  hpp3->Scale(1./5300e6);

  hpp1->Scale(1./4);//delta eta
  hpp2->Scale(1./4);
  hpp3->Scale(1./4);

  hppComb->Add(hpp1,1);
  hppComb->Add(hpp2,1);
  hppComb->Add(hpp3,1);
  hppComb->Print("base");

  hppComb = (TH1F*)hppComb->Rebin(nbins_yaxian,"hppComb",boundaries_yaxian);
  hpp3 = (TH1F*)hpp3->Rebin(nbins_yaxian,"hpp3",boundaries_yaxian);
  hpp2 = (TH1F*)hpp2->Rebin(nbins_yaxian,"hpp2",boundaries_yaxian);
  hpp1 = (TH1F*)hpp1->Rebin(nbins_yaxian,"hpp1",boundaries_yaxian);

  divideBinWidth(hppComb);
  divideBinWidth(hpp1);
  divideBinWidth(hpp2);
  divideBinWidth(hpp3);


  //these were for doing it from the forests directly without the proper JEC's 
  //add the centrality cuts: 

  const int nbins_cent = 6;
  Double_t boundaries_cent[nbins_cent+1] = {0,2,4,12,20,28,36};// multiply by 2.5 to get your actual centrality % (old 2011 data) 
  Double_t ncoll[nbins_cent] = { 1660, 1310, 745, 251, 62.8, 10.8 };
  //const int nbins_cent = 1;
  //Double_t boundaries_cent[nbins_cent+1] = {0,40};
  //Double_t ncoll[nbins_cent] = {362.24}; //use taa instead of ncoll. 
  
  TCut pbpb0 = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIZeroBiasPizel_SingleTrack_v1&&chargedMax/jtpt>0.01";//this is just for the MB file. not really used here so far. 

  TCut pbpb1[nbins_cent];
  TCut pbpb2[nbins_cent];
  TCut pbpb3[nbins_cent];

  TH1F *hpbpb1[nbins_cent],*hpbpb2[nbins_cent],*hpbpb3[nbins_cent];
  TH1F *hpbpbComb[nbins_cent];
  //TH1F* htest = new TH1F("htest","",1000,0,1000);

  TH1F* hTurnon80 = new TH1F("hTurnon80","",150,0,150);
  TH1F* hTurnon65 = new TH1F("hTurnon65","",150,0,150);
  TH1F* hTurnon55 = new TH1F("hTurnon55","",150,0,150);

  TH1F* hTriggerMerged = new TH1F("hTriggerMerged","",150,0,150);
  TH1F* htest80 = new TH1F("htest80","",150,0,150);
  TH1F* htest65 = new TH1F("htest65","",150,0,150);
  TH1F* htest55 = new TH1F("htest55","",150,0,150);

  // check the trigger turn on curve from the MB file. 
  TH1F* hMB = new TH1F("hMB","",150,0,150);
  /*
  jetpbpb2->Project("htest80","jtpt","HLT_HIJet80_v1");
  htest80->Print("base");
  jetpbpb1->Project("htest65","jtpt","HLT_HIJet65_v1&&!HLT_HIJet80_v1");
  htest65->Print("base");
  jetpbpb1->Project("htest55","jtpt","HLT_HIJet55_v1_Prescl*(HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1)");
  htest55->Print("base");
  */
  jetpbpb2_old->Project("hTurnon80","jtpt","HLT_HIJet80_v1&&pcollisionEventSelection&&pHBHENoiseFilter");
  jetpbpb1_old->Project("hTurnon65","jtpt","HLT_HIJet65_v1&&pcollisionEventSelection&&pHBHENoiseFilter");
  jetpbpb1_old->Project("hTurnon55","jtpt","HLT_HIJet55_v1&&pcollisionEventSelection&&pHBHENoiseFilter");

  jetpbpb0_old->Project("hMB","jtpt","pcollisionEventSelection&&pHBHENoiseFilter");

  hTurnon80->Print("base");
  hTurnon65->Print("base");
  hTurnon55->Print("base");

  //hTriggerMerged->Add(htest80);
  //hTriggerMerged->Add(htest65);
  //hTriggerMerged->Add(htest55);

  hTurnon80->Divide(hMB);
  hTurnon65->Divide(hMB);
  hTurnon55->Divide(hMB);

  //plot the trigger turn on curves. 
  TCanvas *ctrig = new TCanvas("ctrig","",800,600);
  hTurnon80->SetMarkerStyle(20);
  hTurnon80->SetMarkerColor(kCyan);
  hTurnon65->SetMarkerStyle(20);
  hTurnon65->SetMarkerColor(kRed);
  hTurnon55->SetMarkerStyle(20);
  hTurnon55->SetMarkerColor(kYellow);

  hTurnon80->SetAxisRange(0,150,"X");
  hTurnon80->SetAxisRange(0,1.3,"Y");
  hTurnon80->SetXTitle("offline Jet p_{T} GeV/c");
  hTurnon80->SetYTitle("Trigger Efficiency");
  hTurnon80->GetXaxis()->SetTitleOffset(1.3);
  hTurnon80->Draw();
  hTurnon65->Draw("same");
  hTurnon55->Draw("same");
				

  TLegend *title3 = myLegend(0.68,0.15,0.92,0.45);

  title3->AddEntry(hTurnon55,"HLT_HIJet55","pl");
  title3->AddEntry(hTurnon65,"HLT_HIJet65","pl");
  title3->AddEntry(hTurnon80,"HLT_HIJet80","pl");
  title3->SetTextSize(0.03);
  title3->Draw();
  
  putCMSPrel(0.1,0.92,0.06);
  drawText("PbPb #int dt = 149.382 #mu b^{-1}, #sqrt{s_{NN}}=2.76 TeV",0.5,0.93,16);

  ctrig->SaveAs("RAA_trigger_turnon.pdf","RECREATE");

  //centrality loop for the pbpb files/histograms 
  for(int i = 0;i<nbins_cent;i++){

    cout<<"centrality boundary = "<<boundaries_cent[i]*2.5<<" - "<<boundaries_cent[i+1]*2.5<<endl;

    pbpb1[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%f&&hiBin<%f",boundaries_cent[i],boundaries_cent[i+1]);
    pbpb2[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet65_v1&&!HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%f&&hiBin<%f",boundaries_cent[i],boundaries_cent[i+1]);
    pbpb3[i] = Form("abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2&&HLT_HIJet55_v1&&!HLT_HIJet65_v1&&!HLT_HIJet80_v1&&(chargedMax/jtpt)>0.01&&hiBin>=%f&&hiBin<%f",boundaries_cent[i],boundaries_cent[i+1]);

    hpbpb1[i] = new TH1F(Form("hpbpb1_%d",i),"",1000,0,1000);
    //hpbpb1[i]->Print("base");
    hpbpb2[i] = new TH1F(Form("hpbpb2_%d",i),"",1000,0,1000);
    hpbpb3[i] = new TH1F(Form("hpbpb3_%d",i),"",1000,0,1000);
    hpbpbComb[i] = new TH1F(Form("hpbpbComb_%d",i),"",1000,0,1000);

    jetpbpb2->Project(Form("hpbpb1_%d",i),"jtpt",pbpb1[i]);
    hpbpb1[i]->Print("base");
    //divideBinWidth(hpbpb1);
    
    jetpbpb1->Project(Form("hpbpb2_%d",i),"jtpt",pbpb2[i]);
    hpbpb2[i]->Print("base");
    //divideBinWidth(hpbpb2);
    
    jetpbpb1->Project(Form("hpbpb3_%d",i),"jtpt","HLT_HIJet55_v1_Prescl"*pbpb3[i]);
    //jetpbpb1->Project("hpbpb3","jtpt","2.34995"*pbpb3);
    hpbpb3[i]->Print("base");
    //divideBinWidth(hpbpb3);

    //scale the PbPb histograms before adding them
    //we have to scale them according to the lumi of the Jet80 file. 
    // HLT file  |   Lumi inverse micro barns 
    // HLT_80    |   149.382 
    // HLT_65    |   3.195
    // HLT_55    |   2.734
    // 

    // the files which im using now is only a fraction of events of that:
    // for the PbPb 55 or 65 file its 0.977
    // for the PbPb 80 file its 0.304
    
    hpbpb1[i]->Scale(1./149.382e6/0.304);//respective lumi seen by the trigger all in inverse micro barns 
    hpbpb2[i]->Scale(1./3.195e6/0.977);
    hpbpb3[i]->Scale(1./2.734e6/0.977);
    
    hpbpb1[i]->Scale(1./4);//delta eta
    hpbpb2[i]->Scale(1./4);
    hpbpb3[i]->Scale(1./4);

    hpbpb1[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i]));//centrality bin width. 
    hpbpb2[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i]));
    hpbpb3[i]->Scale(1./0.025/(boundaries_cent[i+1]-boundaries_cent[i]));
    
    //add the histograms  
    hpbpbComb[i]->Add(hpbpb1[i]);
    hpbpbComb[i]->Add(hpbpb2[i]);
    hpbpbComb[i]->Add(hpbpb3[i]);
    hpbpbComb[i]->Print("base");

    hpbpbComb[i] = (TH1F*)hpbpbComb[i]->Rebin(nbins_yaxian,Form("hpbpbComb_%d",i),boundaries_yaxian);
    hpbpb3[i] = (TH1F*)hpbpb3[i]->Rebin(nbins_yaxian,Form("hpbpb3_%d",i),boundaries_yaxian);
    hpbpb2[i] = (TH1F*)hpbpb2[i]->Rebin(nbins_yaxian,Form("hpbpb2_%d",i),boundaries_yaxian);
    hpbpb1[i] = (TH1F*)hpbpb1[i]->Rebin(nbins_yaxian,Form("hpbpb1_%d",i),boundaries_yaxian);

    divideBinWidth(hpbpbComb[i]);
    divideBinWidth(hpbpb3[i]);
    divideBinWidth(hpbpb2[i]);
    divideBinWidth(hpbpb1[i]);
    
  }
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  c1->Divide(3,2);

  for(int i = 0;i<nbins_cent;i++){
    c1->cd(i+1);
    c1->cd(i+1)->SetLogy();
    hpbpbComb[i]->SetMarkerStyle(20);
    //hpbpbComb[i]->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
    hpbpbComb[i]->SetYTitle("#frac{d^{2} #sigma}{d p_{T} d#eta} #mu barns");
    hpbpbComb[i]->SetXTitle("Jet p_{T} GeV/c");
    
    //TF1 *fPowerLaw = new TF1("fPowerLaw","[0]*pow(x+[1],[2])");
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    //hpbpbComb->Fit("fPowerLaw","","",30,300);
    hpbpbComb[i]->SetAxisRange(50,300,"X");
    //hpbpbComb[i]->SetAxisRange(1e-4,1e-12,"Y");
    hpbpbComb[i]->Draw();
    hpbpb3[i]->SetMarkerStyle(20);
    hpbpb3[i]->SetMarkerColor(kRed);
    hpbpb3[i]->Draw("same");
    hpbpb2[i]->SetMarkerStyle(20);
    hpbpb2[i]->SetMarkerColor(kBlue);
    hpbpb2[i]->Draw("same");
    hpbpb1[i]->SetMarkerStyle(20);
    hpbpb1[i]->SetMarkerColor(kGreen);
    hpbpb1[i]->Draw("same");
    hpbpbComb[i]->Draw("same");

    drawText(Form("%2.0f-%2.0f",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5),0.1,0.9,20);

    if(i == nbins_cent-1){
       
      TLegend *title = myLegend(0.25,0.65,0.55,0.8);
      title->AddEntry(hpbpbComb[nbins_cent-1],"PbPb 2.76 TeV Data Merged","pl");
      title->AddEntry(hpbpb3[nbins_cent-1],"w_{3} * (HLT_55 && !HLT_65 && !HLT_80)","pl");
      title->AddEntry(hpbpb2[nbins_cent-1],"HLT_65 && !HLT_80","pl");
      title->AddEntry(hpbpb1[nbins_cent-1],"HLT_80","pl");
      title->SetTextSize(0.03);
      title->Draw();
      //drawText("PbPb data",0.3,0.65,20);  
      drawText(Form("Anti-k_{T} %s PF Jets R = 0.%d, |#eta|<2, |vz|<15",algo,radius),0.2,0.8,14);
      
    }
    
  }
  
 
  c1->SaveAs(Form("pbpb_2011_vs_pt_combined_nbins_cent_%d.pdf",nbins_cent),"RECREATE");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();
  hppComb->SetMarkerStyle(20);
  //hppComb->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
  hppComb->SetYTitle("#frac{d^2 #sigma}{d p_{T} d#eta} #mu barns");
  hppComb->SetXTitle("Jet p_{T} GeV/c");

  //TF1 *fPowerLaw2 = new TF1("fPowerLaw2","[0]*pow(x+[1],[2])");
  //hppComb->Fit("fPowerLaw2","","",30,300);
  //hppComb->Fit("fPowerLaw2","","",30,300);
  //hppComb->Fit("fPowerLaw2","","",30,300);
  //hppComb->Fit("fPowerLaw2","","",30,300);
  //hppComb->Fit("fPowerLaw2","","",30,300);
  //hppComb->Fit("fPowerLaw2","","",30,300);
  //hppComb->Fit("fPowerLaw2","","",30,300);
  hppComb->SetAxisRange(20,300,"X");
  //hppComb->SetAxisRange(1e-4,1e-12,"Y");
  hppComb->Draw();
  hpp3->SetMarkerStyle(20);
  hpp3->SetMarkerColor(kRed);
  hpp3->Draw("same");
  hpp2->SetMarkerStyle(20);
  hpp2->SetMarkerColor(kBlue);
  hpp2->Draw("same");
  hpp1->SetMarkerStyle(20);
  hpp1->SetMarkerColor(kGreen);
  hpp1->Draw("same");
  hppComb->Draw("same");
  TLegend *title2 = myLegend(0.25,0.65,0.55,0.8);
  title2->AddEntry(hppComb,"PP 2.76 TeV Data  Merged with Doga's residual corrections","pl");
  title2->AddEntry(hpp3,"w_{3} * (HLT_40 && !HLT_60 && !HLT_80)","pl");
  title2->AddEntry(hpp2,"HLT_60 && !HLT_80","pl");
  title2->AddEntry(hpp1,"HLT_80","pl");
  title2->SetTextSize(0.03);
  title2->Draw();
  //drawText("PP data",0.3,0.65,20);  
  drawText(Form("Anti-k_{T} PF Jets R = 0.%d, |#eta|<2, |vz|<15",radius),0.35,0.56,20);
  c2->SaveAs("pp_2011_pt_combined.pdf","RECREATE");

  // calculate the measured RAA here just as a simple calculation 
  hppComb->Scale(1./64); //sigma pp
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->Divide(3,2);

  TH1F* hRAA[nbins_cent];

  for(int i = 0;i<nbins_cent;i++){

    hpbpbComb[i]->Scale(1./ncoll[i]); // ncoll
    hpbpbComb[i]->Scale(1./7.65);//remember what this is. maybe sigma inelastic

    hRAA[i] = (TH1F*)hpbpbComb[i]->Clone("hRAA");
    hRAA[i]->Print("base");
    hRAA[i]->Divide(hppComb);

    c3->cd(i+1);
    hRAA[i]->SetYTitle("R_{AA}");
    hRAA[i]->SetXTitle("p_{T} GeV/c");
    hRAA[i]->SetMarkerStyle(20);
    hRAA[i]->SetMarkerColor(kBlack);
    hRAA[i]->SetAxisRange(0,2,"Y");
    hRAA[i]->Draw();
    drawText(Form("%2.0f-%2.0f",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5),0.1,0.9,20);

    if(i == 5)drawText(Form("R = 0.%d, |#eta|<2, |vz|<15",radius),0.35,0.56,14);
    //drawText("|#eta|<2, |vz|<15 0-100%",0.35,0.76,20);
  }

  c3->SaveAs(Form("RAA_March2014_voronoi_nbins_cent_%d.pdf",nbins_cent),"RECREATE");

  //plot the statistical uncertainty here
  //statistical error/meanvalue as a function of pt for the combined spectra. 
  /*
  TCanvas *c2 = new TCanvas("c2","",800,600);
  //TH1F* hPbPb_Uncert = (TH1F*)hpbpbComb->Clone("hPbPb_Uncert");
 
  TH1F* hPbPb_Uncert = new TH1F("hPbPb_Uncert","",30,0,300);
 
  for(int i = 1;i<=hpbpbComb->GetNbinsX();i++){
    
    double val = hpbpbComb->GetBinContent(i);
    double valErr = hpbpbComb->GetBinError(i);
    double uncert = (double)valErr/val;
    cout<<"uncert = "<<uncert<<endl;
    hPbPb_Uncert->SetBinContent(i,uncert);
    hPbPb_Uncert->SetBinError(i,0);
  }

  hPbPb_Uncert->SetYTitle("uncertainty");
  hPbPb_Uncert->SetXTitle("p_{T} GeV/c");
  hPbPb_Uncert->Draw();
  drawText("PbPb 2011, 55,65 scaled",0.3,0.65,20);  
  drawText("Anti-k_{T} PU PF Jets R = 0.3, |#eta|<2, |vz|<15",0.3,0.56,20);
  c2->SaveAs("pbpb_2013_hlt_merge_scaled_uncert.gif","RECREATE");
  */
  /*
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();
  TH1F* hPPComb = (TH1F*)hppComb->Clone("hPPComb");
  //TH1F* hPPComb_bins = rebin_yaxian(hppComb,"hPPComb_bins");
  
  hPPComb->Scale(1./3.083e11);
  hPPComb->Print("base");
  hPPComb->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
  hPPComb->SetXTitle("Jet p_{T} GeV/c");
  //hPPComb_bins->Scale(1./3.083e11);
  //hPPComb_bins->Print("base");
  //hPPComb_bins->Scale(1./4);
  //divideBinWidth(hPPComb_bins);

  hPPComb->Scale(1./4);
  divideBinWidth(hPPComb);

  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]/pow(x,[1])");
  hPPComb->Fit("fPowerLaw","","",25,500);
  hPPComb->Fit("fPowerLaw","","",25,500);
  hPPComb->Fit("fPowerLaw","","",25,500);
  hPPComb->Fit("fPowerLaw","","",25,500);
  hPPComb->Fit("fPowerLaw","","",25,500);
  hPPComb->SetMarkerColor(kBlue);
  hPPComb->SetMarkerStyle(26);
  hPPComb->SetTitle("PP2013 ak3PF");
  hPPComb->Draw();

  //hPPComb_bins->SetMarkerColor(kRed);
  //hPPComb_bins->SetMarkerStyle(23);
  //hPPComb_bins->Draw("same");

  c2->SaveAs("pp_2013_ak3_pt_evt_frac_merged.gif","RECREATE");

  TCanvas *c5 = new TCanvas("c5","",800,600);
  TH1F* hppFunc = (TH1F*)functionHist(fPowerLaw,hppComb,"Fit Function p_{T} spectra PP 2013 merged");
  TH1F* hPPRatio = (TH1F*)hppComb->Clone("hPPRatio");
  hPPRatio->Divide(hppFunc);
  hPPRatio->SetTitle("Spectra to Fit Ratio");
  hPPRatio->SetXTitle("Jet p_{T} GeV/c");
  hPPRatio->SetYTitle("Measured data/Fit");
  hPPRatio->SetMarkerStyle(8);
  hPPRatio->SetMarkerColor(4);
  hPPRatio->Draw();
  c5->SaveAs("pp_2013_merged_spectra_fit_comp.gif","RECREATE");
  

  
  TFile *fpbpbunfo = TFile::Open("result-2013-akPu3PF-cent-6-isFineBin-0/pbpb_pp_merged_chmx_pt_Unfo_2013_akPu3PF_cent_6_isFineBin_0.root");
  TH1F* hppUnfo = (TH1F*)fpbpbunfo->Get("Unfolded_cent6");
  TH1F* hPPGen = (TH1F*)fpbpbunfo->Get("hGen_cent6");
  hppUnfo->Print("base");
  hPPGen->Print("base");

  hPPGen->Scale(1./4);
  divideBinWidth(hPPGen);

  hppUnfo->Scale(1./3.083e11);
  hppUnfo->Scale(1./4);
  divideBinWidth(hppUnfo);

  hppUnfo->Divide(hPPGen);

  TCanvas *c6 = new TCanvas("c6","",800,600);
  hppUnfo->SetMarkerStyle(21);
  hppUnfo->SetMarkerColor(kRed);
  hPPGen->SetMarkerStyle(21);
  hPPGen->SetMarkerColor(kBlue);
  hppUnfo->Draw();
  // hPPGen->Draw("same");
  c6->SaveAs("pp_2760GeV_unfold_vs_mc.gif","RECREATE");
  */
  
  /*
  TCanvas *c7 = new TCanvas("c7","",800,600);
  c7->SetLogy();
  hPPComb->Draw();
  hPPComb->SetYTitle("");
  hPPComb->SetXTitle("p_{T} GeV/c");
  //hppComb->SetTitle("PP 2013 2.76 TeV ak4PF measured vs unfolded");
  hPPComb->SetMarkerStyle(23);
  hPPComb->SetMarkerColor(kBlue);
  //hppUnfo->SetAxisRange(10,500,"X");
  //hppUnfo->SetMarkerStyle(24);
  //hppUnfo->SetMarkerColor(kRed);
  //hppUnfo->Draw("same");
  
  
  
  TLegend *title5 = myLegend(0.54,0.65,0.85,0.9);
  title5->AddEntry(hppComb,"Measured","pl");
  title5->AddEntry(hppUnfo,"Bayesian iter = 4","pl");
  title5->SetTextSize(0.06);
  title5->Draw();
  gStyle->SetOptStat(0);
  c7->SaveAs("PP2013_measured_vs_unfolded.gif","RECREATE");
  
  //TCanvas 
  */

  //Create output file and save them. 
  TFile f(Form("merge_pbpb_ak%d_%s_HLT_V2_nbins_cent_%d.root",radius,algo,nbins_cent),"RECREATE");
  for(int i = 0;i<nbins_cent;i++){
    hpbpb1[i]->Write();
    hpbpb2[i]->Write();
    hpbpb3[i]->Write();
    hpbpbComb[i]->Write();
    hRAA[i]->Write();
  }
  //hPPComb_bins->Write();
  hpp1->Write();
  hpp2->Write();
  hpp3->Write();
  hppComb->Write();
  //hPPComb->Write();
  //hPbPb_Uncert->Write();
  //hPPComb->Write();
  //hPPGen->Write();
  f.Close();
  
  


}
