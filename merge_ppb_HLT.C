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

static const int nbins_rec = 50;
static const double boundaries_rec[nbins_rec+1] = {
  0,10,20,30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170,
	180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310,
	320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450,
	460, 470, 480, 490, 500
	};


// rebin the spectra
TH1F *rebin2(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_rec,boundaries_rec);
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

  TH1::SetDefaultSumw2();
  
  TFile *fppb0_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppbJetMB_v2.root");
  TFile *fppb1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppbJet40_v2.root");
  TFile *fppb2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_JEC_applied_ppbJet80_v2.root");

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

  //number convention is going to change here to better understand the merging. 
  // 0 - MB
  // 1 - 20
  // 2 - 40
  // 3 - 60
  // 4 - 80
  // 5 - 100
  // 6 - 120

  

  //these histograms are only for merging. 
  TH1F *hppb0 = new TH1F("hppb0","",50,0,500);
  TH1F *hppb1 = new TH1F("hppb1","",50,0,500);
  TH1F *hppb2 = new TH1F("hppb2","",50,0,500);
  TH1F *hppb3 = new TH1F("hppb3","",50,0,500);
  TH1F *hppb4 = new TH1F("hppb4","",50,0,500);
  TH1F *hppb5 = new TH1F("hppb5","",50,0,500);
  TH1F *hppb6 = new TH1F("hppb6","",50,0,500);
  TH1F *hppbComb = new TH1F("hppbComb","",50,0,500);
  
  //define histograms which are used for trigger turn on curves. 
  TH1F* hHLT_120 = new TH1F("hHLT_120","",50,0,500);
  TH1F* hHLT_100 = new TH1F("hHLT_100","",50,0,500);
  TH1F* hHLT_80 = new TH1F("hHLT_80","",50,0,500);
  TH1F* hHLT_60 = new TH1F("hHLT_60","",50,0,500);
  TH1F* hHLT_40 = new TH1F("hHLT_40","",50,0,500);
  TH1F* hHLT_20 = new TH1F("hHLT_20","",50,0,500);
  TH1F* hMB = new TH1F("hMB","",50,0,500);
  TH1F* hMB_0 = new TH1F("hMB_0","",50,0,500);
  TH1F* hMB_1 = new TH1F("hMB_1","",50,0,500);
  TH1F* hMB_2 = new TH1F("hMB_2","",50,0,500);
  TH1F* hBase_100 = new TH1F("hBase_100","",50,0,500);
  TH1F* hBase_40 = new TH1F("hBase_40","",50,0,500);

  /*
  Float_t MB_count = (Float_t)jetppb0_v2->GetEntries();
  cout<<"MB_count = "<<MB_count<<endl; //2.89687e07
    
  
  Float_t presclppb0 = (Float_t)jetppb2_v2->GetEntries("jet100")/jetppb2_v2->GetEntries("jetMB&&jet100");
  Float_t presclppb1 = (Float_t)jetppb2_v2->GetEntries("jet100")/jetppb2_v2->GetEntries("jet20&&jet100");
  Float_t presclppb2 = (Float_t)jetppb2_v2->GetEntries("jet100")/jetppb2_v2->GetEntries("jet40&&jet100");
  Float_t presclppb3 = (Float_t)jetppb2_v2->GetEntries("jet100")/jetppb2_v2->GetEntries("jet60&&jet100");
  Float_t presclppb4 = (Float_t)jetppb2_v2->GetEntries("jet100")/jetppb2_v2->GetEntries("jet80&&jet100");
  //Float_t presclppb5 = (Float_t)jetppb1_v2->GetEntries("jet100")/jetppb1_v2->GetEntries("jet100&&jet100");

  Float_t Presclppb0 = (Float_t)jetppb2_v2->GetEntries("jetMB_p")/jetppb2_v2->GetEntries("jetMB");
  Float_t Presclppb1 = (Float_t)jetppb2_v2->GetEntries("jet20_p")/jetppb2_v2->GetEntries("jet20");
  Float_t Presclppb2 = (Float_t)jetppb2_v2->GetEntries("jet40_p")/jetppb2_v2->GetEntries("jet40");
  Float_t Presclppb3 = (Float_t)jetppb2_v2->GetEntries("jet60_p")/jetppb2_v2->GetEntries("jet60");
  Float_t Presclppb4 = (Float_t)jetppb2_v2->GetEntries("jet80_p")/jetppb2_v2->GetEntries("jet80");
  Float_t Presclppb5 = (Float_t)jetppb2_v2->GetEntries("jet100_p")/jetppb2_v2->GetEntries("jet100");
  Float_t Presclppb6 = (Float_t)jetppb2_v2->GetEntries("jet120_p")/jetppb2_v2->GetEntries("jet120");

  cout<<"ppb prescl0 = "<<presclppb0<<endl;//2202.28
  cout<<"ppb prescl1 = "<<presclppb1<<endl;//1163.58
  cout<<"ppb prescl2 = "<<presclppb2<<endl;//40.9373
  cout<<"ppb prescl3 = "<<presclppb3<<endl;//4.46062
  cout<<"ppb prescl4 = "<<presclppb4<<endl;//1.41429

  cout<<"ppb Prescl0 = "<<Presclppb0<<endl;//
  cout<<"ppb Prescl1 = "<<Presclppb1<<endl;//
  cout<<"ppb Prescl2 = "<<Presclppb2<<endl;//
  cout<<"ppb Prescl3 = "<<Presclppb3<<endl;//
  cout<<"ppb Prescl4 = "<<Presclppb4<<endl;//
  cout<<"ppb Prescl5 = "<<Presclppb5<<endl;//
  cout<<"ppb Prescl6 = "<<Presclppb6<<endl;//
  */
 

  //read the required histograms including the prescale factor 
  TCut ppb0 = "abs(eta_CM)<1&&jetMB&&!jet80&&!jet100&&chMax/pt>0.01";
  TCut ppb1 = "abs(eta_CM)<1&&jet20&&!jet80&&!jet100&&chMax/pt>0.01";
  TCut ppb2 = "abs(eta_CM)<1&&jet40&&!jet80&&!jet100&&chMax/pt>0.01";
  TCut ppb3 = "abs(eta_CM)<1&&jet60&&!jet80&&!jet100&&chMax/pt>0.01";
  TCut ppb4 = "abs(eta_CM)<1&&jet80&&!jet100&&chMax/pt>0.01";//try with or without 60 first and then try without 60
  TCut ppb5 = "abs(eta_CM)<1&&jet100&&chMax/pt>0.01";
  

  jetppb0_v2->Project("hppb0","pt","2202.20"*ppb0);
  hppb0->Print("base");

  //jetppb1_v2->Project("hppb1","pt","1163.58"*ppb1);
  //hppb1->Print("base");

  //jetppb1_v2->Project("hppb2","pt","40.9373"*ppb2);
  //hppb2->Print("base");

  //jetppb1_v2->Project("hppb3","pt","4.46062"*ppb3);
  //hppb3->Print("base");

  jetppb2_v2->Project("hppb4","pt",ppb4);
  hppb4->Print("base");
  //divideBinWidth(hppb4);
  //hppb4->Scale(1./5.25585e10);
  //hppb4->Scale(1./4);

  jetppb2_v2->Project("hppb5","pt",ppb5);
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
  /*
  TCut turnon_MB = "jetMB";
  TCut turnon_20 = "jet20";
  TCut turnon_40 = "jet40";
  TCut turnon_60 = "jet60";
  TCut turnon_80 = "jet80";
  TCut turnon_100 = "jet100";
  //TCut turnon_120 = "jet120";
  
  jetppb0_v2->Project("hMB","pt","2202.2"*turnon_MB);
  hMB->Print("base");
  jetppb0_v2->Project("hMB_2","pt","abs(eta_CM)<1&&chMax/pt>0.01");
  hMB_2->Print("base");
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
  jetppb1_v2->Project("hBase_40","pt");
  hBase_40->Print("base");
  //jetppb2_v2->Project("h120","pt",""*turnon_120);
  */

  jetppb0_v2->Project("hMB","ptlead","jetMB");
  jetppb1_v2->Project("hMB_1","ptlead","jetMB");
  jetppb2_v2->Project("hMB_2","ptlead","jetMB");
  jetppb1_v2->Project("hHLT_20","ptlead","jetMB&&jet20");
  jetppb1_v2->Project("hHLT_40","ptlead","jetMB&&jet40");
  jetppb1_v2->Project("hHLT_60","ptlead","jetMB&&jet60");
  jetppb2_v2->Project("hHLT_80","ptlead","jetMB&&jet80");
  
  
  jetppb2_v2->Project("hHLT_100","ptlead","1"*turnon_100);
  hHLT_100->Print("base");
  jetppb2_v2->Project("hBase_100","ptlead");
  hBase_100->Print("base");

  TH1F* hTurnon_20 = (TH1F*)hHLT_20->Clone("hTurnon_20");
  TH1F* hTurnon_40 = (TH1F*)hHLT_40->Clone("hTurnon_40");
  TH1F* hTurnon_60 = (TH1F*)hHLT_60->Clone("hTurnon_60");
  TH1F* hTurnon_80 = (TH1F*)hHLT_80->Clone("hTurnon_80");
  TH1F* hTurnon_100 = (TH1F*)hHLT_100->Clone("hTurnon_100");

  hTurnon_20->Divide(hMB_1);
  hTurnon_40->Divide(hMB_1);
  hTurnon_60->Divide(hMB_1);
  hTurnon_80->Divide(hMB_2);
  hTurnon_100->Divide(hBase_100);
  
  TFile fout("ppb_turnon_HLT.root","RECREATE");
  hTurnon_20->Write();
  hTurnon_40->Write();
  hTurnon_60->Write();
  hTurnon_80->Write();
  hTurnon_100->Write();
  hMB->Write();
  hHLT_20->Write();
  hHLT_40->Write();
  hHLT_60->Write();
  hHLT_80->Write();
  hHLT_100->Write();
  fout.Close();

  TCanvas *c1 = new TCanvas("c1","",800,600);
  hTurnon_20->SetAxisRange(0,2,"Y");
  hTurnon_20->SetAxisRange(0,140,"X");
  hTurnon_20->SetTitle("HLT Turn on curves");
  hTurnon_20->SetXTitle("p_{T} GeV/c");
  hTurnon_20->SetMarkerStyle(24);
  hTurnon_20->SetMarkerColor(2);
  hTurnon_40->SetMarkerStyle(25);
  hTurnon_40->SetMarkerColor(3);
  hTurnon_60->SetMarkerStyle(26);
  hTurnon_60->SetMarkerColor(4);
  hTurnon_80->SetMarkerStyle(27);
  hTurnon_80->SetMarkerColor(1);
  hTurnon_100->SetMarkerStyle(28);
  hTurnon_100->SetMarkerColor(6);
  hTurnon_20->Draw();
  hTurnon_40->Draw("same");
  hTurnon_60->Draw("same");
  hTurnon_80->Draw("same");
  hTurnon_100->Draw("same");
  
  TLegend *legend = myLegend(0.34,0.65,0.65,0.9);
  legend->AddEntry(hTurnon_20,"HLT_20","pl");
  legend->AddEntry(hTurnon_40,"HLT_40","pl");
  legend->AddEntry(hTurnon_60,"HLT_60","pl");
  legend->AddEntry(hTurnon_80,"HLT_80","pl");
  legend->AddEntry(hTurnon_100,"HLT_100","pl");
  legend->SetTextSize(0.04);
  legend->Draw();
  
  drawText("PPb 2013 Data Prompt Reco",0.15,0.8,20);
  drawText("Anti-k_{T} PU PF Jets R = 0.3,|vz|<15",0.15,0.7,20);

  c1->SaveAs("ppb_2013_HLT_turn_on_curve.jpg","RECREATE");

  //
  TCanvas* c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();
  TH1F* hPPbComb = rebin2(hppbComb,"hPPbComb");
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
  //TH1F* hPPb_Unfo = rebin2(hPPb_Unfo_test,"hPPb_Unfo");
  hPPb_Unfo->SetAxisRange(20,500,"X");
  TH1F* hPPb_Unfo_2 = (TH1F*)hPPb_Unfo->Clone("hPPb_Unfo_2");
  TH1F* hPPb_Unfo_3 = (TH1F*)hPPb_Unfo->Clone("hPPb_Unfo_3");
  TH1F* hPPb_Unfo_BinByBin = (TH1F*)fPPb_Unfo->Get("UnfoldedBinByBin_cent0");
  //TH1F* hPPb_Unfo_BinByBin = rebin2(hPPb_Unfo_BinByBin_test,"hPPb_Unfo_BinByBin");
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
  c3->SaveAs("ppb_merged_MB_eta_cm_1_lowest_measured_vs_unfolded.jpg","RECREATE");

  //create the spectra with event fraction. 
  hPPbComb_2->Scale(1./6.53104e10);
  hPPbComb_2->Scale(1./2);
  hPPb_Unfo_2->Scale(1./6.53104e10);
  hPPb_Unfo_2->Scale(1./2);
  hMB_2->Scale(1./2.89687e07);
  hMB_2->Scale(1./2);
  divideBinWidth(hPPbComb_2);
  divideBinWidth(hPPb_Unfo_2);

  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetLogy();
  hPPbComb_2->SetAxisRange(30,500,"X");
  hPPbComb_2->SetYTitle("#frac{dN}{N_{MB} dp_{T} d#eta}");
  hPPbComb_2->SetXTitle("Jet p_{T} GeV/c");
  hPPbComb_2->SetMarkerStyle(7);
  hPPbComb_2->SetMarkerColor(kBlack);
  hPPbComb_2->Draw();
  hMB_2->SetMarkerStyle(26);
  hMB_2->SetMarkerColor(kRed);
  hMB_2->Draw("same");
  hPPb_Unfo_2->SetMarkerStyle(24);
  hPPb_Unfo_2->SetMarkerColor(kBlack);
  hPPb_Unfo_2->Draw("same");

  TLegend *title2 = myLegend(0.34,0.65,0.65,0.9);
  title2->AddEntry(hPPbComb_2,"PPb(merged) 2013 Measured","pl");
  title2->AddEntry(hMB_2,"PPb MB","pl");
  title2->AddEntry(hPPb_Unfo_2,"PPb(merged) 2013 Bayesian","pl");
  title2->SetTextSize(0.04);
  title2->Draw();
  //c4->SaveAs("ppb_merged_60_lowest_evtfrac_measured_vs_unfolded.gif","RECREATE");  
  c4->SaveAs("ppb_merged_MB_eta_CM_1_lowest_evtfrac_measured_vs_unfolded.jpg","RECREATE");  

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
  drawText("|eta_CM|<1, refpt>=15, chMax/jtpt>0.01, akPu3PF",0.3,0.75,20);
  c6->SaveAs("ppb_2013_mc_reco_gen_eta_CM_1.jpg","RECREATE");


  //TFile f("merge_ppb_60_lowest_HLT_V2.root","RECREATE");
  TFile f("merge_ppb_MB_eta_CM_1_lowest_HLT_V2.root","RECREATE");
  hppb1->Write();
  hppb2->Write();
  //hppb3->Write();
  hppb4->Write();
  hppbComb->Write();
  hPPbComb->Write();
  hPPb_Unfo->Write();
  hPPbComb_2->Write();
  hPPb_Unfo_2->Write();
  hPPbComb_3->Write();
  hPPb_Unfo_3->Write();
  f.Close();
  
  



}
