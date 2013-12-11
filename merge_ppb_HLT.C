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

void merge_ppb_HLT(){

  TH1::SetDefaultSumw2();

  TFile *fppb1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_ppbJet80_v2.root");
  TFile *fppb2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/2013/data/ntuple_2013_ppbJet40_v2.root");

  TTree *jetppb1_v2 = (TTree*)fppb1_v2->Get("jetR3");
  TTree *jetppb2_v2 = (TTree*)fppb2_v2->Get("jetR3");

  TTree *evtppb1_v2 = (TTree*)fppb1_v2->Get("evt");
  TTree *evtppb2_v2 = (TTree*)fppb2_v2->Get("evt");

  jetppb1_v2->AddFriend(evtppb1_v2);
  jetppb2_v2->AddFriend(evtppb2_v2);

  TCut ppb3 = "abs(eta)<2&&jet60&&!jet80&&!jet100&&chMax/pt>0.01";
  //TCut ppb4 = "abs(eta_CM)<1&&jet40&&!jet60&&!jet80&&!jet100&&chMax/pt>0.01";

  TH1F *hppb1 = new TH1F("hppb1","",50,0,500);
  TH1F *hppb2 = new TH1F("hppb2","",50,0,500);
  TH1F *hppb3 = new TH1F("hppb3","",50,0,500);
  //TH1F *hppb4 = new TH1F("hppb4","",50,0,500);
  TH1F *hppbComb = new TH1F("hppbComb","",50,0,500);
  
  Float_t presclppb1 = (Float_t)jetppb1_v2->GetEntries("jet100")/jetppb1_v2->GetEntries("jet80&&jet100");
  Float_t presclppb2 = (Float_t)jetppb1_v2->GetEntries("jet100")/jetppb1_v2->GetEntries("jet60&&jet100");
  Float_t presclppb3 = (Float_t)jetppb1_v2->GetEntries("jet100")/jetppb1_v2->GetEntries("jet40&&jet100");

  cout<<"ppb prescl1 = "<<presclppb1<<endl;//0.400531
  cout<<"ppb prescl2 = "<<presclppb2<<endl;//1.5932
  cout<<"ppb prescl3 = "<<presclppb3<<endl;//15.2208
  
  //these values were wrong since we didnt take the HLT_100 which is not prescaled. 
  //35.5661 - 40 is the lowest
  // 3.71943 for 60 as the lowest 
  
  jetppb1_v2->Project("hppb1","pt","abs(eta)<2&&jet100&&chMax/pt>0.01&&chMax/pt>0.01");
  hppb1->Print("base");
  //divideBinWidth(hppb1);
  //hppb1->Scale(1./6.01453e10); //total no of min bias events for 80
  //hppb1->Scale(1./4);

  jetppb2_v2->Project("hppb2","pt","abs(eta)<2&&jet80&&!jet100&&chMax/pt>0.01&&chMax/pt>0.01");
  hppb2->Print("base");
  //divideBinWidth(hppb2);
  //hppb2->Scale(1./5.25585e10);
  //hppb2->Scale(1./4);

  //jetppb2_v2->Project("hppb3","pt","1.58606"*ppb3);
  jetppb2_v2->Project("hppb3","pt","4.46062"*ppb3);//now its telling me its 4.46101 after using jet60&&jet100
  hppb3->Print("base");
  //divideBinWidth(hppb3);
  //hppb3->Scale(1./4.33834e10);
  //hppb3->Scale(1./4);

  //jetppb2_v2->Project("hppb4","pt","15.2208"*ppb4);
  //hppb4->Print("base");
  //divideBinWidth(hppb4);
  //hppb4->Scale(1./4.33834e10);
  //hppb4->Scale(1./4);

  hppbComb->Add(hppb1,1);
  hppbComb->Add(hppb2,1);
  hppbComb->Add(hppb3,1);
  //hppbComb->Add(hppb4,1);
  hppbComb->Print("base");

 
  TCanvas* c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();
  TH1F* hPPbComb = (TH1F*)hppbComb->Clone("hPPbComb");
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
  //hppb4->SetMarkerStyle(27);
  //hppb4->SetMarkerColor(28);
  //hppb4->Draw("same");
  hppb3->SetMarkerStyle(25);
  hppb3->SetMarkerColor(kBlue);
  hppb3->Draw("same");
  hppb2->SetMarkerStyle(26);
  hppb2->SetMarkerColor(kRed);
  hppb2->Draw("same");
  hppb1->SetMarkerStyle(24);
  hppb1->SetMarkerColor(kGreen);
  hppb1->Draw("same");

  TLegend *title = myLegend(0.54,0.65,0.85,0.9);
  title->AddEntry(hppbComb,"PPb Merged","pl");
  //title->AddEntry(hppb4,"PPb HLT_40","pl");
  title->AddEntry(hppb3,"PPb HLT_60","pl");
  title->AddEntry(hppb2,"PPb HLT_80","pl");
  title->AddEntry(hppb1,"PPb HLT_100","pl");
  title->SetTextSize(0.04);
  title->Draw();
  drawText("PPb 2013 Data PromptReco",0.3,0.65,20);  
  drawText("Anti-k_{T} PU PF Jets R = 0.3, |#eta|<2, |vz|<15",0.3,0.56,20);

  c2->SaveAs("ppb_2013_pt_merged_60_lowest.gif","RECREATE");

  //comparison with the unfolded data: 
  TFile* fPPb_Unfo = TFile::Open("result-2013-ppb-akPu3PF-cent-1/ppb_merge_60_lowest_pp_mc_Unfo_akPu3PF_cent_1.root");
  TH1F* hPPb_Unfo = (TH1F*)fPPb_Unfo->Get("Unfolded_cent0");
  hPPb_Unfo->SetAxisRange(30,500,"X");
  TH1F* hPPb_Unfo_2 = (TH1F*)hPPb_Unfo->Clone("hPPb_Unfo_2");
  TH1F* hPPb_Unfo_3 = (TH1F*)hPPb_Unfo->Clone("hPPb_Unfo_3");
  TH1F* hPPb_Unfo_BinByBin = (TH1F*)fPPb_Unfo->Get("UnfoldedBinByBin_cent0");
  hPPb_Unfo_BinByBin->SetAxisRange(30,500,"X");
  hPPb_Unfo_BinByBin->Print("base");
  hPPb_Unfo->Print("base");
  hPPb_Unfo_2->Print("base");
  hPPb_Unfo_3->Print("base");
  divideBinWidth(hPPb_Unfo);
  divideBinWidth(hPPb_Unfo_BinByBin);
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetLogy();
  hPPbComb->SetAxisRange(30,500,"X");
  hPPbComb->SetYTitle("counts");
  hPPbComb->SetXTitle("Jet p_{T} GeV/c");
  hPPbComb->SetTitle("pPb 2013 Merged");
  hPPbComb->SetMarkerStyle(7);
  hPPbComb->SetMarkerColor(kBlack);
  /*
  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]/pow(x,[1])");
  hPPbComb->Fit("fPowerLaw","","",50,500);
  hPPbComb->Fit("fPowerLaw","","",50,500);
  hPPbComb->Fit("fPowerLaw","","",50,500);
  hPPbComb->Fit("fPowerLaw","","",50,500);
  hPPbComb->Fit("fPowerLaw","","",50,500);
  */
  hPPbComb->Draw();
  hPPb_Unfo->SetMarkerStyle(24);
  hPPb_Unfo->SetMarkerColor(kBlack);
  hPPb_Unfo->Draw("same");
  hPPb_Unfo_BinByBin->SetMarkerStyle(26);
  hPPb_Unfo_BinByBin->SetMarkerColor(kRed);
  hPPb_Unfo_BinByBin->Draw("same");

  TLegend *title1 = myLegend(0.54,0.65,0.85,0.9);
  title1->AddEntry(hPPbComb,"Measured","pl");
  title1->AddEntry(hPPb_Unfo,"Bayesian","pl");
  title1->AddEntry(hPPb_Unfo_BinByBin,"Bin by Bin","pl");
  title1->SetTextSize(0.04);
  title1->Draw();
  c3->SaveAs("ppb_merged_60_lowest_measured_vs_unfolded.gif","RECREATE");

  //create the spectra with event fraction. 
  hPPbComb_2->Scale(1./6.53104e10);
  hPPbComb_2->Scale(1./4);
  hPPb_Unfo_2->Scale(1./6.53104e10);
  hPPb_Unfo_2->Scale(1./4);
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
  hPPb_Unfo_2->SetMarkerStyle(24);
  hPPb_Unfo_2->SetMarkerColor(kBlack);
  hPPb_Unfo_2->Draw("same");

  TLegend *title2 = myLegend(0.54,0.65,0.85,0.9);
  title2->AddEntry(hPPbComb_2,"PPb(merged) 2013 Measured","pl");
  title2->AddEntry(hPPb_Unfo_2,"PPb(merged) 2013 Bayesian","pl");
  title2->SetTextSize(0.04);
  title2->Draw();
  c4->SaveAs("ppb_merged_60_lowest_evtfrac_measured_vs_unfolded.gif","RECREATE");  

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

  TLegend *title3 = myLegend(0.54,0.65,0.85,0.9);
  title3->AddEntry(hPPbComb_2,"PPb(merged) 2013 Measured","pl");
  title3->AddEntry(hPPb_Unfo_2,"PPb(merged) 2013 Bayesian","pl");
  title3->SetTextSize(0.04);
  title3->Draw();
  c5->SaveAs("ppb_merged_60_lowest_crosssection_measured_vs_unfolded.gif","RECREATE");  



  TFile f("merge_ppb_60_lowest_HLT_V2.root","RECREATE");
  hppb1->Write();
  hppb2->Write();
  hppb3->Write();
  hppbComb->Write();
  hPPbComb->Write();
  hPPb_Unfo->Write();
  hPPbComb_2->Write();
  hPPb_Unfo_2->Write();
  hPPbComb_3->Write();
  hPPb_Unfo_3->Write();
  f.Close();

  



}
