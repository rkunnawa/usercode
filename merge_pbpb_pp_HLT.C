// macro to merge different HLT for pbpb(2011) and pp(2013)


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


void merge_pbpb_pp_HLT(){
  
  TH1::SetDefaultSumw2();
  /*
  TFile *fpbpb1 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet80.root");
  TFile *fpbpb2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet65.root");
  TFile *fpbpb3 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet55.root");
  */
  //TFile *fpp1 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet80.root");
  //TFile *fpp2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet40.root");
  
  TFile *fpp1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet80_v2.root");
  TFile *fpp2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet40_v2.root");

  /*
  TTree *jetpbpb1 = (TTree*)fpbpb1->Get("ntjet");
  TTree *jetpbpb2 = (TTree*)fpbpb2->Get("ntjet");
  TTree *jetpbpb3 = (TTree*)fpbpb3->Get("ntjet");
  */
  //TTree *jetpp1 = (TTree*)fpp1->Get("ntjet");
  // TTree *jetpp2 = (TTree*)fpp2->Get("ntjet");

  TTree *jetpp1_v2 = (TTree*)fpp1_v2->Get("jet");
  TTree *jetpp2_v2 = (TTree*)fpp2_v2->Get("jet");

  TTree *evtpp1_v2 = (TTree*)fpp1_v2->Get("evt");
  TTree *evtpp2_v2 = (TTree*)fpp2_v2->Get("evt");

  jetpp1_v2->AddFriend(evtpp1_v2);
  jetpp2_v2->AddFriend(evtpp2_v2);

  //TCut pbpb3 = "abs(eta)<2&&jet55&&!jet65&&!jet80";
  TCut pp3 = "abs(eta)<2&&jet40&&!jet60&&!jet80&&chMax/pt>0.01";
  /*
  TH1F *hpbpb1 = new TH1F("hpbpb1","",50,0,500);
  TH1F *hpbpb2 = new TH1F("hpbpb2","",50,0,500);
  TH1F *hpbpb3 = new TH1F("hpbpb3","",50,0,500);
  TH1F *hpbpbComb = new TH1F("hpbpbComb","",50,0,500);
  */
  TH1F *hpp1 = new TH1F("hpp1","",50,0,500);
  TH1F *hpp2 = new TH1F("hpp2","",50,0,500);
  TH1F *hpp3 = new TH1F("hpp3","",50,0,500);
  TH1F *hppComb = new TH1F("hppComb","",50,0,500);

  //get the prescl factor information. 
  //Float_t presclpbpb3 = (Float_t)jetpbpb1->GetEntries("abs(eta)<2")/jetpbpb1->GetEntries("abs(eta)<2&&jet55");
  //cout<<"pbpb prescl3 = "<<presclpbpb3<<endl;
  Float_t presclpp3 = (Float_t)jetpp1_v2->GetEntries("abs(eta)<2")/jetpp1_v2->GetEntries("abs(eta)<2&&jet40");
  cout<<"pp prescl3 = "<<presclpp3<<endl; //9.24968
  
  //jetpbpb1->Project("hpbpb1","pt","abs(eta)<2&&jet80");
  //hpbpb1->Print("base");
  //divideBinWidth(hpbpb1);

  //jetpbpb2->Project("hpbpb2","pt","abs(eta)<2&&jet65&&!jet80");
  //hpbpb2->Print("base");
  //divideBinWidth(hpbpb2);

  //jetpbpb3->Project("hpbpb3","pt","1.99942"*pbpb3);
  //hpbpb3->Print("base");
  //divideBinWidth(hpbpb3);

  jetpp1_v2->Project("hpp1","pt","abs(eta)<2&&jet80&&chMax/pt>0.01");
  hpp1->Print("base");
  //divideBinWidth(hpp1);
  //hpp1->Scale(1./3.083e11);
  //hpp1->Scale(1./4);

  jetpp2_v2->Project("hpp2","pt","abs(eta)<2&&jet60&&!jet80&&chMax/pt>0.01");
  hpp2->Print("base");
  //divideBinWidth(hpp2);
  //hpp2->Scale(1./3.083e11);
  //hpp2->Scale(1./4);

  jetpp2_v2->Project("hpp3","pt","9.24968"*pp3);
  hpp3->Print("base");
  //divideBinWidth(hpp3);
  //hpp3->Scale(1./3.083e11);
  //hpp3->Scale(1./4);

  //add the histograms
  /*
  hpbpbComb->Add(hpbpb1,1);
  hpbpbComb->Add(hpbpb2,1);
  hpbpbComb->Add(hpbpb3,1);

  hpbpbComb->Print("base");
  */
  hppComb->Add(hpp1,1);
  hppComb->Add(hpp2,1);
  hppComb->Add(hpp3,1);
  hppComb->Print("base");
  /*
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogy();
  hppComb->SetMarkerStyle(29);
  hppComb->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
  hppComb->SetXTitle("Jet p_{T} GeV/c");

  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]*pow(x+[1],[2])");
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Fit("fPowerLaw","","",50,500);
  hppComb->Draw();
  hpp3->SetMarkerStyle(24);
  hpp3->SetMarkerColor(kRed);
  hpp3->Draw("same");
  hpp2->SetMarkerStyle(25);
  hpp2->SetMarkerColor(kBlue);
  hpp2->Draw("same");
  hpp1->SetMarkerStyle(26);
  hpp1->SetMarkerColor(kGreen);
  hpp1->Draw("same");
  TLegend *title = myLegend(0.54,0.65,0.85,0.9);
  title->AddEntry(hppComb,"PP Merged","pl");
  title->AddEntry(hpp3,"PP HLT_40","pl");
  title->AddEntry(hpp2,"PP HLT_60","pl");
  title->AddEntry(hpp1,"PP HLT_80","pl");
  title->SetTextSize(0.06);
  title->Draw();
  drawText("PP 2013 Data PromptReco, pptracking and 2013 pp JEC",0.3,0.65,20);  
  drawText("Anti-k_{T} Particle Flow Jets R = 0.3, |#eta|<2, |vz|<15",0.3,0.56,20);
  c1->SaveAs("pp_2013_pt_combined.gif","RECREATE");

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
    
  */
  
  TCanvas c2;
  c2.SetLogy();
  TH1F* hPPComb = (TH1F*)hppComb->Clone("hPPComb");
  hPPComb->Scale(1./3.083e11);
  hPPComb->Print("base");

  hPPComb->Scale(1./4);
  divideBinWidth(hPPComb);

  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]/pow(x,[1])");
  hPPComb->Fit("fPowerLaw","","",50,500);
  hPPComb->Fit("fPowerLaw","","",50,500);
  hPPComb->Fit("fPowerLaw","","",50,500);
  hPPComb->Fit("fPowerLaw","","",50,500);
  hPPComb->Fit("fPowerLaw","","",50,500);
  hPPComb->Draw();

  c2.SaveAs("pp_2013_pt_evt_frac_merged.gif","RECREATE");
  
  TFile *fpbpbunfo = TFile::Open("result-2013-akPu3PF-cent-6-isFineBin-0/pbpb_pp_merged_chmx_pt_Unfo_2013_akPu3PF_cent_6_isFineBin_0.root");
  TH1F* hppUnfo = (TH1F*)fpbpbunfo->Get("Unfolded_cent6");
  hppUnfo->Print("base");
  
  TCanvas *c7 = new TCanvas("c7","",800,600);
  c7->SetLogy();
  hppComb->Draw();
  hppComb->SetYTitle("counts");
  hppComb->SetXTitle("p_{T} GeV/c");
  hppComb->SetTitle("PP 2013 2.76 TeV ak3PF measured vs unfolded");
  hppComb->SetMarkerStyle(23);
  hppComb->SetMarkerColor(kBlue);
  hppUnfo->SetAxisRange(10,500,"X");
  hppUnfo->SetMarkerStyle(24);
  hppUnfo->SetMarkerColor(kRed);
  hppUnfo->Draw("same");
  TLegend *title5 = myLegend(0.54,0.65,0.85,0.9);
  title5->AddEntry(hppComb,"Measured","pl");
  title5->AddEntry(hppUnfo,"Bayesian iter = 4","pl");
  title5->SetTextSize(0.06);
  title5->Draw();
  gStyle->SetOptStat(0);
  c7->SaveAs("PP2013_measured_vs_unfolded.gif","RECREATE");

  //Create output file and save them. 
  TFile f("merge_HLT_V2.root","RECREATE");
  //hpbpb1->Write();
  //hpbpb2->Write();
  //hpbpb3->Write();
  hpp1->Write();
  hpp2->Write();
  hpp3->Write();
  //hpbpbComb->Write();
  hppComb->Write();
  hPPComb->Write();
  f.Close();
  
  


}
