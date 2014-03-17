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

// rebin the spectra
TH1F *rebin_yaxian(TH1F *h, char *histName)
{
  TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_jetPtBin,boundaries_jetPtBin);
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

void merge_pbpb_pp_HLT(){
  
  TH1::SetDefaultSumw2();
  
  //TFile *fpbpb1 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet80_v2.root");
  //TFile *fpbpb2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet65_v2.root");
  //TFile *fpbpb3 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet55_v2.root");
  
  //TFile *fpp1 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet80.root");
  //TFile *fpp2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet40.root");
  
  TFile *fpp1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet80_v2.root");
  TFile *fpp2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet40_v2.root");

  /*
  TTree *jetpbpb1_v2 = (TTree*)fpbpb1->Get("jetR3");
  TTree *jetpbpb2_v2 = (TTree*)fpbpb2->Get("jetR3");
  TTree *jetpbpb3_v2 = (TTree*)fpbpb3->Get("jetR3");

  TTree *evtpbpb1_v2 = (TTree*)fpbpb1->Get("evt");
  TTree *evtpbpb2_v2 = (TTree*)fpbpb2->Get("evt");
  TTree *evtpbpb3_v2 = (TTree*)fpbpb3->Get("evt");  

  jetpbpb1_v2->AddFriend(evtpbpb1_v2);
  jetpbpb2_v2->AddFriend(evtpbpb2_v2);
  jetpbpb3_v2->AddFriend(evtpbpb3_v2);
  */
  //TTree *jetpp1 = (TTree*)fpp1->Get("ntjet");
  // TTree *jetpp2 = (TTree*)fpp2->Get("ntjet");

  TTree *jetpp1_v2 = (TTree*)fpp1_v2->Get("jetR3");
  TTree *jetpp2_v2 = (TTree*)fpp2_v2->Get("jetR3");

  TTree *evtpp1_v2 = (TTree*)fpp1_v2->Get("evt");
  TTree *evtpp2_v2 = (TTree*)fpp2_v2->Get("evt");

  jetpp1_v2->AddFriend(evtpp1_v2);
  jetpp2_v2->AddFriend(evtpp2_v2);

  //TCut pbpb3 = "abs(eta)<2&&jet55&&!jet65&&!jet80&&chMax/pt>0.01";
  TCut pp3 = "abs(eta)<2&&jet40&&!jet60&&!jet80&&chMax/pt>0.01";
  
  //TH1F *hpbpb1 = new TH1F("hpbpb1","",30,0,300);
  //TH1F *hpbpb2 = new TH1F("hpbpb2","",30,0,300);
  //TH1F *hpbpb3 = new TH1F("hpbpb3","",30,0,300);
  //TH1F *hpbpbComb = new TH1F("hpbpbComb","",30,0,300);
  
  TH1F *hpp1 = new TH1F("hpp1","",nbins_jetPtBin,boundaries_jetPtBin);
  TH1F *hpp2 = new TH1F("hpp2","",nbins_jetPtBin,boundaries_jetPtBin);
  TH1F *hpp3 = new TH1F("hpp3","",nbins_jetPtBin,boundaries_jetPtBin);
  TH1F *hppComb = new TH1F("hppComb","",nbins_jetPtBin,boundaries_jetPtBin);
  
  //get the prescl factor information. 
  //Float_t presclpbpb3 = (Float_t)jetpbpb1_v2->GetEntries("jet80")/jetpbpb1_v2->GetEntries("jet55&&jet80");
  //cout<<"pbpb prescl3 = "<<presclpbpb3<<endl;//1.99871
  Float_t presclpp3 = (Float_t)jetpp1_v2->GetEntries("jet80")/jetpp1_v2->GetEntries("jet40&&jet80");
  cout<<"pp prescl3 = "<<presclpp3<<endl; //9.24968
  /*
  jetpbpb1_v2->Project("hpbpb1","pt","abs(eta)<2&&jet80&&chMax/pt>0.01");
  hpbpb1->Print("base");
  divideBinWidth(hpbpb1);

  jetpbpb2_v2->Project("hpbpb2","pt","abs(eta)<2&&jet65&&!jet80&&chMax/pt>0.01");
  hpbpb2->Print("base");
  divideBinWidth(hpbpb2);

  jetpbpb3_v2->Project("hpbpb3","pt","1.9987"*pbpb3);
  hpbpb3->Print("base");
  divideBinWidth(hpbpb3);
  */
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

  jetpp2_v2->Project("hpp3","pt","9.25038"*pp3);
  //9.25038 - ak5
  //9.24955 - ak4
  //9.24968 - ak3
  hpp3->Print("base");
  //divideBinWidth(hpp3);
  //hpp3->Scale(1./3.083e11);
  //hpp3->Scale(1./4);

  //scale the PbPb histograms before adding them
  //we have to scale them according to the lumi of the Jet80 file. 
  // HLT file  |   Lumi
  // HLT_80    |   150 mub-1
  // HLT_65    |   12.1 mub-1
  // HLT_55    |   0.38 mub-1
  // 
  // therefore scale for HLT_55 = 150/0.38 = 394.73684
  // scale for HLT_65 = 150/12.1 = 12.3967

  //hpbpb2->Scale(12.3867);
  //hpbpb3->Scale(394.7368);

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
  hpbpbComb->SetMarkerStyle(29);
  //hpbpbComb->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
  hpbpbComb->SetYTitle("counts");
  hpbpbComb->SetXTitle("Jet p_{T} GeV/c");

  TF1 *fPowerLaw = new TF1("fPowerLaw","[0]*pow(x+[1],[2])");
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Fit("fPowerLaw","","",30,300);
  hpbpbComb->Draw();
  hpbpb3->SetMarkerStyle(24);
  hpbpb3->SetMarkerColor(kRed);
  hpbpb3->Draw("same");
  hpbpb2->SetMarkerStyle(25);
  hpbpb2->SetMarkerColor(kBlue);
  hpbpb2->Draw("same");
  hpbpb1->SetMarkerStyle(26);
  hpbpb1->SetMarkerColor(kGreen);
  hpbpb1->Draw("same");
  TLegend *title = myLegend(0.54,0.65,0.85,0.9);
  title->AddEntry(hpbpbComb,"PbPb Merged","pl");
  title->AddEntry(hpbpb3,"w_{3} * (HLT_55 && !HLT_65 && !HLT_80)","pl");
  title->AddEntry(hpbpb2,"HLT_65 && !HLT_80","pl");
  title->AddEntry(hpbpb1,"HLT_80","pl");
  title->SetTextSize(0.03);
  title->Draw();
  drawText("PbPb 2011, 55,65 scaled",0.3,0.65,20);  
  drawText("Anti-k_{T} PU PF Jets R = 0.3, |#eta|<2, |vz|<15",0.3,0.56,20);
  c1->SaveAs("pbpb_2013_pt_combined.gif","RECREATE");
  */

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
  TFile f("merge_pp_ak3_HLT_V2.root","RECREATE");
  //hpbpb1->Write();
  //hpbpb2->Write();
  //hpbpb3->Write();
  //hPPComb_bins->Write();
  hpp1->Write();
  hpp2->Write();
  hpp3->Write();
  //hpbpbComb->Write();
  hppComb->Write();
  //hPPComb->Write();
  //hPbPb_Uncert->Write();
  hPPComb->Write();
  hPPGen->Write();
  f.Close();
  
  


}
