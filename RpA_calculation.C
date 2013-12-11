// macro to read in the output from Unfold_RpPb_V0.C and calculate the RpA 


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


static const int nbins_recrebin = 11;
static const double boundaries_recrebin[nbins_recrebin+1] = {
  100,110,120,130,140,150,160,170,180,200,240,300
};

static const int nbins_recrebin_2 = 18;
static const double boundaries_recrebin_2[nbins_recrebin_2+1] = {
  50,60,70,80,90,100,110,120,130,140,150,160,170,180,200,240,300,360,420
};

const Double_t jetPtBin[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;


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

// rebin the spectra
TH1F *rebin(TH1F *h, char *histName)
{
	TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_recrebin_2,boundaries_recrebin_2);
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

TH1F *rebin2(TH1F *h, char *histName)
{
	TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nJetPtBin,jetPtBin);
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

void RpA_calculation(){

  TH1::SetDefaultSumw2();

  TFile* fin = TFile::Open("result-2013-ppb-akPu3PF-cent-1/ppb_merge_60_lowest_pp_mc_Unfo_akPu3PF_cent_1.root");
  TFile* fin_2 = TFile::Open("merge_ppb_60_lowest_HLT_V2.root");
  //TFile* fin_3 = TFile::Open("pPb_MB_pt.root");

  //get the required histograms.
  //formulae:: RpA = 1/N_evt d^2N/(dp_T d eta) * sigma_total
  //                 ---------------------------------------
  //                 N_coll * d^2 sigma/(dp_T d eta)

  //values: N_evt       = 6.53104 e10
  //        sigma_total = 70 e-3 b (70 mb)
  //        N_coll      = 6.9 

  // the pPb histogram from the file is already scaled with these factors 

  TH1F* hRpA_Meas = (TH1F*)fin_2->Get("hPPbComb_2");
  TH1F* hRpA_Unfo = (TH1F*)fin_2->Get("hPPb_Unfo_2");
  //TH1F* hpPb_MB = (TH1F*)fin_3->Get("hpPb_MB");

  hRpA_Meas->Print("base");
  hRpA_Unfo->Print("base");
  //hpPb_MB->Print("base");
  
  //this is undividing by bin width
  //hRpA_Meas->Scale(10);
  //hRpA_Unfo->Scale(10);
  
  //first undivide by bin width and then rebin, and then divide by bin width. 

  TH1F* hRpA_Meas_rebin = rebin(hRpA_Meas,"hRpA_Meas_rebin");
  TH1F* hRpA_Meas_rebin2 = rebin2(hRpA_Meas,"hRpA_Meas_rebin2");
  TH1F* hRpA_Unfo_rebin = rebin(hRpA_Unfo,"hRpA_Unfo_rebin");
  TH1F* hRpA_Unfo_rebin2 = rebin2(hRpA_Unfo,"hRpA_Unfo_rebin2");

  //divide by bin width
  //divideBinWidth(hRpA_Meas_rebin);
  //divideBinWidth(hRpA_Meas_rebin2);
  //divideBinWidth(hRpA_Unfo_rebin);
  //divideBinWidth(hRpA_Unfo_rebin2);

  //plot the spectra with measured vs unfolded. 
  TH1F* hpPbmeas = (TH1F*)hRpA_Meas->Clone("hpPbmeas");
  TH1F* hpPbunfo = (TH1F*)hRpA_Unfo->Clone("hpPbunfo");

  //divideBinWidth(hpPbmeas);
  //divideBinWidth(hpPbunfo);

  TCanvas *c = new TCanvas("c","",1000,800);
  c->Divide(2,1);
  c->cd(1);
  c->cd(1)->SetLogy();
  hpPbunfo->SetTitle("pPb 2013 5.02 TeV Merged, Unfolded p_{T} Spectra");
  hpPbunfo->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
  hpPbunfo->GetYaxis()->SetTitleOffset(1.4);
  hpPbunfo->SetXTitle("Jet p_{T} GeV/c");
  hpPbunfo->GetXaxis()->SetRangeUser(50,500);
  hpPbunfo->Draw();
  hpPbmeas->GetXaxis()->SetRangeUser(50,500);
  hpPbmeas->SetMarkerStyle(22);
  hpPbmeas->SetMarkerColor(kBlack);
  hpPbmeas->Draw("same");
  TLegend *titl = myLegend(0.54,0.65,0.85,0.9);
  titl->AddEntry(hpPbmeas,"pPb2013 Meas akPu3PF","pl");
  titl->AddEntry(hpPbunfo,"pPb2013 Unfo Bayesian","pl");
  titl->SetTextSize(0.03);
  titl->Draw();
  drawText("Anti-k_{T}Pu PF R = 0.3",0.43,0.6,22);
  drawText("|#eta|<2, |vz|<15",0.47,0.5,22);
  
  c->cd(2);
  TH1F* hpPbRatio = (TH1F*)hpPbmeas->Clone("hpPbRatio");
  hpPbRatio->Divide(hpPbunfo);
  hpPbRatio->SetXTitle("Jet p_{T} GeV/c");
  hpPbRatio->SetYTitle("Ratio Measured/Unfolded");
  hpPbRatio->SetTitle("pPb2013 akPu3PF merged");
  hpPbRatio->GetYaxis()->SetRangeUser(0,2);
  hpPbRatio->Draw();
  c->SaveAs("pPb_2013_akPu3PF_merged_unfolded_pt.gif","RECREATE");

  TH1F* hPP_Refe = (TH1F*)fin->Get("hGen_cent1");
  hPP_Refe->Print("base");
  hPP_Refe->Scale(1./70);
  hPP_Refe->Scale(6.9);
  hPP_Refe->Scale(1./10);//p_T width
  hPP_Refe->Scale(1./2);//RAPIDITY

  TH1F* hPP_Refe_rebin = rebin(hPP_Refe,"hPP_Refe_rebin");
  TH1F* hPP_Refe_rebin2 = rebin2(hPP_Refe,"hPP_Refe_rebin2");


  TCanvas *c1 = new TCanvas("c1","",800,600);
  hRpA_Meas->Divide(hPP_Refe);
  hRpA_Unfo->Divide(hPP_Refe);
  
  hRpA_Meas->SetMarkerStyle(23);
  hRpA_Unfo->SetMarkerStyle(24);

  hRpA_Meas->SetMarkerColor(kBlack);
  hRpA_Unfo->SetMarkerColor(kRed);

  hRpA_Meas->SetYTitle("RpA");
  hRpA_Meas->SetXTitle("Jet p_{T} GeV");
  hRpA_Meas->Draw();
  hRpA_Unfo->Draw("same");

  c1->SaveAs("RpA_merge_60_lowest_calculation_bin_v0.gif","RECREATE");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  hRpA_Meas_rebin->Divide(hPP_Refe_rebin);
  hRpA_Unfo_rebin->Divide(hPP_Refe_rebin);
  
  hRpA_Meas_rebin->SetMarkerStyle(23);
  hRpA_Unfo_rebin->SetMarkerStyle(24);

  hRpA_Meas_rebin->SetMarkerColor(kBlack);
  hRpA_Unfo_rebin->SetMarkerColor(kRed);

  hRpA_Meas_rebin->SetYTitle("RpA");
  hRpA_Meas_rebin->SetXTitle("Jet p_{T} GeV");
  hRpA_Meas_rebin->Draw();
  hRpA_Unfo_rebin->Draw("same");

  c2->SaveAs("RpA_merge_60_lowest_calculation_bin_v1.gif","RECREATE");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  hRpA_Meas_rebin2->Divide(hPP_Refe_rebin2);
  hRpA_Unfo_rebin2->Divide(hPP_Refe_rebin2);
  
  hRpA_Meas_rebin2->SetMarkerStyle(23);
  hRpA_Unfo_rebin2->SetMarkerStyle(24);

  hRpA_Meas_rebin2->SetMarkerColor(kBlack);
  hRpA_Unfo_rebin2->SetMarkerColor(kRed);

  hRpA_Meas_rebin2->SetYTitle("RpA");
  hRpA_Meas_rebin2->SetXTitle("Jet p_{T} GeV");
  hRpA_Meas_rebin2->Draw();
  hRpA_Unfo_rebin2->Draw("same");

  c3->SaveAs("RpA_merge_60_lowest_calculation_bin_v2.gif","RECREATE");



}
