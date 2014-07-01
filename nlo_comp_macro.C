#include <iostream>
#include <stdio.h>
#include <fstream>
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

static const int nbins_yaxian_large = 29;
static const double boundaries_yaxian_large[nbins_yaxian_large+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638,790,967};

//static const int nbins_yaxian_large = 25;
//static const double boundaries_yaxian_large[nbins_yaxian_large+1] = {47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638,790,967};

void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.425,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.425);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
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
	TH1F *hRebin = new TH1F(Form("%s_rebin",h->GetName()),Form("rebin %s",h->GetTitle()),nbins_recrebin,boundaries_recrebin);
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


void nlo_comp_macro(int radius = 3){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  
  TFile *fPP = TFile::Open("result-2013-akVs3PF-cent-1-isFineBin-0/pbpb_pp_merged_chmx_pt_isMC_0_Unfo_2013_akVs3PF_cent_1_isFineBin_0.root");
  //TFile *fNLO_err = TFile::Open("fnl4350a_cteq");
  TFile *fNLO_nnpdf = TFile::Open("fnl4350a_nnpdf21-nlo_aspdf_new.root");
  TFile *fNLO_cteq = TFile::Open("fnl4350a_cteq66-nlo_aspdf_all_new.root");
  TFile *fNLO_ct10n = TFile::Open("fnl4350a_ct10n-nlo_aspdf_new.root");
  TFile *fNLO_hera = TFile::Open("fnl4350a_hera15all-nlo_aspdf_new.root");
  TFile *fPP_data_R3 = TFile::Open("pp_2013_2760TeV_data_ak3PF.root");
  TFile *fPP_data_R4 = TFile::Open("pp_2013_2760TeV_data_ak4PF.root");
  TFile *fPP_data_R5 = TFile::Open("pp_2013_2760TeV_data_ak5PF.root");
  
  TFile *faditya = TFile::Open("comparison.root");

  //alright lets get the unfolded data here: remember we need it for eta range -2 to +2 
  TFile* fPP_unfo_R3 = TFile::Open("pp_2013_2760_abs_eta_2_mc_ak3PF.root");
  TFile* fPP_unfo_R4 = TFile::Open("pp_2013_2760_abs_eta_2_mc_ak4PF.root");
  TFile* fPP_unfo_R5 = TFile::Open("pp_2013_2760_abs_eta_2_mc_ak5PF.root");

  TFile fout("pp_2760GeV_nlo_histos.root","RECREATE");
  fout.cd();

  TH1F* hPP_nnpdf_NLO = (TH1F*)fNLO_nnpdf->Get("h100200");
  TH1F* hPP_cteq_NLO = (TH1F*)fNLO_cteq->Get("h100200");
  TH1F* hPP_ct10n_NLO = (TH1F*)fNLO_ct10n->Get("h100200");
  TH1F* hPP_hera_NLO = (TH1F*)fNLO_hera->Get("h100200");

  //TH1F* hPP_data_R_3 = (TH1F*)fPP_data_R3->Get("hppComb");
  //TH1F* hPP_data_R_4 = (TH1F*)fPP_data_R4->Get("hppComb");
  //TH1F* hPP_data_R_5 = (TH1F*)fPP_data_R5->Get("hppComb");

  TH1F* hPP_data_R_3 = (TH1F*)fPP_unfo_R3->Get("Unfolded_cent1");
  TH1F* hPP_data_R_4 = (TH1F*)fPP_unfo_R4->Get("Unfolded_cent1");
  TH1F* hPP_data_R_5 = (TH1F*)fPP_unfo_R5->Get("Unfolded_cent1");

  hPP_data_R_3->Scale(1./5300e6);
  hPP_data_R_3->Scale(1./4);
  divideBinWidth(hPP_data_R_3);

  hPP_data_R_4->Scale(1./5300e6);
  hPP_data_R_4->Scale(1./4);
  divideBinWidth(hPP_data_R_4);

  hPP_data_R_5->Scale(1./5300e6);
  hPP_data_R_5->Scale(1./4);
  divideBinWidth(hPP_data_R_5);

  // NLO histograms without any R# at the end correspond to R=0.3 the standard. others are named accordingly

  TH1F* hPP_nnpdf_NLO_R4 = (TH1F*)fNLO_nnpdf->Get("h100300");
  TH1F* hPP_nnpdf_NLO_R2 = (TH1F*)fNLO_nnpdf->Get("h100100");

  TH1F* hPP_cteq_NLO_R4 = (TH1F*)fNLO_cteq->Get("h100300");
  TH1F* hPP_ct10n_NLO_R4 = (TH1F*)fNLO_ct10n->Get("h100300");
  TH1F* hPP_hera_NLO_R4 = (TH1F*)fNLO_hera->Get("h100300");

  TH1F* hPP_cteq_NLO_R2 = (TH1F*)fNLO_cteq->Get("h100100");
  TH1F* hPP_ct10n_NLO_R2 = (TH1F*)fNLO_ct10n->Get("h100100");
  TH1F* hPP_hera_NLO_R2 = (TH1F*)fNLO_hera->Get("h100100");

  TH1F* hPP_err = (TH1F*)fNLO_cteq->Get("h100203");
  TH1F* hPP_err_R4 = (TH1F*)fNLO_cteq->Get("h100303");
  TH1F* hPP_err_R2 = (TH1F*)fNLO_cteq->Get("h100103");
  
  for(int i = 0;i<hPP_nnpdf_NLO->GetNbinsX();i++){
    
    Float_t valErr = hPP_err->GetBinError(i);
    hPP_nnpdf_NLO->SetBinError(i,valErr);
    hPP_cteq_NLO->SetBinError(i,valErr);
    hPP_hera_NLO->SetBinError(i,valErr);
    hPP_ct10n_NLO->SetBinError(i,valErr);

    Float_t valErr_R4 = hPP_err_R4->GetBinError(i);
    hPP_nnpdf_NLO_R4->SetBinError(i,valErr_R4);
    hPP_cteq_NLO_R4->SetBinError(i,valErr_R4);
    hPP_ct10n_NLO_R4->SetBinError(i,valErr_R4);
    hPP_hera_NLO_R4->SetBinError(i,valErr_R4);

    Float_t valErr_R2 = hPP_err_R2->GetBinError(i);
    hPP_nnpdf_NLO_R2->SetBinError(i,valErr_R2);
    hPP_cteq_NLO_R2->SetBinError(i,valErr_R2);
    hPP_ct10n_NLO_R2->SetBinError(i,valErr_R2);
    hPP_hera_NLO_R2->SetBinError(i,valErr_R2);

  }

  hPP_nnpdf_NLO->SetName("hPP_nnpdf_NLO");
  hPP_cteq_NLO->SetName("hPP_cteq_NLO");
  hPP_ct10n_NLO->SetName("hPP_ct10n_NLO");
  hPP_hera_NLO->SetName("hPP_hera_NLO");

  hPP_nnpdf_NLO_R2->SetName("hPP_nnpdf_NLO_R2");
  hPP_nnpdf_NLO_R2->Print("base");
  hPP_cteq_NLO_R2->SetName("hPP_cteq_NLO_R2");
  hPP_cteq_NLO_R2->Print("base");
  hPP_ct10n_NLO_R2->SetName("hPP_ct10n_NLO_R2");
  hPP_ct10n_NLO_R2->Print("base");
  hPP_hera_NLO_R2->SetName("hPP_hera_NLO_R2");

  hPP_nnpdf_NLO_R4->SetName("hPP_nnpdf_NLO_R4");
  hPP_cteq_NLO_R4->SetName("hPP_cteq_NLO_R4");
  hPP_ct10n_NLO_R4->SetName("hPP_ct10n_NLO_R4");
  hPP_hera_NLO_R4->SetName("hPP_hera_NLO_R4");

  hPP_nnpdf_NLO->Write();
  hPP_cteq_NLO->Write();
  hPP_ct10n_NLO->Write();
  hPP_hera_NLO->Write();
  hPP_nnpdf_NLO_R2->Write();
  hPP_cteq_NLO_R2->Write();
  hPP_ct10n_NLO_R2->Write();
  hPP_hera_NLO_R2->Write();
  hPP_nnpdf_NLO_R4->Write();
  hPP_cteq_NLO_R4->Write();
  hPP_ct10n_NLO_R4->Write();
  hPP_hera_NLO_R4->Write();

  fout.Write();

  TH1F* hPPrebin = (TH1F*)hPP_data_R_3->Clone("hPPrebin");
  //TH1F* hPPrebin_test = (TH1F*)faditya->Get("Corrected Jet Spectrum Aditya");
  //TH1F* hPPrebin = (TH1F*)hPPrebin_test->Rebin(nbins_yaxian_large,"hPPrebin",boundaries_yaxian_large);
  //TH1F* hPPunfo = (TH1F*)fPP->Get("Unfolded_cent6");
  TH1F* hPPgen = (TH1F*)fPP->Get("hGen_cent1");
  //hPPrebin->Scale(64);//remove the sigma scaling from the previous macros
  //dont need this now since we are taking it from a dedicated macro which gives us diff cross section 
  hPP_data_R_3->Scale(1e9);
  hPP_data_R_4->Scale(1e9);
  hPP_data_R_5->Scale(1e9);

  //hPP_nnpdf_NLO->Scale(1./4);
  //hPP_cteq_NLO->Scale(1./4);
  //hPP_hera_NLO->Scale(1./4);

  hPPrebin->Scale(1e9);
  //hPPrebin->Scale(1./5300e6);
  //hPPrebin->Scale(1./4);
  //divideBinWidth(hPPrebin);
 
  hPPgen->Scale(1./4);
  hPPgen->Scale(1e9);

  TH1F* hRatio_nnpdf = (TH1F*)hPP_nnpdf_NLO->Rebin(nbins_yaxian_large,"hRatio_nnpdf",boundaries_yaxian_large);
  hRatio_nnpdf->Divide(hPPrebin);

  TH1F* hRatio_cteq = (TH1F*)hPP_cteq_NLO->Rebin(nbins_yaxian_large,"hRatio_cteq",boundaries_yaxian_large);
  hRatio_cteq->Divide(hPPrebin);

  TH1F* hRatio_ct10n = (TH1F*)hPP_ct10n_NLO->Rebin(nbins_yaxian_large,"hRatio_ct10n",boundaries_yaxian_large);
  hRatio_ct10n->Divide(hPPrebin);

  TH1F* hRatio_hera = (TH1F*)hPP_hera_NLO->Rebin(nbins_yaxian_large,"hRatio_hera",boundaries_yaxian_large);
  hRatio_hera->Divide(hPPrebin);

  TH1F* hRatio_ppgen = (TH1F*)hPPgen->Rebin(nbins_yaxian_large,"hRatio_ppgen",boundaries_yaxian_large);
  hRatio_ppgen->Divide(hPPrebin);

  TH1F* hRatio_nnpdf_R_2_4 = (TH1F*)hPP_nnpdf_NLO_R2->Rebin(nbins_yaxian_large,"hRatio_nnpdf_R_2_4",boundaries_yaxian_large);
  hRatio_nnpdf_R_2_4->Divide(hPP_nnpdf_NLO_R4);

  TH1F* hRatio_nnpdf_R_3_4 = (TH1F*)hPP_nnpdf_NLO->Rebin(nbins_yaxian_large,"hRatio_nnpdf_R_3_4",boundaries_yaxian_large);
  hRatio_nnpdf_R_3_4->Divide(hPP_nnpdf_NLO_R4);

  TH1F* hRatio_cteq_R_2_4 = (TH1F*)hPP_cteq_NLO_R2->Rebin(nbins_yaxian_large,"hRatio_cteq_R_2_4",boundaries_yaxian_large);
  hRatio_cteq_R_2_4->Divide(hPP_cteq_NLO_R4);

  TH1F* hRatio_cteq_R_3_4 = (TH1F*)hPP_cteq_NLO->Rebin(nbins_yaxian_large,"hRatio_cteq_R_3_4",boundaries_yaxian_large);
  hRatio_cteq_R_3_4->Divide(hPP_cteq_NLO_R4);

  TH1F* hRatio_ct10n_R_2_4 = (TH1F*)hPP_ct10n_NLO_R2->Rebin(nbins_yaxian_large,"hRatio_ct10n_R_2_4",boundaries_yaxian_large);
  hRatio_ct10n_R_2_4->Divide(hPP_ct10n_NLO_R4);

  TH1F* hRatio_ct10n_R_3_4 = (TH1F*)hPP_ct10n_NLO->Rebin(nbins_yaxian_large,"hRatio_ct10n_R_3_4",boundaries_yaxian_large);
  hRatio_ct10n_R_3_4->Divide(hPP_ct10n_NLO_R4);

  TH1F* hRatio_hera_R_2_4 = (TH1F*)hPP_hera_NLO_R2->Rebin(nbins_yaxian_large,"hRatio_hera_R_2_4",boundaries_yaxian_large);
  hRatio_hera_R_2_4->Divide(hPP_hera_NLO_R4);

  TH1F* hRatio_hera_R_3_4 = (TH1F*)hPP_hera_NLO->Rebin(nbins_yaxian_large,"hRatio_hera_R_3_4",boundaries_yaxian_large);
  hRatio_hera_R_3_4->Divide(hPP_hera_NLO_R4);

  TH1F* hRatio_data_nnpdf_R_3 = (TH1F*)hPP_nnpdf_NLO->Rebin(nbins_yaxian_large,"hRatio_data_nnpdf_R_3",boundaries_yaxian_large);
  hRatio_data_nnpdf_R_3->Divide(hPP_data_R_3);
  
  TH1F* hRatio_data_nnpdf_R_4 = (TH1F*)hPP_nnpdf_NLO_R4->Rebin(nbins_yaxian_large,"hRatio_data_nnpdf_R_4",boundaries_yaxian_large);
  hRatio_data_nnpdf_R_4->Divide(hPP_data_R_4);

  TH1F* hRatio_data_R_3_4 = (TH1F*)hPP_data_R_3->Clone("hRatio_data_R_3_4");
  hRatio_data_R_3_4->Divide(hPP_data_R_4);

  TH1F* hRatio_data_R_3_5 = (TH1F*)hPP_data_R_3->Clone("hRatio_data_R_3_5");
  hRatio_data_R_3_5->Divide(hPP_data_R_5);

  TH1F* hRatio_data_R_4_5 = (TH1F*)hPP_data_R_4->Clone("hRatio_data_R_4_5");
  hRatio_data_R_4_5->Divide(hPP_data_R_5);



  TCanvas *c1 = new TCanvas("c1","",800,600);
  formatCanvas(c1);
  c1->cd(1);
  c1->cd(1)->SetLogy();
  hPP_nnpdf_NLO->SetMarkerColor(kRed);
  hPP_nnpdf_NLO->SetMarkerStyle(20);
  hPP_cteq_NLO->SetMarkerColor(kBlue);
  hPP_cteq_NLO->SetMarkerStyle(20);
  hPP_ct10n_NLO->SetMarkerColor(9);//purple
  hPP_ct10n_NLO->SetMarkerStyle(20);
  hPP_hera_NLO->SetMarkerColor(kGreen);
  hPP_hera_NLO->SetMarkerStyle(20);
  hPPrebin->SetMarkerColor(kBlack);
  hPPrebin->SetMarkerStyle(8);
  hPPgen->SetMarkerColor(kOrange);
  hPPgen->SetMarkerStyle(8);
  //hPPgen->SetMarkerColor(kRed);
  //hPPgen->SetMarkerStyle(8);
  
  hPP_nnpdf_NLO->SetYTitle("#frac{d^{2} #sigma}{d p_{T} d #eta} (pb#frac{GeV}{c})");
  hPP_nnpdf_NLO->SetXTitle("p_{T} (GeV/c)");
  hPP_nnpdf_NLO->SetAxisRange(22,500,"X");
  hPP_nnpdf_NLO->SetTitle(" ");
  hPP_nnpdf_NLO->Draw("p");
  hPP_cteq_NLO->Draw("same p");
  hPP_ct10n_NLO->Draw("same p");
  hPP_hera_NLO->Draw("same p");
  hPPrebin->Draw("same p");
  hPPgen->Draw("same p");

  TLegend * title = myLegend(0.47, 0.50,0.67, 0.8);
  title->AddEntry(hPP_nnpdf_NLO,"NLO nnpdf","pl");
  title->AddEntry(hPP_cteq_NLO,"NLO cteq","pl");
  title->AddEntry(hPP_ct10n_NLO,"NLO ct10n","pl");
  title->AddEntry(hPP_hera_NLO,"NLO hera","pl");
  title->AddEntry(hPPgen,"pp MC spectra","pl");
  title->AddEntry(hPPrebin,"pp unfolded 2013 data","pl");
  title->SetTextSize(0.04);
  title->Draw();

  putCMSPrel(0.1,0.92,0.06);
  drawText("pp 2013, #sqrt{s}=2.76(TeV), #int L dt = 5.3 (pb)^{-1}",0.35,0.92,16);
  drawText(Form("anti k_{T} R = 0.3",radius),0.47,0.83,16);

  c1->cd(2);
  hRatio_nnpdf->SetYTitle(" X / pp data");
  hRatio_nnpdf->SetXTitle("p_{T} (GeV/c)");
  hRatio_nnpdf->SetTitle(" ");
  hRatio_nnpdf->SetAxisRange(0,2,"Y");
  hRatio_nnpdf->SetMarkerColor(kRed);
  hRatio_nnpdf->SetMarkerStyle(20);
  hRatio_hera->SetMarkerColor(kGreen);
  hRatio_hera->SetMarkerStyle(20);
  hRatio_cteq->SetMarkerColor(kBlue);
  hRatio_cteq->SetMarkerStyle(20);
  hRatio_ct10n->SetMarkerColor(9);
  hRatio_ct10n->SetMarkerStyle(20);
  hRatio_ppgen->SetMarkerColor(kOrange);
  hRatio_ppgen->SetMarkerStyle(20);
  hRatio_nnpdf->Draw("p");
  hRatio_nnpdf->SetAxisRange(22,500,"X");
  hRatio_hera->Draw("same p");
  hRatio_cteq->Draw("same p");
  hRatio_ppgen->Draw("same p");

  c1->SaveAs(Form("pp_2760GeV_NLO_ak%dPF_vs_MC_gen_spectra_doga.pdf",radius),"RECREATE");

  //get the information from the ratio per bins - to use to scale down the NLO in 5.02 TeV. 
  ofstream R_nnpdf,R_hera,R_cteq;
  R_nnpdf.open(Form("ratio_nnpdf_vs_pp_data_2760_ak%d.txt",radius));
  R_hera.open(Form("ratio_hera_vs_pp_data_2760_ak%d.txt",radius));
  R_cteq.open(Form("ratio_cteq_vs_pp_data_2760_ak%d.txt",radius)); 
  for(int i = 0;i<hRatio_nnpdf->GetNbinsX();i++){
    R_nnpdf<<i<<"\t"<<hRatio_nnpdf->GetBinContent(i)<<endl;
    R_hera<<i<<"\t"<<hRatio_hera->GetBinContent(i)<<endl;
    R_cteq<<i<<"\t"<<hRatio_cteq->GetBinContent(i)<<endl;
  }
  
  R_nnpdf.close();
  R_hera.close();
  R_cteq.close();

  
  //draw the results for NLO comparison within different radius at the same energy 
  TCanvas *c2 = new TCanvas("c2","",800,600);
  formatCanvas(c2);
  
  c2->cd(1);
  hPP_nnpdf_NLO->SetMarkerStyle(22);
  hPP_nnpdf_NLO->SetMarkerColor(3);
  hPP_nnpdf_NLO->Draw("p");
  hPP_nnpdf_NLO_R2->SetMarkerStyle(22);
  hPP_nnpdf_NLO_R2->SetMarkerColor(2);
  hPP_nnpdf_NLO_R2->Draw("same p");
  hPP_nnpdf_NLO_R4->SetMarkerStyle(22);
  hPP_nnpdf_NLO_R4->SetMarkerColor(4);
  hPP_nnpdf_NLO_R4->Draw("same p");
  hPP_data_R_4->SetMarkerStyle(23);
  hPP_data_R_4->SetMarkerColor(4);
  hPP_data_R_4->Draw("same p");
  hPP_data_R_3->SetMarkerStyle(23);
  hPP_data_R_3->SetMarkerColor(3);
  hPP_data_R_3->Draw("same p");
  hPP_data_R_5->SetMarkerStyle(23);
  hPP_data_R_5->SetMarkerColor(9);
  hPP_data_R_5->Draw("same p");

  TLegend * title2 = myLegend(0.47, 0.50,0.67, 0.8);
  title2->AddEntry(hPP_nnpdf_NLO_R2,"NNPDF21 R=0.2","pl");
  title2->AddEntry(hPP_nnpdf_NLO,"NNPDF21 R=0.3","pl");
  title2->AddEntry(hPP_nnpdf_NLO_R4,"NNPDF21 R=0.4","pl");
  title2->AddEntry(hPP_data_R_3,"Data R=0.3","pl");
  title2->AddEntry(hPP_data_R_4,"Data R=0.4","pl");
  title2->AddEntry(hPP_data_R_5,"Data R=0.5","pl");
  title2->SetTextSize(0.04);
  title2->Draw();

  putCMSPrel(0.1,0.92,0.06);
  drawText("pp 2013, #sqrt{s}=2.76(TeV), #int L dt = 5.3 (pb)^{-1}",0.35,0.92,16);
  drawText(Form("anti k_{T}, Data vs Theory",radius),0.47,0.83,16);

  c2->cd(2);
  hRatio_nnpdf_R_2_4->SetMarkerStyle(29);
  hRatio_nnpdf_R_2_4->SetMarkerColor(2);
  hRatio_nnpdf_R_2_4->SetAxisRange(0.6,1.6,"Y");
  hRatio_nnpdf_R_2_4->SetAxisRange(22,500,"X");
  hRatio_nnpdf_R_2_4->SetTitle(" ");
  hRatio_nnpdf_R_2_4->SetYTitle("Ratios");
  hRatio_nnpdf_R_2_4->SetXTitle("p_{T}(GeV/c)");
  hRatio_nnpdf_R_2_4->Draw("p");
  hRatio_nnpdf_R_3_4->SetMarkerStyle(29);
  hRatio_nnpdf_R_3_4->SetMarkerColor(3);
  hRatio_nnpdf_R_3_4->Draw("same p");
  hRatio_data_nnpdf_R_3->SetMarkerStyle(33);
  hRatio_data_nnpdf_R_3->SetMarkerColor(4);
  hRatio_data_nnpdf_R_3->Draw("same p");
  hRatio_data_nnpdf_R_4->SetMarkerStyle(33);
  hRatio_data_nnpdf_R_4->SetMarkerColor(7);
  hRatio_data_nnpdf_R_4->Draw("same p");
  hRatio_data_R_3_4->SetMarkerStyle(34);
  hRatio_data_R_3_4->SetMarkerColor(6);
  hRatio_data_R_3_4->Draw("same p");
  hRatio_data_R_3_5->SetMarkerStyle(34);
  hRatio_data_R_3_5->SetMarkerColor(7);
  hRatio_data_R_3_5->Draw("same p");
  hRatio_data_R_4_5->SetMarkerStyle(34);
  hRatio_data_R_4_5->SetMarkerColor(8);
  hRatio_data_R_4_5->Draw("same p");

  c2->cd(1);
  TLegend * title3 = myLegend(0.67, 0.40,0.77, 0.8);
  title3->AddEntry(hRatio_nnpdf_R_2_4,"NLO R=0.2/R=0.4","pl");
  title3->AddEntry(hRatio_nnpdf_R_3_4,"NLO R=0.3/R=0.4","pl");
  title3->AddEntry(hRatio_data_nnpdf_R_3,"R=0.3 NNPDF21/Data","pl");
  title3->AddEntry(hRatio_data_nnpdf_R_4,"R=0.4 NNPDF21/Data","pl");
  title3->AddEntry(hRatio_data_R_3_4,"Data R=0.3/R=0.4","pl");
  title3->AddEntry(hRatio_data_R_3_5,"Data R=0.3/R=0.5","pl");
  title3->AddEntry(hRatio_data_R_4_5,"Data R=0.4/R=0.5","pl");
  title3->SetTextSize(0.04);
  title3->Draw();
  
  c2->SaveAs("pp_2760GeV_data_NLO_radius_comparison.pdf","RECREATE");

  fout.Write();
  fout.Close();

  /*
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  c2->Divide(2,1);
  c2->cd(1);
  c2->cd(1)->SetLogy();
  hPPunfo->SetTitle("PP 2013 2.76 TeV Merged, Unfolded p_{T} Spectra");
  hPPunfo->SetYTitle("#frac{dN}{N_{MB} d p_{T} d #eta}");
  hPPunfo->GetYaxis()->SetTitleOffset(1.4);
  hPPunfo->SetXTitle("Jet p_{T} GeV/c");
  hPPunfo->GetXaxis()->SetRangeUser(50,500);
  hPPunfo->Draw();
  hPPmeas->GetXaxis()->SetRangeUser(50,500);
  hPPmeas->SetMarkerStyle(22);
  hPPmeas->SetMarkerColor(kBlack);
  hPPmeas->Draw("same");
  TLegend *titl = myLegend(0.54,0.65,0.85,0.9);
  titl->AddEntry(hPPmeas,"PP2013 Meas ak3PF","pl");
  titl->AddEntry(hPPunfo,"PP2013 Unfo Bayesian","pl");
  titl->SetTextSize(0.03);
  titl->Draw();
  drawText("Anti-k_{T}PF R = 0.3",0.43,0.6,22);
  drawText("|#eta|<2, |vz|<15",0.47,0.5,22);
  
  c2->cd(2);
  TH1F* hPPRatio = (TH1F*)hPPmeas->Clone("hPPRatio");
  hPPRatio->Divide(hPPunfo);
  hPPRatio->SetXTitle("Jet p_{T} GeV/c");
  hPPRatio->SetYTitle("Ratio Measured/Unfolded");
  hPPRatio->SetTitle("PP2013 ak3PF merged");
  hPPRatio->GetYaxis()->SetRangeUser(0,2);
  hPPRatio->Draw();
  c2->SaveAs("pp_2013_ak3_merged_unfolded_pt.pdf","RECREATE");
  
  */
  
  /*
  TH1F *hNLO_err = (TH1F*)fNLO->Get("h100203");
  for(int i = 0;i<hNLO_err->GetNbinsX();i++){
    Float_T valErr = hNLO_Err->GetBinError(i);
    hNLO->SetBinError(i,valErr);
  }
  */
  /*
  //h100300 - ak4PF, h100200 - ak3PF
  TH1F* hNLO_2 = (TH1F*)hNLO->Clone("hNLO_2");
  hNLO->Print("base");

  TCanvas *cComp = new TCanvas("cComp","",800,600);
  cComp->Divide(2,1);
  cComp->cd(1);

  cComp->cd(1)->SetLogy();
  hPPrebin->SetMarkerStyle(22);
  hPPrebin->SetMarkerColor(kRed);

  hNLO->SetMarkerStyle(20);
  hNLO->SetMarkerColor(kBlack);
  hNLO->SetTitle("PP 2.76 TeV");
  hNLO->SetYTitle("#sigma pb");
  hNLO->SetXTitle("Jet p_{T} GeV/c");
  hNLO->Draw("pl");

  hPPMCrebin->SetMarkerStyle(21);
  hPPMCrebin->SetMarkerColor(kBlue);

  hPPMCrebin->Draw("same");

  hPPrebin->Draw("same");

  TLegend *title = myLegend(0.34,0.65,0.54,0.75);
  title->AddEntry(hNLO,"NLO nnpdf21","pl");
  title->AddEntry(hPPrebin,"PP ak3PF unfo","pl");
  title->AddEntry(hPPMCrebin,"PP ak3PF MC Gen","pl");
  title->SetTextSize(0.06);
  title->Draw();
  

  cComp->cd(2);
  TH1F *hPP = (TH1F*)hPPrebin->Clone("hPP");
  TH1F* hPPMC_v2 = (TH1F*)hPPMCrebin->Clone("hPPMC_v2");

  hPP->Print("base");
  hPPMC_v2->Print("base");
  
  hPP->Divide(hNLO);
  hPPMC_v2->Divide(hNLO);
  hPP->SetTitle("Ratio of PP ak3PF unfolded and MC Gen to NLO");
  //hPP->SetTitle("Ratio of PP ak4PF measured to NLO");
  hPP->SetYTitle("#frac{#sigma_{PP}}{#sigma_{NLO}}");
  hPP->SetXTitle("Jet p_{T} GeV/c");
  hPP->Draw();
  hPPMC_v2->Draw("same");

  cComp->SaveAs("fastNLO_comparison/ratio_pp_ak3_merged_NLO_nnpdf21nlo.pdf","RECREATE");

  TCanvas c1;
  c1.SetLogy();
  hPPrebin->SetAxisRange(50,500,"X");
  hPPrebin->SetYTitle("#sigma (pb)");
  hPPrebin->SetXTitle("Jet p_{T} GeV/c");
  hPPrebin->SetTitle("PP ak3PF unfolded");
  hPPrebin->Draw();
  c1.SaveAs("fastNLO_comparison/PP_2013_ak3_merged_Unfolded_crosssection.pdf","RECREATE");
  c1.SaveAs("fastNLO_comparison/PP_2013_ak3_merged_Unfolded_crosssection.C","RECREATE");
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetGrid();
  c3->SetLogy();
  //TGraphErrors graph_Expected("./fastNLO_comparison/files/s2760_R0.3fine.tex","%lg %lg");
  //graph_Expected.SetTitle("Non Perturbative corrections (extrapolation from Atlas) s2760;p_{T} GeV/c;#sigma nb");
  //graph_Expected.DrawClone("E3AL");
  

  TFile *fNPC = TFile::Open("fastNLO_comparison/files/npc_extrapolation_ivan.root");
  TH1F* hNPC = (TH1F*)fNPC->Get("hNPC");
  hNPC->Scale(1000);// 1000 for the pb from nb 
  //hNPC->Print("base");
  //TH1F* hNPC_rebin1 = rebin2(hNPC,"hNPC_rebin1");
  //divideBinWidth(hNPC_rebin1);
  hNPC->SetTitle("Non Perturbative corrections (extrapolation from Atlas) s2760");
  hNPC->SetXTitle("p_{T} GeV/c");
  hNPC->SetYTitle("#sigma (pb)");
  hNPC->Draw();
  hPPrebin->Draw("same");

  c3->SaveAs("fastNLO_comparison/NPC_atlas_ak3.pdf","RECREATE");


  TCanvas *c4 = new TCanvas("c4","",800,600);
  c4->SetLogy();
  //TH1F* hNPC_rebin = (TH1F*)hNPC->Clone("hNPC_rebin");
  TH1F* hNPC_rebin = rebin2(hNPC,"hNPC_rebin");
  divideBinWidth(hNPC_rebin);
  hPPrebin_2->SetTitle("PP Cross sections Data and Theory comparisons");
  hPPrebin_2->SetYTitle("#sigma (pb)");
  hPPrebin_2->SetXTitle("p_{T} GeV/c");
  //hNPC_rebin->SetAxisRange(1e4,1e-2,"Y");  
  hPPrebin_2->SetMarkerStyle(23);
  hPPrebin_2->SetMarkerColor(kBlack);
  hPPrebin_2->SetAxisRange(50,450,"X");
  hPPrebin_2->Draw("E");
  //hPPrebin_2->SetAxisRange(50,450,"X");
  hNPC_rebin->SetMarkerStyle(21);
  hNPC_rebin->SetMarkerColor(kRed);
  hNPC_rebin->Draw("same");
  hNLO_2->SetMarkerStyle(25);
  hNLO_2->SetMarkerColor(kBlue);
  hNLO_2->Draw("same");

  TLegend *title2 = myLegend(0.54,0.65,0.85,0.9);
  title2->AddEntry(hNPC_rebin,"Ivan NPC - Atlas R=0.3","pl");
  title2->AddEntry(hPPrebin_2,"PP2013 ak3PF unfolded","pl");
  title2->AddEntry(hNLO_2,"CMS NLO nnpdf21 R=0.3","l");

  title2->SetTextSize(0.04);
  title2->Draw();
  gStyle->SetOptStat(0);
 
  c4->SaveAs("fastNLO_comparison/pp_ak3_nlo_overlay_hist.pdf","RECREATE");

  TCanvas *c6 = new TCanvas("c6","",800,600);
  TH1F* hPPratio = (TH1F*)hPPrebin_2->Clone("hPPratio");

  //add the error bars directly before dividing them and then set that as the error. 
  //FLoat_t delta_PP = 0;
  //Float_t delta_NPC = 0;
  //Float_t delta
  
  hPPratio->Divide(hNPC_rebin);
  hPPratio->SetTitle("Ratio of PP 2013 ak3 unfolded w/ Ivan's NPC Atlas");
  hPPratio->SetYTitle(" ");
  hPPratio->SetXTitle("p_{T} GeV/c");
  hPPratio->Draw();
  c6->SaveAs("pp_ak3_npc_ratio.pdf","RECREATE");

  TCanvas *c5 = new TCanvas("c5","",800,600);
  hNPC->Draw();
  hNPC_rebin->Draw("same");
  c5->SetLogy();
  c5->SaveAs("Ivan_plot_rebin_ak3.pdf","RECREATE");
 
  */


}
