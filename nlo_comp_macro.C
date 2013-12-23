

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


void nlo_comp_macro(){

  TH1::SetDefaultSumw2();
  
  TFile *fPP = TFile::Open("result-2013-akPu4PF-cent-6-isFineBin-0/pbpb_pp_merged_chmx_pt_isMC_0_Unfo_2013_akPu4PF_cent_6_isFineBin_0.root");
  //TFile *fPP = TFile::Open("merge_pbpb_pp_ak5_HLT_V2.root");
  TFile *fNLO = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_nnpdf21-nlo_aspdf.root");
  //TFile *fNLO = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_cteq66-nlo_aspdf.root");
  TFile *fPP_meas = TFile::Open("merge_pbpb_pp_ak4_HLT_V2.root");//ak4
  //TFile *fPP_meas = TFile::Open("merge_HLT_V2.root");//ak3

  TH1F *htest_unfo = (TH1F*)fPP->Get("Unfolded_cent6");
  TH1F *htest_meas = (TH1F*)fPP_meas->Get("hppComb");
  TH1F *hPPunfo = (TH1F*)htest_unfo->Clone("hPPunfo");
  TH1F *hPPmeas = (TH1F*)htest_meas->Clone("hPPmeas");
  hPPmeas->Print("base");
  TH1F *hPPrebin_meas = rebin(htest_meas,"hPPrebin_meas");
  TH1F *hPPrebin = rebin(htest_unfo,"hPPrebin");
  //TH1F *htest = (TH1F*)fPP->Get("hppComb");
  //TH1F *hPPrebin = rebin(htest,"hPPrebin");

  hPPrebin->Scale(1./4);
  divideBinWidth(hPPrebin);
  hPPrebin->Scale(1./5.3);
  hPPmeas->Scale(1./4);
  hPPmeas->Scale(1./3.083e11);
  hPPunfo->Scale(1./4);
  hPPunfo->Scale(1./3.083e11);
  TH1F *hPPrebin_2 = rebin2(htest_unfo,"hPPrebin_2");
  hPPrebin_2->Scale(1./4);
  hPPrebin_2->Scale(1./5.3);
  divideBinWidth(hPPrebin_2);
  //TH1F *hPPrebin_2 = (TH1F*)hPPrebin->Clone("hPPrebin_2");
  gStyle->SetOptStat(0);

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
  titl->AddEntry(hPPmeas,"PP2013 Meas ak4PF","pl");
  titl->AddEntry(hPPunfo,"PP2013 Unfo Bayesian","pl");
  titl->SetTextSize(0.03);
  titl->Draw();
  drawText("Anti-k_{T}PF R = 0.4",0.43,0.6,22);
  drawText("|#eta|<2, |vz|<15",0.47,0.5,22);
  
  c2->cd(2);
  TH1F* hPPRatio = (TH1F*)hPPmeas->Clone("hPPRatio");
  hPPRatio->Divide(hPPunfo);
  hPPRatio->SetXTitle("Jet p_{T} GeV/c");
  hPPRatio->SetYTitle("Ratio Measured/Unfolded");
  hPPRatio->SetTitle("PP2013 ak4PF merged");
  hPPRatio->GetYaxis()->SetRangeUser(0,2);
  hPPRatio->Draw();
  c2->SaveAs("pp_2013_ak4_merged_unfolded_pt.gif","RECREATE");
  

  
  TH1F *hNLO = (TH1F*)fNLO->Get("h100300");
  //h100300 - ak4PF, h100200 - ak3PF
  TH1F* hNLO_2 = (TH1F*)hNLO->Clone("hNLO_2");

  TCanvas *cComp = new TCanvas("cComp","",800,600);
  cComp->Divide(2,1);
  cComp->cd(1);

  cComp->cd(1)->SetLogy();
  hPPrebin->SetMarkerStyle(22);
  hPPrebin->SetMarkerColor(kRed);

  hNLO->SetMarkerStyle(20);
  hNLO->SetMarkerColor(kBlack);
  hNLO->SetTitle("");
  hNLO->SetYTitle("#sigma pb");
  hNLO->SetXTitle("Jet p_{T} GeV/c");
  hNLO->Draw("pl");
  hPPrebin->Draw("same");

  TLegend *title = myLegend(0.54,0.65,0.85,0.9);
  title->AddEntry(hNLO,"NLO","pl");
  title->AddEntry(hPPrebin,"PP ak4PF unfo","pl");
  title->SetTextSize(0.06);
  title->Draw();
  

  cComp->cd(2);
  TH1F *hPP = (TH1F*)hPPrebin->Clone("hPP");
  hPP->Divide(hNLO);
  hPP->SetTitle("Ratio of PP ak4PF unfolded to NLO");
  //hPP->SetTitle("Ratio of PP ak4PF measured to NLO");
  hPP->SetYTitle("#frac{#sigma_{PP}}{#sigma_{NLO}}");
  hPP->SetXTitle("Jet p_{T} GeV/c");
  hPP->Draw();

  cComp->SaveAs("fastNLO_comparison/ratio_pp_ak4_merged_NLO_nnpdf21nlo.gif","RECREATE");

  TCanvas c1;
  c1.SetLogy();
  hPPrebin->SetAxisRange(50,500,"X");
  hPPrebin->SetYTitle("#sigma (pb)");
  hPPrebin->SetXTitle("Jet p_{T} GeV/c");
  hPPrebin->SetTitle("PP ak4PF unfolded");
  hPPrebin->Draw();
  c1.SaveAs("fastNLO_comparison/PP_2013_ak4_merged_Unfolded_crosssection.gif","RECREATE");
  c1.SaveAs("fastNLO_comparison/PP_2013_ak4_merged_Unfolded_crosssection.C","RECREATE");
  
  TCanvas *c3 = new TCanvas("c3","",800,600);
  c3->SetGrid();
  c3->SetLogy();
  //TGraphErrors graph_Expected("./fastNLO_comparison/files/s2760_R0.3fine.tex","%lg %lg");
  //graph_Expected.SetTitle("Non Perturbative corrections (extrapolation from Atlas) s2760;p_{T} GeV/c;#sigma nb");
  //graph_Expected.DrawClone("E3AL");
  

  TFile *fNPC = TFile::Open("fastNLO_comparison/files/npc_extrapolation_ivan_ak4.root");
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

  c3->SaveAs("fastNLO_comparison/NPC_atlas_ak4.gif","RECREATE");


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
  title2->AddEntry(hNPC_rebin,"Ivan NPC - Atlas R=0.4","pl");
  title2->AddEntry(hPPrebin_2,"PP2013 ak4PF unfolded","pl");
  title2->AddEntry(hNLO_2,"CMS NLO nnpdf21 R=0.4","l");

  title2->SetTextSize(0.04);
  title2->Draw();
  gStyle->SetOptStat(0);
 
  c4->SaveAs("fastNLO_comparison/pp_ak4_nlo_overlay_hist.pdf","RECREATE");

  TCanvas *c6 = new TCanvas("c6","",800,600);
  TH1F* hPPratio = (TH1F*)hPPrebin_2->Clone("hPPratio");

  //add the error bars directly before dividing them and then set that as the error. 
  //FLoat_t delta_PP = 0;
  //Float_t delta_NPC = 0;
  //Float_t delta
  
  hPPratio->Divide(hNPC_rebin);
  hPPratio->SetTitle("Ratio of PP 2013 ak4 unfolded w/ Ivan's NPC Atlas");
  hPPratio->SetYTitle(" ");
  hPPratio->SetXTitle("p_{T} GeV/c");
  hPPratio->Draw();
  c6->SaveAs("pp_ak4_npc_ratio.pdf","RECREATE");

  TCanvas *c5 = new TCanvas("c5","",800,600);
  hNPC->Draw();
  hNPC_rebin->Draw("same");
  c5->SetLogy();
  c5->SaveAs("Ivan_plot_rebin_ak4.pdf","RECREATE");
 
  


}
