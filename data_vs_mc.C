// Raghav Kunnawalkam Elayavalli
// created April 11th 2014

// macro to check data vs mc. 


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

#include "headers/utilities_V0.h"

void formatCanvas(TCanvas *c){
  c->Divide(1,2,0.01,0.01);
  c->cd(1);
  c->GetPad(1)->SetLogy();
  c->GetPad(1)->SetPad(0.,0.425,1.,1.);
  c->GetPad(2)->SetPad(0.,0.0,1.,0.425);
  c->GetPad(2)->SetBottomMargin(0.3);
  c->GetPad(2)->SetGridy(1);
}

static const int nbins_yaxian = 29;
static const double boundaries_yaxian[nbins_yaxian+1] = {22, 27, 33, 39, 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638,790,967};

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

void * convertToInvYield(TH1 *hist) {
    for(int i = 1; i<=hist->GetNbinsX(); i++) {
        double content = hist->GetBinContent(i);
        double pt = hist->GetBinCenter(i);
        double error = hist->GetBinError(i);
        
        double new_content = content/(2.*TMath::Pi()*pt);
        double new_error = error/(2.*TMath::Pi()*pt);
        
        hist->SetBinContent(i,new_content);
        hist->SetBinError(i,new_error);
    }
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



void data_vs_mc(int radius = 2){

  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);

  // data files - pp 
  TFile *fpp1_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet80_v2.root");
  TFile *fpp2_v2 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_JEC_applied_ppJet40_v2.root");
  
  //do it for the pp - need to check up on this. 
  TTree *jetpp1_v2 = (TTree*)fpp1_v2->Get(Form("jetR%d",radius));
  TTree *jetpp2_v2 = (TTree*)fpp2_v2->Get(Form("jetR%d",radius));

  TTree *evtpp1_v2 = (TTree*)fpp1_v2->Get("evt");
  TTree *evtpp2_v2 = (TTree*)fpp2_v2->Get("evt");

  jetpp1_v2->AddFriend(evtpp1_v2);
  jetpp2_v2->AddFriend(evtpp2_v2);

  //get all the pp spectra here: 
  TCut pp3 = "abs(eta)>1.0&&abs(eta)<1.5&&jet40&&!jet60&&!jet80&&chMax/pt>0.01";
  
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
  
  jetpp1_v2->Project("hpp1","pt","abs(eta)>1.0&&abs(eta)<1.5&&jet80&&(chMax/pt)>0.01");
  hpp1->Print("base");
  
  jetpp2_v2->Project("hpp2","pt","abs(eta)>1.0&&abs(eta)<1.5&&jet60&&!jet80&&(chMax/pt)>0.01");
  hpp2->Print("base");
  
  jetpp2_v2->Project("hpp3","pt","9.25038"*pp3);
  // 9.25038 was the value. 
  //jetpp2_v2->Project("hpp3","pt","jet40_p"*pp3);
  hpp3->Print("base");

  Double_t deta = 1; //change this for every iteration. 
  
  hpp1->Scale(1./5300e6);//pp lumi
  hpp2->Scale(1./5300e6);
  hpp3->Scale(1./5300e6);

  hpp1->Scale(1./deta);//delta eta
  hpp2->Scale(1./deta);
  hpp3->Scale(1./deta);

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

  //hppComb->Scale(1e5);
 
  /*
  //load the MC files here: 
  TFile *fpp_mc = TFile::Open(Form("pp_2013_2760_abs_eta_10_15_mc_ak%dPF.root",radius));
  TH1F* hpp_mc = (TH1F*)fpp_mc->Get("hGen_cent1");

  hpp_mc->Print("base");
  hpp_mc->Scale(1./deta);//divide by delta eta
  divideBinWidth(hpp_mc);
  //hpp_mc->Scale(1e9);

  //looks like we dont need this now since the pp paper talks about differential cross section. 
  //convertToInvYield(hppComb);
  //convertToInvYield(hpp_mc);

  //draw the canvas here: 
  TCanvas *c = new TCanvas("c","",800,600);
  formatCanvas(c);
  c->cd(1);
  hpp_mc->SetMarkerStyle(24);
  hpp_mc->SetMarkerColor(kBlack);
  hppComb->SetMarkerStyle(20);
  hppComb->SetMarkerColor(kRed);

  hppComb->SetYTitle("#frac{d^{2} #sigma}{dp_{T} d#eta} ((mb)^{-1}/GeV/c)");
  hppComb->SetXTitle("Jet p_{T} (GeV/c)");
  hppComb->SetAxisRange(0,500,"X");
  hppComb->Draw();
  hpp_mc->Draw("same");

  TLegend *title = myLegend(0.47,0.55,0.67,0.8);
  title->AddEntry(hppComb,"Data","pl");
  title->AddEntry(hpp_mc,"Pythia","pl");
  title->SetTextSize(0.04);
  title->Draw();

  putCMSPrel(0.1,0.92,0.06);
  drawText("pp 2013, #int L dt = 5.3 (pb)^{-1}, 1.0<|#eta|<1.5",0.35,0.92,16);
  drawText(Form("anti k_{T} R = 0.%d, #sqrt{s} = 2.76 TeV",radius),0.47,0.83,16);

  c->cd(2);
  TH1F* hRatio = (TH1F*)hpp_mc->Clone("hRatio");
  hRatio->Divide(hppComb);
  hRatio->SetYTitle("#frac{Pythia}{Data}");
  hRatio->SetXTitle("Jet p_{T} (GeV/c)");
  hRatio->SetTitle(" ");
  hRatio->SetAxisRange(0,500,"X");
  hRatio->Draw();

  c->SaveAs(Form("pp_2013_data_vs_mc_10_15_ak%dPF.pdf",radius),"RECREATE");
  */

  TFile f(Form("pp_2013_2760TeV_data_ak%dPF.root",radius),"RECREATE");
  hppComb->Write();
  hpp_mc->Write();
  hRatio->Write();
  f.Write();
  f.Close();



}
