//macro for unfolded data to NLO comparison



#include <iostream>
#include <stdio.h>

#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <TStopwatch.h>

#include "headers/utilities.h"

using namespace std;

void data_nlo_comp(){

	
  //start timer
  TStopwatch timer;
  timer.Start();

  TH1::SetDefaultSumw2();

  

  //read the required input files
  TFile *fData = TFile::Open("result-2013-akPu3PF-cent6-isFineBin-0/pbpb_pp_chmx_pt_Unfo_2013_akPu3PF_cent6_isFineBin_0.root");

  //nlo
  TFile *fabkm095 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_abkm095-nlo_aspdf.root");
  TFile *fct10n = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_ct10n-nlo_aspdf.root");
  TFile *fhera15all = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_hera15all-nlo_aspdf.root");
  TFile *fnnpdf21 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_nnpdf21-nlo_aspdf.root");
  TFile *fabm115 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_abm115-nlo_aspdf.root");
  TFile *fcteq66 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_cteq66-nlo_aspdf.root");
  TFile *fmstw2008 = TFile::Open("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/mc/fnl4350_mstw2008-nlo_aspdf.root");
  
  //binning for the nlo
  //Int_t nbins = 11;
  //Double_t boundaries[nbins+1] = {100,110,120,130,140,150,160,170,180,200,240,300};

  //read the histograms
  TH1F *hPPUnfo = (TH1F*)fData->Get("Unfolded_cent6;1");
  hPPUnfo->Print("base");
  TH1F *hPPUnfo_rebin_0 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  TH1F *hPPUnfo_rebin_1 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  TH1F *hPPUnfo_rebin_2 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  TH1F *hPPUnfo_rebin_3 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  TH1F *hPPUnfo_rebin_4 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  TH1F *hPPUnfo_rebin_5 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  TH1F *hPPUnfo_rebin_6 = rebin(hPPUnfo,"PP 2013 Bayesian Unfolded p_T spectra");
  
  TH1F *habkm095 = (TH1F*)fabkm095->Get("h100200");
  habkm095->Print("base");
  
  TH1F *hct10n = (TH1F*)fct10n->Get("h100200");
  hcy10n->Print("base");

  TH1F *hhera15all = (TH1F*)fhera15all->Get("h100200");
  habkm095->Print("base");

  TH1F *hnnpdf21 = (TH1F*)fnnpdf21->Get("h100200");
  hnnpdf21->Print("base");

  TH1F *habm115 = (TH1F*)fabm115->Get("h100200");
  habm115->Print("base");

  TH1F *hcteq66 = (TH1F*)fcteq66->Get("h100200");
  hcteq66->Print("base");

  TH1F *hmstw2008 = (TH1F*)fmstw2008->Get("h100200");
  hmstw2008->Print("base");

  hPPUnfo->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo);

  hPPUnfo_rebin_0->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_0);

  hPPUnfo_rebin_1->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_1);

  hPPUnfo_rebin_2->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_2);

  hPPUnfo_rebin_3->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_3);

  hPPUnfo_rebin_4->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_4);

  hPPUnfo_rebin_5->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_5);

  hPPUnfo_rebin_6->Scale((5.3*e12)/(3.585*e11*4));
  divideBinWidth(hPPUnfo_rebin_6);

  


  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;




}
