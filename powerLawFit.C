// macro to do the power law fit for data, MC and unfolded spectra 

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
#include <TStopwatch.h>

#include "RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"

#include "RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "headers/utilities.h"

using namespace std;



void powerLawFit(){
	
  //start timer
  TStopwatch timer;
  timer.Start();

  TH1::SetDefaultSumw2();

  //get the histograms from the output file of Unfold_RAA_V0.C 
  TFile *fin = TFile::Open("result-2013-akPu3PF-cent-6-isFineBin-0/pbpb_pp_chmx_pt_Unfo_2013_akPu3PF_cent_6_isFineBin_0.root");

  int  year = 2013;
  int algo =3;
  int isFineBin = 0;

  for(int i = 0;i<=nbins_cent;i++){

    if(i<nbins_cent){ //for PbPb
      cout<<endl<<endl<<endl<<endl<<"inside the power law fit loop in centrality i = "<<i<<endl<<endl<<endl<<endl;
      // Data
      
      TCanvas *cPowerLawData = new TCanvas("cPowerLawData","Data",800,600);
      cPowerLawData->Divide(2,1);
      cPowerLawData->cd(1);
      cPowerLawData->cd(1)->SetLogy();
      //TH1F *hTempData = (TH1F*)uhist[i]->hMeas->Clone(Form("Measured Pt_spectra PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      TH1F *hTempData = (TH1F*)fin->Get(Form("hMeas_cent%d",i));
      hTempData->SetTitle(Form("Measured Pt cpectra PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      TH1F *hRatioData = (TH1F*)hTempData->Clone("hRatioData");
      hTempData->SetAxisRange(100,330,"X");
      hTempData->SetMarkerStyle(22);
      hTempData->SetMarkerColor(kBlue);
      TF1 *fPowerLaw1 = new TF1("fPowerLaw1","[0]*pow(x+[2],[1])");
      fPowerLaw1->SetParameters(1e14,-5,0);
      //fPowerLaw->SetParameters(1e10,-8.8,40);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      TH1F *hFuncData = (TH1F*)functionHist(fPowerLaw1,hTempData,Form("Fit function Pt data PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      hTempData->SetXTitle("Jet p_{T} GeV/c");
      fPowerLaw1->SetMarkerStyle(8);
      fPowerLaw1->SetMarkerColor(4);
      fPowerLaw1->Draw("same");
      cPowerLawData->cd(2);
      hRatioData->SetAxisRange(100,330,"X");
      hRatioData->Divide(hFuncData);
      hRatioData->SetTitle("Spectra to Fit Ratio");
      hRatioData->SetXTitle("jet p_{T} GeV/c");
      hRatioData->SetYTitle("data/fit");
      hRatioData->SetMarkerColor(4);
      hRatioData->SetMarkerStyle(8);
      hRatioData->Draw();
      //cPowerLawData->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PbPb_Data_cent_%d_%d_.gif",year,algoName[algo],nbins_cent,isFineBin,boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5),"RECREATE");
      cPowerLawData->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PbPb_Data_cent_%d.gif",year,algoName[algo],nbins_cent,isFineBin,i),"RECREATE");
      cPowerLawData->Destructor();
     
      
      // MC have to check if it turns out ok, if not then have to change the starting parameters for the fit function. 
      TCanvas *cPowerLawMC = new TCanvas("cPowerLawMC","MC",800,600);
      cPowerLawMC->Divide(2,1);
      cPowerLawMC->cd(1);
      cPowerLawMC->cd(1)->SetLogy();
      //TH1F *hTempMC = (TH1F*)uhist[i]->hGen->Clone(Form("Generated Pt_spectra PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      TH1F *hTempMC = (TH1F*)fin->Get(Form("hGen_cent%d",i));
      hTempMC->SetTitle(Form("Generated Pt cpectra PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      TH1F *hRatioMC = (TH1F*)hTempMC->Clone("hRatioMC");
      hTempMC->SetAxisRange(100,330,"X");
      hTempMC->SetMarkerStyle(22);
      hTempMC->SetMarkerColor(kBlue);
      // trying matt's idea of [0]/(pow(x,[1]))
      TF1 *fPowerLaw2 = new TF1("fPowerLaw2","[0]/(pow(x,[1]))");
      // best fit so far
      //TF1 *fPowerLaw2 = new TF1("fPowerLaw2","[0]*pow(x+[2],[1])");
      //fPowerLaw2->SetParameters(1e10,-8.8,40);
      // sevil's idea
      //fPowerLaw2->SetParLimits(0,1e-5,1e-2);
      //fPowerLaw2->SetParLimits(1,-5,0);
      //fPowerLaw2->SetParLimits(2,-10,0);
      hTempMC->Fit("fPowerLaw2","LL","",100,330);
      hTempMC->Fit("fPowerLaw2","LL","",100,330);
      hTempMC->Fit("fPowerLaw2","LL","",100,330);
      hTempMC->Fit("fPowerLaw2","LL","",100,330);
      hTempMC->Fit("fPowerLaw2","LL","",100,330);
      hTempMC->Fit("fPowerLaw2","LL","",100,330);
      TH1F *hFuncMC = (TH1F*)functionHist(fPowerLaw2,hTempMC,Form("Fit function Pt MC PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      hTempMC->SetXTitle("Jet p_{T} GeV/c");
      fPowerLaw2->SetMarkerStyle(8);
      fPowerLaw2->SetMarkerColor(4);
      fPowerLaw2->Draw("same");
      cPowerLawMC->cd(2);
      hRatioMC->SetAxisRange(100,330,"X");
      hRatioMC->Divide(hFuncMC);
      hRatioMC->SetTitle("Spectra to Fit Ratio");
      hRatioMC->SetXTitle("jet p_{T} GeV/c");
      hRatioMC->SetYTitle("mc/fit");
      hRatioMC->SetMarkerColor(4);
      hRatioMC->SetMarkerStyle(8);
      hRatioMC->Draw();
      //cPowerLawMC->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PbPb_MC_cent_%d_%d_.gif",year,algoName[algo],nbins_cent,isFineBin,boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5),"RECREATE");
      cPowerLawMC->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PbPb_MC_cent_%d.gif",year,algoName[algo],nbins_cent,isFineBin,i),"RECREATE");
      cPowerLawMC->Destructor();
      
      
      // Reco
      TCanvas *cPowerLawReco = new TCanvas("cPowerLawReco","Reco",800,600);
      cPowerLawReco->Divide(2,1);
      cPowerLawReco->cd(1);
      cPowerLawReco->cd(1)->SetLogy();
      //TH1F *hTempReco = (TH1F*)uhist[i]->hReco->Clone(Form("Unfolded Pt_spectra PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      TH1F *hTempReco = (TH1F*)fin->Get(Form("Unfolded_cent%d",i));
      hTempReco->SetTitle(Form("Unfolded Pt cpectra PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      TH1F *hRatioReco = (TH1F*)hTempReco->Clone("hRatioReco");
      hTempReco->SetAxisRange(100,330,"X");
      hTempReco->SetMarkerStyle(22);
      hTempReco->SetMarkerColor(kBlue);
      TF1 *fPowerLaw3 = new TF1("fPowerLaw3","[0]*pow(x+[2],[1])");
      fPowerLaw3->SetParameters(1e14,-5,0);
      //fPowerLaw->SetParameters(1e10,-8.8,40);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      TH1F *hFuncReco = (TH1F*)functionHist(fPowerLaw3,hTempReco,Form("Fit function Pt Reco PbPb cent %2.0f-%2.0f%%",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      hTempReco->SetXTitle("Jet p_{T} GeV/c");
      fPowerLaw3->SetMarkerStyle(8);
      fPowerLaw3->SetMarkerColor(4);
      fPowerLaw3->Draw("same");
      cPowerLawReco->cd(2);
      hRatioReco->SetAxisRange(100,330,"X");
      hRatioReco->Divide(hFuncReco);
      hRatioReco->SetTitle("Spectra to Fit Ratio");
      hRatioReco->SetXTitle("jet p_{T} GeV/c");
      hRatioReco->SetYTitle("reco/fit");
      hRatioReco->SetMarkerColor(4);
      hRatioReco->SetMarkerStyle(8);
      hRatioReco->Draw();
      cPowerLawReco->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PbPb_Reco_cent_%d.gif",year,algoName[algo],nbins_cent,isFineBin,i),"RECREATE");
      cPowerLawReco->Destructor();
      
    }

    if(i == nbins_cent){ // for PP 
      cout<<endl<<endl<<endl<<endl<<"inside the power law fit loop in centrality i = "<<i<<endl<<endl<<endl<<endl;
      
      //Data
      TCanvas *cPowerLawData = new TCanvas("cPowerLawData","Data",800,600);
      cPowerLawData->Divide(2,1);
      cPowerLawData->cd(1);
      cPowerLawData->cd(1)->SetLogy();
      //TH1F *hTempData = (TH1F*)uhist[i]->hMeas->Clone("Measured Pt_spectra PP");
      TH1F *hTempData = (TH1F*)fin->Get(Form("hMeas_cent%d",i));
      hTempData->SetTitle("Measured Pt cpectra PP");
      TH1F *hRatioData = (TH1F*)hTempData->Clone("hRatioData");
      hTempData->SetAxisRange(100,330,"X");
      hTempData->SetMarkerStyle(22);
      hTempData->SetMarkerColor(kBlue);
      TF1 *fPowerLaw1 = new TF1("fPowerLaw1","[0]*pow(x+[2],[1])");
      fPowerLaw1->SetParameters(1e14,-5,0);
      //fPowerLaw->SetParameters(1e10,-8.8,40);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      hTempData->Fit("fPowerLaw1","LL","",100,330);
      TH1F *hFuncData = (TH1F*)functionHist(fPowerLaw1,hTempData,"Fit function Pt data PP");
      hTempData->SetXTitle("Jet p_{T} GeV/c");
      fPowerLaw1->SetMarkerStyle(8);
      fPowerLaw1->SetMarkerColor(4);
      fPowerLaw1->Draw("same");
      cPowerLawData->cd(2);
      hRatioData->SetAxisRange(100,330,"X");
      hRatioData->Divide(hFuncData);
      hRatioData->SetTitle("Spectra to Fit Ratio");
      hRatioData->SetXTitle("jet p_{T} GeV/c");
      hRatioData->SetYTitle("data/fit");
      hRatioData->SetMarkerColor(4);
      hRatioData->SetMarkerStyle(8);
      hRatioData->Draw();
      cPowerLawData->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PP_Data.gif",year,algoName[algo],nbins_cent,isFineBin),"RECREATE");
      
      // MC
      TCanvas *cPowerLawMC = new TCanvas("cPowerLawMC","MC",800,600);
      cPowerLawMC->Divide(2,1);
      cPowerLawMC->cd(1);
      cPowerLawMC->cd(1)->SetLogy();
      //TH1F *hTempMC = (TH1F*)uhist[i]->hGen->Clone("Generated Pt_spectra PP");
      TH1F *hTempMC = (TH1F*)fin->Get(Form("hGen_cent%d",i));
      hTempMC->SetTitle("Generated p_T spectra PP");
      TH1F *hRatioMC = (TH1F*)hTempMC->Clone("hRatioMC");
      hTempMC->SetAxisRange(100,330,"X");
      hTempMC->SetMarkerStyle(22);
      hTempMC->SetMarkerColor(kBlue);
      TF1 *fPowerLaw2 = new TF1("fPowerLaw2","[0]/pow(x,[1])");
      fPowerLaw2->SetParameters(1e-6,5e-2);
      //fPowerLaw2->SetParameters(1e14,-5,0);
      //fPowerLaw2->SetParameters(1e10,-8.8,40);
      //fPowerLaw2->SetParameters(1e10,-8.8,40);
      
      //fPowerLaw2->SetParLimits(0,1e-5,1e-2);
      //fPowerLaw2->SetParLimits(1,-5,0);
      //fPowerLaw2->SetParLimits(2,-10,0);
      hTempMC->Fit("fPowerLaw2","","",100,330);
      hTempMC->Fit("fPowerLaw2","","",100,330);
      hTempMC->Fit("fPowerLaw2","","",100,330);
      hTempMC->Fit("fPowerLaw2","","",100,330);
      hTempMC->Fit("fPowerLaw2","","",100,330);
      hTempMC->Fit("fPowerLaw2","","",100,330);
      TH1F *hFuncMC = (TH1F*)functionHist(fPowerLaw2,hTempMC,"Fit function Pt MC PP");
      hTempMC->SetXTitle("Jet p_{T} GeV/c");
      fPowerLaw2->SetMarkerStyle(8);
      fPowerLaw2->SetMarkerColor(4);
      fPowerLaw2->Draw("same");
      cPowerLawMC->cd(2);
      hRatioMC->SetAxisRange(100,330,"X");
      hRatioMC->Divide(hFuncMC);
      hRatioMC->SetTitle("Spectra to Fit Ratio");
      hRatioMC->SetXTitle("jet p_{T} GeV/c");
      hRatioMC->SetYTitle("mc/fit");
      hRatioMC->SetMarkerColor(4);
      hRatioMC->SetMarkerStyle(8);
      hRatioMC->Draw();
      cPowerLawMC->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PP_MC.gif",year,algoName[algo],nbins_cent,isFineBin),"RECREATE");
      
      // Reco
      TCanvas *cPowerLawReco = new TCanvas("cPowerLawReco","Reco",800,600);
      cPowerLawReco->Divide(2,1);
      cPowerLawReco->cd(1);
      cPowerLawReco->cd(1)->SetLogy();
      //TH1F *hTempReco = (TH1F*)uhist[i]->hReco->Clone("Unfolded Pt_spectra PP");
      TH1F *hTempReco = (TH1F*)fin->Get(Form("Unfolded_cent%d",i));
      hTempReco->SetTitle("Unfolded p_T cpectra PP");
      TH1F *hRatioReco = (TH1F*)hTempReco->Clone("hRatioReco");
      hTempReco->SetAxisRange(100,330,"X");
      hTempReco->SetMarkerStyle(22);
      hTempReco->SetMarkerColor(kBlue);
      TF1 *fPowerLaw3 = new TF1("fPowerLaw3","[0]*pow(x+[2],[1])");
      fPowerLaw3->SetParameters(1e14,-5,0);
      //fPowerLaw->SetParameters(1e10,-8.8,40);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      hTempReco->Fit("fPowerLaw3","LL","",100,330);
      TH1F *hFuncReco = (TH1F*)functionHist(fPowerLaw1,hTempReco,"Fit function Pt Reco PP");
      hTempReco->SetXTitle("Jet p_{T} GeV/c");
      fPowerLaw3->SetMarkerStyle(8);
      fPowerLaw3->SetMarkerColor(4);
      fPowerLaw3->Draw("same");
      cPowerLawReco->cd(2);
      hRatioReco->SetAxisRange(100,330,"X");
      hRatioReco->Divide(hFuncReco);
      hRatioReco->SetTitle("Spectra to Fit Ratio");
      hRatioReco->SetXTitle("jet p_{T} GeV/c");
      hRatioReco->SetYTitle("reco/fit");
      hRatioReco->SetMarkerColor(4);
      hRatioReco->SetMarkerStyle(8);
      hRatioReco->Draw();
      cPowerLawReco->SaveAs(Form("result-%d-%s-cent-%d-isFineBin-%d/PowerLaw_PP_Reco.gif",year,algoName[algo],nbins_cent,isFineBin),"RECREATE");
      
    }
  }
	

  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}
