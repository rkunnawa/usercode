// macro to calculate the RAA with systematics and excluding the Duplicate events. 

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
#include <cstdlib>
#include <cmath>

//#include "RooUnfold-1.1.1/src/RooUnfoldResponse.h"
//#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"

//#include "RooUnfold-1.1.1/src/RooUnfoldSvd.h"
//#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#include "utilities.h"
//#include "headers/bayesianUnfold.h"
//#include "headers/prior.h"


using namespace std;

// jet RpPb calculation - 9/23/2013

void Unfold_RpPb_V2(int method,int algo,int isMC = 0){
  

  //#ifdef __CINT__
  //gSystem->Load("RooUnfold-1.1.1/libRooUnfold");
  //#endif

  // start the timer
  TStopwatch timer;
  timer.Start();
  
  gStyle->SetErrorX(0.5);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.13);


  cout<<" ------------         Unfolding - Raghav 09 23 13          ----------           "<<endl;
  cout<<" ==============================================================================="<<endl;
	
  int nBayesianIter = 4;
  char chmet1[100];
  if(method==1) {
    sprintf(chmet1,"Bayes_unfo");
  } else if(method==2) {
    sprintf(chmet1,"Svd_unfo ");
  } else if(method==3) {
    sprintf(chmet1,"BinByBin_unfo");
  }
	
  printf("Method : %s \n",chmet1);
  printf("AlgoPP : %s \n",algoNamePP[algo]);
	
  cout << "==================================== TRAIN ====================================" << endl;
	
  // ================ Bin Size ======================================================================
	
	
  // ================ PbPb PtBin ======================================================================
	
	
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

	
  //*******************lumi number for the sample all the units are in inverse mb ***************************
  float ppblumi=30.9; // this is for the 80 GeV HLT
  float pplumi = 5330; //2013 80 GeV HLT

  //*************************************************************************

  // Output file
	
  TFile *ppb_Unfo = new TFile("RpPb_calculation_data_vs_mc.root","RECREATE");	
  // Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];// Should i make 3 copies of them

  // Define RooUnfold response
  //RooUnfoldResponse *response[nbins_cent+1];
  //RooUnfoldResponse res(uhist[0]->hResMeas,uhist[0]->hResTrue);

	
  // Initialize Histograms
	
  for (int i=0;i<=nbins_cent;i++) {
    uhist[i] = new UnfoldingHistos(i);
    //response[i] = new RooUnfoldResponse(uhist[i]->hResMeas,uhist[i]->hResTrue);
  }

  //load the MC files : First PP

	
  // load the required pp signal sample 5.02 GeV files as a cross check for pPb.
	
  const int nbins_pthat_pptrack = 8;
  Double_t boundaries_pthat_pptrack[nbins_pthat_pptrack+1];
  char *fileName_pthat_pptrack[nbins_pthat_pptrack+1];
  Double_t xsection_pptrack[nbins_pthat_pptrack+1];
	
  boundaries_pthat_pptrack[0] = 15;
  fileName_pthat_pptrack[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt15/HiForest_v77_v2_merged01/pt15_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[0] = 5.335e-01;
	
  boundaries_pthat_pptrack[1] = 30;
  fileName_pthat_pptrack[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt30/HiForest_v77_v2_merged01/pt30_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[1] = 3.378e-02;
	
  boundaries_pthat_pptrack[2] = 50;
  fileName_pthat_pptrack[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt50/HiForest_v77_v2_merged01/pt50_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[2] = 3.778e-03;
	
  boundaries_pthat_pptrack[3] = 80;
  fileName_pthat_pptrack[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt80/HiForest_v77_v2_merged01/pt80_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[3] = 4.412e-04;
	
  boundaries_pthat_pptrack[4] = 120;
  fileName_pthat_pptrack[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt120/HiForest_v77_v2_merged01/pt120_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[4] = 6.147e-05;
	
  boundaries_pthat_pptrack[5] = 170;
  fileName_pthat_pptrack[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt170/HiForest_v77_v2_merged01/pt170_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[5] = 1.018e-05;
	
  boundaries_pthat_pptrack[6] = 220;
  fileName_pthat_pptrack[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt220/HiForest_v77_v2_merged02/pt220_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[6] = 2.477e-06;
	
  boundaries_pthat_pptrack[7] = 280;
  fileName_pthat_pptrack[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt280/HiForest_v77_v2_merged01/pt280_HP04_hiforest77_hiSignal.root";
  xsection_pptrack[7] = 6.160e-07;
	
  boundaries_pthat_pptrack[8] = 1000;
  xsection_pptrack[8] = 0;
	

		
  // Vertex reweighting for pp
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-15,15);
  fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
  //TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-15,15);// these two are taken from pawan's slides for the pp 40-60 file. 
  //fVzPP->SetParameters(8.38193e-01,-0.00376098,0.00262389,-5.44672e-05,2.31855e-05);

  JetData *dataPP[nbins_pthat_pptrack];
  for (int i=0;i<nbins_pthat_pptrack;i++) dataPP[i] = new JetData(fileName_pthat_pptrack[i],Form("%sJetAnalyzer/t",algoNamePP[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));

  TH1F *hPtHatPP = new TH1F("hPtHatPP","",nbins_pthat_pptrack,boundaries_pthat_pptrack);
  TH1F *hPtHatRawPP = new TH1F("hPtHatRawPP","",nbins_pthat_pptrack,boundaries_pthat_pptrack);
  for (int i=0;i<nbins_pthat_pptrack;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat_pptrack,boundaries_pthat_pptrack);
    dataPP[i]->tJet->Project("hPtHatTmp","pthat","abs(vz)<15&&abs(jteta)<1");
    hPtHatRawPP->Add(hPtHatTmp);
    delete hPtHatTmp;
  }

  TH1F *hVzPPMC = new TH1F("hVzPPMC","",60,-15,15);
  hVzPPMC->Sumw2();
	
	
  for (int i=0;i<nbins_pthat_pptrack;i++) {
    if (xsection_pptrack[i]==0) continue;
    //float scale=(xsectionPP[i]-xsectionPP[i+1])/dataPP[i]->tJet->GetEntries(Form("pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1]));
    cout <<"Loading PP pthat"<<boundaries_pthat_pptrack[i]
	 <<" sample, cross section = "<<xsection_pptrack[i]
	 << Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat_pptrack[i],boundaries_pthat_pptrack[i+1])<<endl;
    for (Long64_t jentry2=0; jentry2<dataPP[i]->tJet->GetEntries();jentry2++) {
      dataPP[i]->tEvt->GetEntry(jentry2);
      dataPP[i]->tJet->GetEntry(jentry2);
      dataPP[i]->tGenJet->GetEntry(jentry2);
      if(dataPP[i]->pthat<boundaries_pthat_pptrack[i] || dataPP[i]->pthat>boundaries_pthat_pptrack[i+1]) continue;
      if(dataPP[i]->bin<=28) continue;
      int pthatBin = hPtHatPP->FindBin(dataPP[i]->pthat);
      float scale = (xsection_pptrack[pthatBin-1]-xsection_pptrack[pthatBin])/hPtHatRawPP->GetBinContent(pthatBin);
      //float scale = xsection_pptrack[pthatBin-1]-xsection_pptrack[pthatBin];
      if(fabs(dataPP[i]->vz)>15) continue;
      double weight_cent=1;
      double weight_pt=1;
      double weight_vz=1;
	    
      weight_vz = fVzPP->Eval(dataPP[i]->vz);
      //double newvz = dataPP[i]->vz + 0.32;//shift for pp2013 60-40 GeV
      //weight_vz = 1./fVzPP->Eval(newvz);// this is including the new shift in the MC taken from pawan's slides 
      if (weight_vz>5||weight_vz<0.5) {
      cout <<dataPP[i]->vz<<" "<<weight_vz<<endl;
      //weight_vz = 1;
      }
      hPtHatPP->Fill(dataPP[i]->pthat,scale*weight_vz);
      int hasLeadingJet = 0;
      hVzPPMC->Fill(dataPP[i]->vz,scale*weight_vz);
      /*
	for (int k= 0; k < dataPP[i]->njets; k++) {
	if ( dataPP[i]->jteta[k]  > 2. || dataPP[i]->jteta[k] < -2. ) continue;
	if ( dataPP[i]->jtpt[k]>100) {
	hasLeadingJet = 1;
	}
	break;
	      
	}
	if (hasLeadingJet == 0) continue;
      */
      for (int k= 0; k < dataPP[i]->njets; k++) {
	int subEvt=-1;
	if ( dataPP[i]->refpt[k]  < 30. ) continue;
	if ( dataPP[i]->jteta[k]  > 2. || dataPP[i]->jteta[k] < -2. ) continue;
	if ( dataPP[i]->trackMax[k]/dataPP[i]->jtpt[k]<0.01) continue;
	//if (uhist[nbins_cent]->hMeasMatch!=0) {
	//   int ptBinNumber = uhist[nbins_cent]->hMeasMatch->FindBin(dataPP[i]->jtpt[k]);
	//   int ratio = uhist[nbins_cent]->hMeasMatch->GetBinContent(ptBinNumber);
	//if (ratio!=0) weight_pt = 1./ratio;
	//}
	      
	if (!isMC||jentry2<dataPP[i]->tJet->GetEntries()/2.) {
	  //response[nbins_cent]->Fill(dataPP[i]->jtpt[k],dataPP[i]->refpt[k],scale*weight_vz);
	  uhist[nbins_cent]-> hMatrix->Fill(dataPP[i]->refpt[k],dataPP[i]->jtpt[k],scale*weight_vz);
	  uhist[nbins_cent]->hGen->Fill(dataPP[i]->jtpt[k],scale*weight_vz);
	}
	if (isMC&&jentry2>dataPP[i]->tJet->GetEntries()/2.) {
		
	  uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);
	  uhist[nbins_cent]-> hMeas->Fill(dataPP[i]->jtpt[k],scale*weight_vz);
	}
      }
    }
  }
  cout<<endl<<endl<<endl<<"test pp mc loaded histograms"<<endl<<endl<<endl;
  uhist[nbins_cent]->hMatrix->Print("base");
  //uhist[nbins_cent]->hResMeas->Print("base");
  //uhist[nbins_cent]->hResTrue->Print("base");
  uhist[nbins_cent]->hGen->Print("base");

  cout<<endl<<endl<<endl<<"finished loading PP mc"<<endl<<endl<<endl;


  //load the PbPb hyjing+pythia mc files for unfolding. 

  /* 	
  //these are the files required for the unfolding. 
  const int nbins_pthat_ppb = 9;
  Double_t boundaries_pthat_ppb[nbins_pthat_ppb+1];
  char *fileName_pthat_ppb[nbins_pthat_ppb+1];
  Double_t xsection_ppb[nbins_pthat_ppb+1];
	
  boundaries_pthat_ppb[0] = 15;
  fileName_pthat_ppb[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt15/HiForest_v77_merged01/pt15_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[0] = 5.335e-01;// make sure that you have the correct one here.
	
  boundaries_pthat_ppb[1] = 30;
  fileName_pthat_ppb[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt30/HiForest_v77_merged01/pt30_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[1] = 3.378e-02;
	
  boundaries_pthat_ppb[2] = 50;
  fileName_pthat_ppb[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt50/HiForest_v77_merged01/pt50_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[2] = 3.778e-03;
	
  boundaries_pthat_ppb[3] = 80;
  fileName_pthat_ppb[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt80/HiForest_v77_merged01/pt80_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[3] = 4.412e-04;
	
  boundaries_pthat_ppb[4] = 120;
  fileName_pthat_ppb[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt120/HiForest_v77_merged01/pt120_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[4] = 6.147e-05;
	
  boundaries_pthat_ppb[5] = 170;
  fileName_pthat_ppb[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt170/HiForest_v77_merged01/pt170_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[5] = 1.018e-05;
	
  boundaries_pthat_ppb[6] = 220;
  fileName_pthat_ppb[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt220/HiForest_v77_merged01/pt220_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[6] = 2.477e-06;
	
  boundaries_pthat_ppb[7] = 280;
  fileName_pthat_ppb[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt280/HiForest_v77_merged01/pt280_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[7] = 6.160e-07;
	
  boundaries_pthat_ppb[8] = 370;
  fileName_pthat_ppb[8] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt370/HiForest_v77_merged01/pt370_HP04_prod16_v77_merged_forest_0.root";
  xsection_ppb[8] = 1.088e-07;
	
  boundaries_pthat_ppb[8] = 1000;
  xsection_ppb[8] = 0;
	
  */

  // now load the data files 

  // ppb ak3pf

  TFile *infData;
  cout<<"loading PbPb data"<<endl;
  //infData = new TFile("pbpb_jet80_pt_dup_removed.root");
  //uhist[nbins_cent-1]->hMeas = (TH1F*)infData->Get("h");

  // for the merged file:
  //infData = new TFile("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/pPb/pPb_meas_pt_histos.root");
  infData = new TFile("merge_test_output.root");
  //for this file pPb_meas_pt.root the total number of events = 4825168, total no of events with selection cuts - 1265115
  
  TH1F *htemp = (TH1F*)infData->Get("hCombined");
  //uhist[nbins_cent-1]->hMeas = rebin(htemp,Form("hMeas_cent%d",nbins_cent-1));
  uhist[nbins_cent-1]->hMeasMatch = 0;
  Float_t nEntries = htemp->GetEntries();
  cout<<"finished Loading pPb data"<<endl;
  cout<<"test pPb data histogram"<<endl;
  uhist[nbins_cent-1]->hMeas->Print("base");

  /*	
  // adding the trigger efficiency here: according to Yen-Jie
  Float_t bincontent = 0., binerror = 0.;
  for(int i = 0;i<nbins_cent;i++){
    for(int j = 0;j<nbins_recrebin;j++){
      bincontent = uhist[i]->hMeas->GetBinContent(j);
      binerror = uhist[i]->hMeas->GetBinError(j);
      bincontent /= trigEffInc[j];
      binerror /= trigEffInc[j];
      uhist[i]->SetBinContent(bincontent);
      uhist[i]->SetBinError(binerror);
    }
  }
  */
  cout<<endl<<endl<<"test pPb data after the trigger efficiency applied"<<endl<<endl<<endl;
  uhist[nbins_cent-1]->hMeas->Print("base");


  cout<<"Now Proceeding to Calculation of RpPb"<<endl;

  TCanvas *c1;
  c1 = new TCanvas("c1","",800,600);
  c1->Divide(2,1);
  
	
  //necessary histograms to calculate RAA

  TH1F *hGenPP          = (TH1F*)uhist[nbins_cent]->hGen->Clone("hGenPP");
  Float_t nEntries_mc = hGenPP->GetEntries();
  TH1F *hMCGenPP = (TH1F*)uhist[nbins_cent]->hGen->Clone("hMCGenPP");

  divideBinWidth(hGenPP);
  divideBinWidth(hMCGenPP);

  hGenPP->Sumw2();
  hMCGenPP->Sumw2();
  
  hMCGenPP->Scale(1./(70*1000));

  c1->cd(1);
  hGenPP->Draw();

  TH1F *hMeasRAA         = (TH1F*)uhist[nbins_cent-1]->hMeas->Clone("hMeasRAA");
  TH1F *hDataMeasPPb = (TH1F*)uhist[nbins_cent-1]->hMeas->Clone("hDataMeasPPb");

  
  divideBinWidth(hMeasRAA);
  divideBinWidth(hDataMeasPPb);

  hMeasRAA->Sumw2();
  hDataMeasPPb->Sumw2();

  c1->cd(2);
  hMeasRAA->Draw();

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->Divide(2,1);

  //probably this pp gen scaling should have some values as well. 

  //hGenPP                ->Scale(1/nEntries_mc);
  c2->cd(1);
  hGenPP->Draw();

  //correct scaling factors for RpPb: nEntries, 0.85224 (ratio of entries which have abs(eta)<1 for ) and should ideally have sigma_inelastic 70mb but i dont know where to put it. 

  //hMeasRAA             ->Scale(1./6.9/0.85224/70/1000000);
  hMeasRAA             ->Scale(1./6.9);
  hDataMeasPPb         ->Scale(1./6.9);
  c2->cd(2);
  hDataMeasPPb->Draw();

  hMeasRAA->Divide(hGenPP);

  TCanvas * cRAA = new TCanvas("cRAA","RAA",1200,800);

  hMeasRAA->SetYTitle("Jet R_{pPb}");  
  hMeasRAA->SetXTitle("Jet p_{T} GeV/c");
  hMeasRAA->SetTitle("");

  hMeasRAA->SetAxisRange(30,300,"X");
  //hMeasRAA->SetAxisRange(0,5,"Y");
  hMeasRAA ->SetLineColor(kBlack);
  hMeasRAA ->SetMarkerColor(kBlack);


  cout <<" Plotting different methods RpPb in akPU3PF " << endl;


  hMeasRAA->Draw();
  hMeasRAA->Print("base");

  TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
  leg->SetTextSize(0.05);

  leg->AddEntry(hMeasRAA,"No unfolding","pl");
  leg->Draw();
  putCMSPrel(0.2,0.83,0.06);
  drawText("HLT_80_100 Anti-k_{T}PU3PF,|eta|<1,|vz|<15",0.2,0.23,20);
  
  cRAA->SaveAs("RpPb_meas_preliminary_data_vs_mc.pdf","RECREATE");
  cRAA->SaveAs("RpPb_meas_preliminary_data_vs_mc.root","RECREATE");
  cRAA->SaveAs("RpPb_meas_preliminary_data_vs_mc.C","RECREATE");

  ppb_Unfo->Write();
  ppb_Unfo->Close();
  
  TFile f("RpPb_preliminary_histos.root","RECREATE");
  hMeasRAA->Write();
  hMCGenPP->Write();
  hDataMeasPPb->Write();

  f.Close();
	
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;


}

