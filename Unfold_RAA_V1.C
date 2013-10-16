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

#include "headers/utilities.h"
#include "headers/bayesianUnfold.h"
#include "headers/prior.h"

#include "RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"


using namespace std;

// jet RAA calculation - 9/23/2013

void Unfold_RAA_V1(int method = 1,int algo = 3,int isMC = 0){
  

#ifdef __CINT__
  gSystem->Load("RooUnfold-1.1.1/libRooUnfold.so");

  gStyle->SetErrorX(0.5);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.13);

#endif

  // start the timer
  TStopwatch timer;
  timer.Start();
 


  cout<<" ------------         Unfolding - Raghav 10 09 13          ----------           "<<endl;
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
  printf("AlgoPbPb:%s \n",algoName[algo]);

	
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

	
  //*******************lumi number for the sample all the units are in inverse mb ***************************
  float pbpblumi=150.; // this is for the 80 GeV HLT
  float pplumi = 5330; //2013 80 GeV HLT

  //*************************************************************************

  // Output file
	
  TFile *pbpb_Unfo = new TFile(Form("/net/hidsk0001/d00/scratch/rkunnawa/rootfiles/unfold_pbpb_cent1_pp_HLT_80_newvz_dup_removed_isMC_%d_2013_%s_%s.root", isMC , algoName[algo] ,chmet1),"RECREATE");	
  
  // Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];

  // Define RooUnfold response
  RooUnfoldResponse *response[nbins_cent+1];
	
  // Initialize Histograms
	
  for (int i=0;i<=nbins_cent;i++) {
    uhist[i] = new UnfoldingHistos(i);
    response[i] = new RooUnfoldResponse(uhist[i]->hResMeas,uhist[i]->hResTrue);
  }

  //load the MC files : First PbPb
  //TFile *pbpb_MC = new TFile("/net/hidsk0001/d00/scratch/rkunnawa/rootfiles/PbPb_2011_MC_276_pt.root","RECREATE");
  //TH1F *hpbpb_MC = new TH1F("hpbpb_MC","pt spectra from different pt hat files pbpb 2011 2.76 TeV JEC",nbins_rec,boundaries_rec);
  // Pthat binning
  //const int nbins_pthat = 8;
  const int nbins_pthat = 9;
  Double_t boundaries_pthat[nbins_pthat+1];
  char *fileName_pthat[nbins_pthat+1];
  Double_t xsection[nbins_pthat+1];
	
  ////// New MC samples

  boundaries_pthat[0]=30;
  fileName_pthat[0]="/mnt/hadoop/cms/store/user/kjung/PbPbMCProd/HydjetDrum_Dijet30_Tracking/merged/Dijet30_HydjetDrum_v27_Full_mergedV1.root";
  xsection[0]= 1.079e-02;
	
  boundaries_pthat[1]=50;
  fileName_pthat[1]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet50_HydjetDrum_v27_mergedV1.root";
  xsection[1]= 1.021e-03;
	
  boundaries_pthat[2]=80;
  fileName_pthat[2]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet80_HydjetDrum_v27_mergedV1.root";
  xsection[2]= 9.913e-05;
	
  boundaries_pthat[3]=100;
  fileName_pthat[3]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet100_HydjetDrum_v27_mergedV1.root ";
  xsection[3]= 3.069e-05 ;
	
  boundaries_pthat[4]=120;
  fileName_pthat[4]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/Dijet120_HydjetDrum_v28_mergedV1.root";
  xsection[4]= 1.128e-05;
	
  boundaries_pthat[5]=170;
  fileName_pthat[5]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet170_HydjetDrum_v27_mergedV1.root";
  xsection[5]= 1.470e-06;
	
  boundaries_pthat[6]=200;
  fileName_pthat[6]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/Dijet200_HydjetDrum_v28_mergedV1.root";
  xsection[6]= 5.310e-07;
	
  boundaries_pthat[7]=250;
  fileName_pthat[7]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/Dijet250_HydjetDrum_v28_mergedV1.root";
  xsection[7]= 1.192e-7;
	
  boundaries_pthat[8]=300;
  fileName_pthat[8]="/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/Dijet300_HydjetDrum_v28_mergedV1.root";
  xsection[8]= 3.176e-08;
	
  xsection[9] = 0;
  boundaries_pthat[9]=1000;

  // load the info from these files: 

  JetData *data[nbins_pthat];
  for (int i=0;i<nbins_pthat;i++) data[i] = new JetData(fileName_pthat[i],Form("%sJetAnalyzer/t",algoName[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));

  TH1F *hPtHat = new TH1F("hPtHat","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatRaw = new TH1F("hPtHatRaw","",nbins_pthat,boundaries_pthat);
	
  for (int i=0;i<nbins_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat,boundaries_pthat);
    data[i]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRaw->Add(hPtHatTmp);
    delete hPtHatTmp;
  }

  // Reweight to describe MC
  TF1 *fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);
  TF1 *fCentralityWeight = new TF1("fCentralityWeight","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
  fCentralityWeight->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);

	
	
  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
  TH1F *hCentData = new TH1F("hCentData","",40,0,40);
  TH1F *hCentMC = new TH1F("hCentMC","",40,0,40);
  
  TH1F *hVzData = new TH1F("hVzData","",60,-15,15);
  TH1F *hVzMC = new TH1F("hVzMC","",60,-15,15);

  hCent->Sumw2();
  hCentData->Sumw2();
  hCentMC->Sumw2();
  hVzData->Sumw2();
  hVzMC->Sumw2();
  

	
  // Fill PbPb MC
  for (int i=0;i<nbins_pthat;i++) {
    if (xsection[i]==0) continue;
    cout <<"Loading pthat"<<boundaries_pthat[i]
	 <<" sample, cross section = "<<xsection[i]
	 << Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[i],boundaries_pthat[i+1])<<endl;
    for (Long64_t jentry2=0; jentry2<data[i]->tJet->GetEntries();jentry2++) {
      data[i]->tEvt->GetEntry(jentry2);
      data[i]->tJet->GetEntry(jentry2);
      data[i]->tGenJet->GetEntry(jentry2);
      //if(data[i]->pthat<boundaries_pthat[i] || data[i]->pthat>boundaries_pthat[i+1]) continue;
      int pthatBin = hPtHat->FindBin(data[i]->pthat);
      float scale = (xsection[pthatBin-1]-xsection[pthatBin])/hPtHatRaw->GetBinContent(pthatBin);
      if(fabs(data[i]->vz)>15) continue;
      int cBin = hCent->FindBin(data[i]->bin)-1;
      double weight_cent=1;
      double weight_pt=1;
      double weight_vz=1;
	    
      weight_cent = fCentralityWeight->Eval(data[i]->bin);
      weight_vz = fVz->Eval(data[i]->vz);
      hCentMC->Fill(data[i]->bin,scale*weight_cent*weight_vz);
      hVzMC->Fill(data[i]->vz,scale*weight_cent*weight_vz);
      if (cBin>=nbins_cent) continue;
      if (cBin==-1) continue;
      hPtHat->Fill(data[i]->pthat,scale*weight_cent*weight_vz);
      /*
	int hasLeadingJet = 0;
	for (int k= 0; k < data[i]->njets; k++) {
	if ( data[i]->jteta[k]  > 2. || data[i]->jteta[k] < -2. ) continue;
	if ( data[i]->jtpt[k]>100) {
	hasLeadingJet = 1;
	}
	break;
	      
	}
	if (hasLeadingJet == 0) continue;
      */
      for (int k= 0; k < data[i]->njets; k++) {
	int subEvt=-1;
	if ( data[i]->refpt[k]  < 30. ) continue;
	if ( data[i]->jteta[k]  > 2. || data[i]->jteta[k] < -2. ) continue;
	if ( data[i]->trackMax[k]/data[i]->jtpt[k]<0.01) continue;
	for (int l= 0; l< data[i]->ngen;l++) {
	  if (data[i]->refpt[k]==data[i]->genpt[l]) {
	    subEvt = data[i]->gensubid[l];
	    break;
	  }
	}
	if (subEvt!=0) continue;
	//if (uhist[cBin]->hMeasMatch!=0) {
	//   int ptBinNumber = uhist[cBin]->hMeasMatch->FindBin(data[i]->jtpt[k]);
	//   int ratio = uhist[cBin]->hMeasMatch->GetBinContent(ptBinNumber);
	//if (ratio!=0) weight_pt = 1./ratio;
	//}
	if (!isMC||jentry2<data[i]->tJet->GetEntries()/2.) {
	  response[cBin]->Fill(data[i]->jtpt[k],data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	  uhist[cBin]-> hMatrix->Fill(data[i]->refpt[k],data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	  uhist[cBin]->hGen->Fill(data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	  //pbpb_MC->Fill(data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	}
	if (isMC&&jentry2>data[i]->tJet->GetEntries()/2.) {
	  uhist[cBin]-> hGen->Fill(data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	  uhist[cBin]-> hMeas->Fill(data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	  uhist[cBin]-> hMeasJECSys->Fill(data[i]->jtpt[k]*(1.+0.02/nbins_cent*(nbins_cent-i)),scale*weight_cent*weight_pt*weight_vz);
		
					
	}
      }
	    
    }
		
    
  }
  
  //hpbpb_MC->Write();
  //pbpb_MC->Close();

  cout<<endl<<endl<<endl<<"test PbPb mc loaded histograms"<<endl<<endl<<endl;

  uhist[nbins_cent-1]->hMatrix->Print("base");
  //uhist[nbins_cent-1]->hResMeas->Print("base");
  //uhist[nbins_cent-1]->hResTrue->Print("base");
  uhist[nbins_cent-1]->hGen->Print("base");

  cout<<endl<<endl<<"Finished loading PbPb mc"<<endl<<endl<<endl;
	
  uhist[nbins_cent-1]->hGen->Draw();

 // Fill PP MC
	
  //these files are from the HiForest Pythia signal twiki. 
  // ================ pp PtBin ======================================================================
  const int nbinsPP_pthat = 8;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
	
  boundariesPP_pthat[0]=30;
  fileNamePP_pthat[0]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet30_merged.root";
  xsectionPP[0]= 1.079e-02;
	
  boundariesPP_pthat[1]=50;
  fileNamePP_pthat[1]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet50_merged.root";
  xsectionPP[1]= 1.021e-03;
	
  boundariesPP_pthat[2]=80;
  fileNamePP_pthat[2]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet80_merged.root";
  xsectionPP[2]= 9.913e-05;
	
  boundariesPP_pthat[3]=120;
  fileNamePP_pthat[3]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet120_merged.root";
  xsectionPP[3]= 1.128e-05;
	
  boundariesPP_pthat[4]=170;
  fileNamePP_pthat[4]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet170_merged.root";
  xsectionPP[4]= 1.470e-06;
	
  boundariesPP_pthat[5]=200;
  fileNamePP_pthat[5]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet200_merged.root";
  xsectionPP[5]= 5.310e-07;
	
  boundariesPP_pthat[6]=250;
  fileNamePP_pthat[6]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet250_merged.root";
  xsectionPP[6]= 1.192e-07;
	
  boundariesPP_pthat[7]=300;
  fileNamePP_pthat[7]="/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged/pp276Dijet300_merged.root";
  xsectionPP[7]= 3.176e-08;
	
  xsectionPP[8] = 0;
  boundariesPP_pthat[8]=1000;

  /*

  // these files are from HiForestPA2013, pp signal sample prod 16 at 2.76 TeV 
  // ================ pp PtBin ======================================================================
  const int nbinsPP_pthat = 8;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
	
  //files changed by raghav to include the pp2013 data present /data/user/belt/hiforest2013/
  boundariesPP_pthat[0]=15;
  fileNamePP_pthat[0]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt15/HiForest_v77_merged01/pt15_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[0]= 2.034e-01;
	
  boundariesPP_pthat[1]=30;
  fileNamePP_pthat[1]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt30/HiForest_v77_merged01/pt30_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[1]= 1.075e-02;
	
  boundariesPP_pthat[2]=50;
  fileNamePP_pthat[2]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt50/HiForest_v77_merged01/pt50_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[2]= 1.025e-03;
	
  boundariesPP_pthat[3]=80;
  fileNamePP_pthat[3]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt80/HiForest_v77_merged01/pt80_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[3]= 9.865e-05;
	
  boundariesPP_pthat[4]=120;
  fileNamePP_pthat[4]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt120/HiForest_v77_merged01/pt120_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[4]= 1.129e-05;
	
  boundariesPP_pthat[5]=170;
  fileNamePP_pthat[5]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt170/HiForest_v77_merged01/pt170_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[5]= 1.465e-06;
	
  boundariesPP_pthat[6]=220;
  fileNamePP_pthat[6]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt220/HiForest_v77_merged01/pt220_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[6]= 2.837e-07;
	
  boundariesPP_pthat[7]=280;
  fileNamePP_pthat[7]="/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt280/HiForest_v77_merged01/pt280_P01_prod16_v77_merged_forest_0.root";
  xsectionPP[7]= 5.323e-08;
	
  xsectionPP[8] = 0;
  boundariesPP_pthat[8]=1000;

  */

		
  // Vertex reweighting for pp
 
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-15,15);// these two are taken from pawan's slides for the pp 40-60 file. 
  fVzPP->SetParameters(8.38193e-01,-0.00376098,0.00262389,-5.44672e-05,2.31855e-05);

  JetData *dataPP[nbinsPP_pthat];
  for (int i=0;i<nbinsPP_pthat;i++) dataPP[i] = new JetData(fileNamePP_pthat[i],Form("%sJetAnalyzer/t",algoNamePP[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));

  TH1F *hPtHatPP = new TH1F("hPtHatPP","",nbinsPP_pthat,boundariesPP_pthat);
  TH1F *hPtHatRawPP = new TH1F("hPtHatRawPP","",nbinsPP_pthat,boundariesPP_pthat);
  for (int i=0;i<nbinsPP_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
    dataPP[i]->tJet->Project("hPtHatTmp","pthat","abs(vz)<15");
    hPtHatRawPP->Add(hPtHatTmp);
    delete hPtHatTmp;
  }

  TH1F *hVzPPMC = new TH1F("hVzPPMC","",60,-15,15);
  hVzPPMC->Sumw2();
	
	
  for (int i=0;i<nbinsPP_pthat;i++) {
    if (xsectionPP[i]==0) continue;
    //float scale=(xsectionPP[i]-xsectionPP[i+1])/dataPP[i]->tJet->GetEntries(Form("pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1]));
    cout <<"Loading PP pthat"<<boundariesPP_pthat[i]
	 <<" sample, cross section = "<<xsectionPP[i]
	 << Form(" pthat>%.0f&&pthat<%.0f",boundariesPP_pthat[i],boundariesPP_pthat[i+1])<<endl;
    for (Long64_t jentry2=0; jentry2<dataPP[i]->tJet->GetEntries();jentry2++) {
      dataPP[i]->tEvt->GetEntry(jentry2);
      dataPP[i]->tJet->GetEntry(jentry2);
      dataPP[i]->tGenJet->GetEntry(jentry2);
      if(dataPP[i]->pthat<boundariesPP_pthat[i] || dataPP[i]->pthat>boundariesPP_pthat[i+1]) continue;
      if(dataPP[i]->bin<=28) continue;
      int pthatBin = hPtHatPP->FindBin(dataPP[i]->pthat);
      float scale = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/hPtHatRawPP->GetBinContent(pthatBin);
      if(fabs(dataPP[i]->vz)>15) continue;
      double weight_cent=1;
      double weight_pt=1;
      double weight_vz=1;
	    
      double newvz = dataPP[i]->vz + 0.32;//shift for pp2013 60-40 GeV
      weight_vz = 1./fVzPP->Eval(newvz);// this is including the new shift in the MC taken from pawan's slides 
      //if (weight_vz>5||weight_vz<0.5){
      //cout <<dataPP[i]->vz<<" "<<weight_vz<<endl;
      // weight_vz = 1;
      //}
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
	  response[nbins_cent]->Fill(dataPP[i]->jtpt[k],dataPP[i]->refpt[k],scale*weight_vz);
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



  TFile *infPbPb = TFile::Open("/net/hidsk0001/d00/scratch/rkunnawa/pbpb_jet80_pt_dup_removed.root");
  cout<<"loading PbPb Data"<<endl;
  uhist[nbins_cent-1]->hMeas = (TH1F*)infPbPb->Get("h");
 
  TFile *infPP = TFile::Open("/net/hidsk0001/d00/scratch/rkunnawa/pp_jet80_pt_dup_removed.root");
  cout<<"loading PP data"<<endl;
  uhist[nbins_cent]->hMeas = (TH1F*)infPP->Get("h");

  uhist[nbins_cent-1]->hMeas->Print("base");
  uhist[nbins_cent]->hMeas->Print("base");

  cout<<endl<<endl<<endl<<"finished loading Data"<<endl<<endl<<endl;

  //

  
  // starting the response matrix calculation
  cout <<"Response Matrix..."<<endl;
	
  TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",1200,800);
  cMatrix->Divide(3,3);
  TCanvas* cResponseNorm = new TCanvas("cResponseNorm","Normalized Response Matrix",1200,800);
  cResponseNorm->Divide(3,3);
	
  for (int i=0;i<=nbins_cent;i++){
    cMatrix->cd(i+1);
		
    TF1 *f = new TF1("f","[0]*pow(x+[2],[1])");
    f->SetParameters(1e10,-8.8,40);
    TH1F *hGenSpectraCorr = (TH1F*)uhist[i]->hMatrix->ProjectionX()->Clone(Form("hGenSpectraCorr_cent%d",i));
    hGenSpectraCorr->Fit("f"," ");
    hGenSpectraCorr->Fit("f","","");
    hGenSpectraCorr->Fit("f","LL ","");
    TH1F *fHist = functionHist(f,hGenSpectraCorr,Form("fHist_cent%d",i));
    hGenSpectraCorr->Divide(fHist);
    for (int y=1;y<=uhist[i]->hMatrix->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	if (uhist[i]->hMatrix->GetBinContent(x,y)<=0*uhist[i]->hMatrix->GetBinError(x,y)) {
	  uhist[i]->hMatrix->SetBinContent(x,y,0);
	  uhist[i]->hMatrix->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hMatrix->GetBinContent(x,y);
      }
	    
      for (int x=1;x<=uhist[i]->hMatrix->GetNbinsX();x++) {
	double ratio = 1;
	if (hGenSpectraCorr->GetBinContent(x)!=0) ratio = 1e5/hGenSpectraCorr->GetBinContent(x);
	uhist[i]->hMatrix->SetBinContent(x,y,uhist[i]->hMatrix->GetBinContent(x,y)*ratio);
	uhist[i]->hMatrix->SetBinError(x,y,uhist[i]->hMatrix->GetBinError(x,y)*ratio);
      }
    }
    //uhist[i]->hMatrix->Smooth(0);
		
    uhist[i]->hResponse = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponse_cent%d",i));
    for (int y=1;y<=uhist[i]->hResponse->GetNbinsY();y++) {
      double sum=0;
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {
	if (uhist[i]->hResponse->GetBinContent(x,y)<=0*uhist[i]->hResponse->GetBinError(x,y)) {
	  uhist[i]->hResponse->SetBinContent(x,y,0);
	  uhist[i]->hResponse->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponse->GetBinContent(x,y);
      }
	    
      for (int x=1;x<=uhist[i]->hResponse->GetNbinsX();x++) {
	if (sum==0) continue;
	double ratio = uhist[i]->hMeas->GetBinContent(y)/sum;
	if (uhist[i]->hMeas->GetBinContent(y)==0) ratio = 1e-100/sum;
	//uhist[i]->hResponse->SetBinContent(x,y,uhist[i]->hResponse->GetBinContent(x,y)*ratio);
	//uhist[i]->hResponse->SetBinError(x,y,uhist[i]->hResponse->GetBinError(x,y)*ratio);
      }
    }
		
    uhist[i]->hResponseNorm = (TH2F*)uhist[i]->hMatrix->Clone(Form("hResponseNorm_cent%d",i));
    for (int x=1;x<=uhist[i]->hResponseNorm->GetNbinsX();x++) {
      double sum=0;
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {
	if (uhist[i]->hResponseNorm->GetBinContent(x,y)<=0*uhist[i]->hResponseNorm->GetBinError(x,y)) {
	  uhist[i]->hResponseNorm->SetBinContent(x,y,0);
	  uhist[i]->hResponseNorm->SetBinError(x,y,0);
	}
	sum+=uhist[i]->hResponseNorm->GetBinContent(x,y);
      }
	    
      for (int y=1;y<=uhist[i]->hResponseNorm->GetNbinsY();y++) {
	if (sum==0) continue;
	double ratio = 1./sum;
	uhist[i]->hResponseNorm->SetBinContent(x,y,uhist[i]->hResponseNorm->GetBinContent(x,y)*ratio);
	uhist[i]->hResponseNorm->SetBinError(x,y,uhist[i]->hResponseNorm->GetBinError(x,y)*ratio);
      }
    }
		
		
    uhist[i]->hResponse->Draw("col");
    uhist[i]->hResponse->Draw("colz");
    uhist[i]->hResponse->SetAxisRange(1e-10,1,"Z");
	  
    cResponseNorm->cd(i+1);
    uhist[i]->hResponseNorm->Draw("colz");
    uhist[i]->hResponseNorm->SetAxisRange(1e-10,1,"Z");
  }

  cout<<"checking for response matrix PbPb"<<endl;
  uhist[nbins_cent-1]->hResponse->Print("base");
  uhist[nbins_cent]->hResponse->Print("base");

  cout<<"checking for normalized response matrix PbPb"<<endl;
  uhist[nbins_cent-1]->hResponseNorm->Print("base");
  cout<<"checking for normalized response matrix PP"<<endl;
  uhist[nbins_cent]->hResponseNorm->Print("base");

  


  // Unfolding starts: 

	
  cout << "==================================== TEST =====================================" << endl;
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  char chmet[100];
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
  TCanvas * cPbPb = new TCanvas("cPbPb","PbPb",1200,800);
  cPbPb->Divide(2,1);
	
  for (int i=0;i<=nbins_cent;i++) {
    cPbPb->cd(i+1)->SetLogy();
	  
	
    // Do Bin-by-bin
    cout<<"doing bin by bin unfolding"<<endl;

    TH1F* hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY();
    TH1F* hMCGen          = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",100,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));

    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
    //      uhist[i]->hRecoBinByBin = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
		


    // Do unfolding
    prior myPrior(uhist[i]->hMatrix,uhist[i]->hMeas,0);
    myPrior.unfold(uhist[i]->hMeas,1);
    TH1F *hPrior = (TH1F*)hMCGen->Clone("hPrior");
    removeZero(hPrior);
    //hPrior->Scale(uhist[i]->hMeas->Integral(0,1000)/hPrior->Integral(0,1000));
    
    TH1F *hReweighted = (TH1F*)(TH1F*)uhist[i]->hResponse->ProjectionY(Form("hReweighted_cent%d",i));

    //bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrix,hPrior,0);
    //myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    //bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrix,hPrior,0);
    //myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);

    bayesianUnfold myUnfolding(uhist[i]->hMatrix,myPrior.hPrior,0);
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter);

    cout <<"Unfolding bin "<<i<<endl;


    delete hBinByBinCorRaw;
    delete hMCGen;

    // Iteration Systematics
    for (int j=2;j<7;j++){
      bayesianUnfold myUnfoldingSys(uhist[i]->hMatrix,hPrior,0);
      myUnfoldingSys.unfold(uhist[i]->hMeas,j);
      uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
    }

    uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    // uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJECSys_cent%i",i));
    //uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    //uhist[i]->hRecoBinByBin = (TH1F*) unfold2.Hreco();
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));	
    
    //cleanup(uhist[i]->hReco);
    uhist[i]->hMeas->SetMarkerStyle(20);
    uhist[i]->hMeas->SetMarkerColor(2);
    uhist[i]->hReco->SetMarkerStyle(25);
    uhist[i]->hReco->Draw("");
    uhist[i]->hReco->SetAxisRange(100,300);
    //      TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
    //      hReproduced->SetMarkerColor(4);
    //      hReproduced->SetMarkerStyle(24);
    uhist[i]->hMeas->Draw("same");
	

  }


  cout<<"check the output histograms of unfolding"<<endl<<endl<<endl;
  cout<<"PbPb unfolding"<<endl;
  uhist[nbins_cent-1]->hReco->Print("base");
  cout<<"PP unfolding"<<endl;
  uhist[nbins_cent]->hReco->Print("base");



  cout<<"Now Proceeding to Calculation of RAA"<<endl;

  TCanvas * cRAA = new TCanvas("cRAA","RAA",1200,800);

	
  //necessary histograms to calculate RAA

  TH1F *hRecoBinByBinPP = (TH1F*)uhist[nbins_cent]->hRecoBinByBin->Clone("hRecoBinByBinPP");
  TH1F *hRecoPP          = (TH1F*)uhist[nbins_cent]->hReco->Clone("hRecoPP");
  TH1F *hMeasPP          = (TH1F*)uhist[nbins_cent]->hMeas->Clone("hMeasPP");
  //histos for measured pp data and mc, unfolded pp data, 
  TH1F *hMCGenPP = (TH1F*)uhist[nbins_cent]->hGen->Clone("hMCGenPP");
  TH1F *hDataMeasPP = (TH1F*)uhist[nbins_cent]->hMeas->Clone("hDataMeasPP");
  TH1F *hDataUnfoPP = (TH1F*)uhist[nbins_cent]->hReco->Clone("hDataUnfoPP");

  //divideBinWidth(hRecoBinByBinPP);
  //divideBinWidth(hRecoPP);
  //divideBinWidth(hMeasPP);
  
  hRecoBinByBinPP->Sumw2();
  hRecoPP->Sumw2();
  hMeasPP->Sumw2();
  hMCGenPP->Sumw2();
  hDataMeasPP->Sumw2();
  hDataUnfoPP->Sumw2();
  
  divideBinWidth(hMCGenPP);
  divideBinWidth(hDataMeasPP);
  divideBinWidth(hDataUnfoPP);

  // Scale PP histograms - 1./pplumi - lumiosity /64 - sigma inelastic /1000000 / 0.989 percentage of events which have |vz|<15. for this dup events removed it 1. 

  hRecoPP                ->Scale(1./pplumi/1000000);
  hRecoBinByBinPP       ->Scale(1./pplumi/1000000);
  hMeasPP                ->Scale(1./pplumi/1000000);
	
  //float CorFac[6] = {1.0331,1.0331,1.0300,1.0259,1.0217,1.0114};
	

  TH1F *hRecoBinByBinRAA = (TH1F*)uhist[nbins_cent-1]->hRecoBinByBin->Clone("hRecoBinByBinRAA");
  TH1F *hRecoRAA         = (TH1F*)uhist[nbins_cent-1]->hReco->Clone("hRecoRAA");
  TH1F *hMeasRAA         = (TH1F*)uhist[nbins_cent-1]->hMeas->Clone("hMeasRAA");
  TH1F *hMCGenPbPb = (TH1F*)uhist[nbins_cent-1]->hGen->Clone("hMCGenPbPb");
  TH1F *hDataMeasPbPb = (TH1F*)uhist[nbins_cent-1]->hMeas->Clone("hDataMeasPbPb");
  TH1F *hDataUnfoPbPb = (TH1F*)uhist[nbins_cent-1]->hReco->Clone("hDataUnfoPbPb");
  
	 
  //divideBinWidth(hRecoBinByBinRAA);
  //divideBinWidth(hRecoRAA);
  //divideBinWidth(hMeasRAA);
  
  hRecoBinByBinRAA->Sumw2();
  hRecoRAA->Sumw2();
  hMeasRAA->Sumw2();
  hMCGenPbPb->Sumw2();
  hDataMeasPbPb->Sumw2();
  hDataUnfoPbPb->Sumw2();
  
  divideBinWidth(hMCGenPbPb);
  divideBinWidth(hDataMeasPbPb);
  divideBinWidth(hDataUnfoPbPb);


  // Scale PbPb Hisotograms  and Calculate RAA, changed 8/19/2013 due to yaxian's requests 1./ correction factor (from MC)/7.65 - dont know what that is / 5.66 taa /40 = boundaries_cent[i+1] - boundaries_cent[i] /scaling for taa/ 1147500000 - min bias events, 0.97129 - ratio of events (abs vz<15)/events. for the dup events removed we have ratio of |vz|<15 is 0.94114

  //hRecoBinByBinRAA    ->Scale(1./0.94114/5.67/0.025/40/1147500000);
  //hRecoRAA             ->Scale(1./0.94114/5.67/0.025/40/1147500000);
  //hMeasRAA             ->Scale(1./0.94114/5.67/0.025/40/1147500000);
  //for merged it is as follows 
  //hRecoBinByBinRAA    ->Scale(1./5.67/0.025/40/1147500000);
  //hRecoRAA             ->Scale(1./5.67/0.025/40/1147500000);
  //hMeasRAA             ->Scale(1./5.67/0.025/40/1147500000);


  hRecoBinByBinRAA->Scale(1./(362.24*1147500000));
  hRecoRAA->Scale(1./(362.24*1147500000));
  hMeasRAA->Scale(1./(362.24*1147500000)); // ncoll * total min bias events. 

  hRecoRAA->Divide(hRecoPP);
  hRecoBinByBinRAA->Divide(hRecoBinByBinPP);
  hMeasRAA->Divide(hMeasPP);


  hRecoRAA->SetYTitle("Jet R_{AA}");  
  hRecoRAA->SetXTitle("Jet p_{T} GeV/c");
  hRecoRAA->SetTitle("");

  hRecoRAA->SetAxisRange(50,300,"X");
  hRecoRAA->SetAxisRange(0,2,"Y");
  hMeasRAA ->SetLineColor(kGreen);
  hMeasRAA ->SetMarkerColor(kGreen);
  hRecoRAA->SetMarkerStyle(24);
  hRecoRAA->SetLineColor(4);
  hRecoRAA->SetMarkerColor(4);
  hRecoBinByBinRAA->SetMarkerStyle(33);
  hRecoBinByBinRAA->SetLineColor(kRed);
  hRecoBinByBinRAA->SetMarkerColor(kRed); 

  cout <<" Plotting different methods RAA in akPu3PF " << endl;

  hRecoRAA->Draw("");
  hRecoRAA->Print("base");
  hMeasRAA->Draw("same");
  hMeasRAA->Print("base");
  hRecoBinByBinRAA->Draw("same");

  TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
  leg->SetTextSize(0.05);
  leg->AddEntry(hRecoRAA,"Bayesian","pl");
  leg->AddEntry(hRecoBinByBinRAA,"Bin-by-bin","pl");
  leg->AddEntry(hMeasRAA,"No unfolding","pl");
  leg->Draw();
  putCMSPrel(0.2,0.83,0.06);
  drawText("HLT_80 Anti-k_{T}Pu3PF,|eta|<2,|vz|<15 0-100%",0.2,0.23,20);
  
  cRAA->SaveAs("raa_calculation_HLT_80_pbpb_pp2013_newvz_dup_removed.pdf","RECREATE");
  cRAA->SaveAs("raa_calculation_HLT_80_pbpb_pp2013_newvz_dup_removed.root","RECREATE");
  cRAA->SaveAs("raa_calculation_HLT_80_pbpb_pp2013_newvz_dup_removed.C","RECREATE");


  //hRecoRAA->SetAxisRange(30,500,"Y");


  pbpb_Unfo->Write();
  pbpb_Unfo->Close();
  
  TFile *foutput = new TFile("RAA_PbPb2011_PP2013_new_vz_HLT_80_dup_removed.root","RECREATE");
  hRecoRAA->Write();
  hMeasRAA->Write();
  hRecoBinByBinRAA->Write();
  hMCGenPP->Write();
  hDataMeasPP->Write();
  hDataUnfoPP->Write();
  hMCGenPbPb->Write();
  hDataMeasPbPb->Write();
  hDataUnfoPbPb->Write();  

  TCanvas *c1 = new TCanvas("c1","",800,600);

  c1->SetLogy();

  hDataMeasPbPb->SetTitle("Measured vs Unfolded PbPb HLT_80 with dup events removed");
  hDataMeasPbPb->SetYTitle("counts/binwidth");
  hDataMeasPbPb->SetXTitle("Jet p_{T}");
  hDataMeasPbPb->SetMarkerStyle(24);
  hDataMeasPbPb->SetMarkerColor(kBlue);

  hDataUnfoPbPb->SetMarkerStyle(33);
  hDataUnfoPbPb->SetMarkerColor(kRed);
  //hDataUnfoPbPb->SetAxisRange(30,500,"Y");
  //hDataMeasPbPb->SetAxisRange(30,500,"Y");

  hDataMeasPbPb->Draw();
  hDataUnfoPbPb->Draw("same");


  TLegend *leg1 = myLegend(0.6,0.65,0.95,0.9);
  leg1->SetTextSize(0.05);
  leg1->AddEntry(hDataUnfoPbPb,"Bayesian","pl");
  leg1->AddEntry(hDataMeasPbPb,"No unfolding","pl");
  leg1->Draw();
  putCMSPrel(0.2,0.83,0.06);
  drawText("HLT_80 Anti-k_{T}Pu3PF,|eta|<2,|vz|<15 0-100%",0.2,0.23,20);

  c1->SaveAs("PbPb2011_HLT_80_measured_vs_unfold.pdf","RECREATE");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogy();

  hDataMeasPP->SetTitle("Measured vs Unfolded PP HLT_80 with dup events removed");
  hDataMeasPP->SetYTitle("counts/binwidth");
  hDataMeasPP->SetXTitle("Jet p_{T}");
  hDataMeasPP->SetMarkerStyle(24);
  hDataMeasPP->SetMarkerColor(kBlue);

  hDataUnfoPP->SetMarkerStyle(33);
  hDataUnfoPP->SetMarkerColor(kRed);
  //hDataUnfoPP->SetAxisRange(30,500,"Y");
  //hDataMeasPP->SetAxisRange(30,500,"Y");
  hDataMeasPP->Draw();
  hDataUnfoPP->Draw("same");


  TLegend *leg2 = myLegend(0.6,0.65,0.95,0.9);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(hDataUnfoPP,"Bayesian","pl");
  leg2->AddEntry(hDataMeasPP,"No unfolding","pl");
  leg2->Draw();
  putCMSPrel(0.2,0.83,0.06);
  drawText("HLT_80 Anti-k_{T}3PF,|eta|<2,|vz|<15",0.2,0.23,20);  

  c2->SaveAs("PP2013_HLT_80_measured_vs_unfold.pdf","RECREATE");

  foutput->Close();
	
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;


}

