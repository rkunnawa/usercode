// RpPb calculation after Unfolding, 
// the PP reference in this one is with the MC.  we dont have any data at that energy 5.02 TeV 
// this is done without chMax/jtpt cut 

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

#include "headers/utilities_V0.h"
#include "headers/bayesianUnfold.h"
#include "headers/prior.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
// Update Raghav Kunnawalkam Elayavalli 10/29/2013
//==============================================================================

void Unfold_RpPb_V0(int method = 1,int algo = 3,bool useSpectraFromFile = 0, bool useMatrixFromFile = 0, int doToy = 0, int isMC = 0,char *spectraFileName = "", int doJECSys = 0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
{
  //#ifdef __CINT__
  //gSystem->Load("libRooUnfold");
  //#endif
  
  //start timer
  TStopwatch timer;
  timer.Start();
	
  int useFixedIterativeSys = 1;
  int isPyquen = 1;
	
  bool yinglu = 0;// 0 - file location in MIT
	
  gStyle->SetErrorX(0.5);
  gStyle->SetPaintTextFormat("3.2f");
  gStyle->SetOptLogz(1);
  gStyle->SetPadRightMargin(0.13);
  cout<<" ------------         Unfolding Raghav 10 17 13          ----------           "<<endl;
  cout<<" ==============================================================================="<<endl;
	
  int nBayesianIter = 4;// 4 or 6
  char chmet1[100];
  if(method==1) {
    sprintf(chmet1,"Bayes unfo");
  } else if(method==2) {
    sprintf(chmet1,"Svd unfo ");
  } else if(method==3) {
    sprintf(chmet1,"BinByBin unfo");
  }
	
  printf("Method : %s \n",chmet1);
	
  cout << "==================================== TRAIN ====================================" << endl;
	
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
	
  // Pthat binning
  //const int nbins_pthat = 8;
  const int nbins_pthat = 9;
  Double_t boundaries_pthat[nbins_pthat+1];
  char *fileName_pthat[nbins_pthat+1];
  Double_t xsection[nbins_pthat+1];
	
  char *fileName_pthat_pq;
  fileName_pthat_pq="/hadoop/store/user/belt/hiForest2/v27_28/pyquenFull80_HYDJET.root";
  //fileName_pthat_pq="/hadoop/store/user/belt/hiForest2/pyquen_hydjet_hiforest_test.root";
	
	
  ////// New MC sample
  //this MC is without the PU subtraction, but the data we are using is akPu3PF, 
  /*
  if (!yinglu){
		
    boundaries_pthat[0] = 15;
    fileName_pthat[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt15/HiForest_v77_merged01/pt15_HP04_prod16_v77_merged_forest_0.root";
    xsection[0] = 5.335e-01;// make sure that you have the correct one here.
		
    boundaries_pthat[1] = 30;
    fileName_pthat[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt30/HiForest_v77_merged01/pt30_HP04_prod16_v77_merged_forest_0.root";
    xsection[1] = 3.378e-02;
		
    boundaries_pthat[2] = 50;
    fileName_pthat[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt50/HiForest_v77_merged01/pt50_HP04_prod16_v77_merged_forest_0.root";
    xsection[2] = 3.778e-03;
		
    boundaries_pthat[3] = 80;
    fileName_pthat[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt80/HiForest_v77_merged01/pt80_HP04_prod16_v77_merged_forest_0.root";
    xsection[3] = 4.412e-04;
		
    boundaries_pthat[4] = 120;
    fileName_pthat[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt120/HiForest_v77_merged01/pt120_HP04_prod16_v77_merged_forest_0.root";
    xsection[4] = 6.147e-05;
		
    boundaries_pthat[5] = 170;
    fileName_pthat[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt170/HiForest_v77_merged01/pt170_HP04_prod16_v77_merged_forest_0.root";
    xsection[5] = 1.018e-05;
		
    boundaries_pthat[6] = 220;
    fileName_pthat[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt220/HiForest_v77_merged01/pt220_HP04_prod16_v77_merged_forest_0.root";
    xsection[6] = 2.477e-06;
		
    boundaries_pthat[7] = 280;
    fileName_pthat[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt280/HiForest_v77_merged01/pt280_HP04_prod16_v77_merged_forest_0.root";
    xsection[7] = 6.160e-07;
		
    boundaries_pthat[8] = 370;
    fileName_pthat[8] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt370/HiForest_v77_merged01/pt370_HP04_prod16_v77_merged_forest_0.root";
    xsection[8] = 1.088e-07;
		
    boundaries_pthat[9] = 1000;
    xsection[9] = 0;
		

  }
  */


  //MC JEC with PU subtraction algorithm

  if(!yinglu){
    boundaries_pthat[0] = 15;
    fileName_pthat[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt15_HP04_prod25_v85_merged_forest_0.root";
    xsection[0] = 5.335e-01;

    boundaries_pthat[1] = 30;
    fileName_pthat[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt30_HP04_prod25_v85_merged_forest_0.root";
    xsection[1] = 3.378e-02;

    boundaries_pthat[2] = 50;
    fileName_pthat[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt50_HP04_prod25_v85_merged_forest_0.root";
    xsection[2] = 3.778e-02;

    boundaries_pthat[3] = 80;
    fileName_pthat[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt80_HP04_prod25_v85_merged_forest_0.root";
    xsection[3] = 4.412e-04;

    boundaries_pthat[4] = 120;
    fileName_pthat[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt120_HP04_prod25_v85_merged_forest_0.root";
    xsection[4] = 6.147e-05;

    boundaries_pthat[5] = 170;
    fileName_pthat[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt170_HP04_prod25_v85_merged_forest_0.root";
    xsection[5] = 1.018e-05;

    boundaries_pthat[6] = 220;
    fileName_pthat[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt220_HP04_prod25_v85_merged_forest_0.root";
    xsection[6] = 2.477e-06;

    boundaries_pthat[7] = 280;
    fileName_pthat[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt280_HP04_prod25_v85_merged_forest_0.root";
    xsection[7] = 6.160e-07;

    boundaries_pthat[8] = 370;
    fileName_pthat[8] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod25/HiForest_v85_merged01/pt370_HP04_prod25_v85_merged_forest_0.root";
    xsection[8] = 1.088e-07;

    boundaries_pthat[9] = 1000;
    xsection[9] = 0;

  }
       
	
	
	
  // ================ pp PtBin ======================================================================
  const int nbinsPP_pthat = 8;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
	

  if (!yinglu) {
		
    boundariesPP_pthat[0] = 15;
    fileNamePP_pthat[0] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt15/HiForest_v77_v2_merged01/pt15_HP04_hiforest77_hiSignal.root";
    xsectionPP[0] = 5.335e-01;// make sure that you have the correct one here.
		
    boundariesPP_pthat[1] = 30;
    fileNamePP_pthat[1] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt30/HiForest_v77_v2_merged01/pt30_HP04_hiforest77_hiSignal.root";
    xsectionPP[1] = 3.378e-02;
		
    boundariesPP_pthat[2] = 50;
    fileNamePP_pthat[2] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt50/HiForest_v77_v2_merged01/pt50_HP04_hiforest77_hiSignal.root";
    xsectionPP[2] = 3.778e-03;
		
    boundariesPP_pthat[3] = 80;
    fileNamePP_pthat[3] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt80/HiForest_v77_v2_merged01/pt80_HP04_hiforest77_hiSignal.root";
    xsectionPP[3] = 4.412e-04;
		
    boundariesPP_pthat[4] = 120;
    fileNamePP_pthat[4] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt120/HiForest_v77_v2_merged01/pt120_HP04_hiforest77_hiSignal.root";
    xsectionPP[4] = 6.147e-05;
		
    boundariesPP_pthat[5] = 170;
    fileNamePP_pthat[5] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt170/HiForest_v77_v2_merged01/pt170_HP04_hiforest77_hiSignal.root";
    xsectionPP[5] = 1.018e-05;
		
    boundariesPP_pthat[6] = 220;
    fileNamePP_pthat[6] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt220/HiForest_v77_v2_merged01/pt220_HP04_hiforest77_hiSignal.root";
    xsectionPP[6] = 2.477e-06;
		
    boundariesPP_pthat[7] = 280;
    fileNamePP_pthat[7] = "/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt280/HiForest_v77_v2_merged01/pt280_HP04_hiforest77_hiSignal.root";
    xsectionPP[7] = 6.160e-07;
		
    boundariesPP_pthat[8] = 1000;
    xsectionPP[8] = 0;

  }
	
  //*******************lumi number for the sample***************************
  float lumi=30.9;//ppb
  float pplumi=5300.;
  //float lumi=129.;
  //float pplumi=212.;
  //*************************************************************************
	
  // Output file
  TFile *ppb_Unfo = new TFile(Form("result-2013-ppb-%s-cent-%d/ppb_merge_MB_eta_CM_1_lowest_pp_mc_Unfo_%s_cent_%d.root",algoName[algo],nbins_cent,algoName[algo],nbins_cent),"RECREATE");
  // Histograms used by RooUnfold
  UnfoldingHistos *uhist[nbins_cent+1];
	
  // Define RooUnfold response
  RooUnfoldResponse *response[nbins_cent+1];
	
  // Initialize Histograms
	
  for (int i=0;i<=nbins_cent;i++) {
    uhist[i] = new UnfoldingHistos(i);
    response[i] = new RooUnfoldResponse(uhist[i]->hResMeas,uhist[i]->hResTrue);
  }
	
	
  // Initialize reweighting functions
	
  // Vertex & centrality reweighting for PbPb
  //TF1 *fVz;
  //TF1* fCentralityWeight;
  TCut dataSelection;
  TCut dataSelectionPPb;
  TCut dataSelectionPP;
  TCut TriggerSelectionPP;
  TCut TriggerSelectionPPb;
	
  if (isMC) {
    // MC closure test, no reweighting
    //fVz = new TF1("fVz","1");
    //fCentralityWeight =  new TF1("fCentralityWeight","1");
    dataSelection = "abs(vz)<15&&trackMax/jtpt>0.01&&abs(jteta)<2";
  } else {
    // Reweight to describe data
    // 
    //fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
    //fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);
    //fCentralityWeight = new TF1("fCentralityWeight","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
    //fCentralityWeight->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);
    //dataSelectionPbPb = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2";
    //dataSelectionPP = "abs(vz)<15&&pPAcollisionEventSelectionPA&&pHBHENoiseFilter&&abs(jteta)<2";
    //TriggerSelectionPP = "HLT_PAJet80_NoJetID_v1";
    //TriggerSelectionPbPb ="HLT_HIJet80_v1";
		
    dataSelectionPPb = "abs(eta)<2";//not used here since we are taking it form the merging file 
    dataSelectionPP = "abs(eta)<2";
		
    TriggerSelectionPP = "jet80";
    TriggerSelectionPPb = "jet80";
		
  }
  
  // Vertex reweighting for pp //this one was the old one that i had
  //TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  //fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
	
  //this is from yaxian's file 
  TF1 * fVz = new TF1("fVz","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);  
  fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
  
   	
  // Read data file
  TFile *infData;
  infData = new TFile("merge_ppb_MB_eta_CM_1_lowest_HLT_V2.root");
  TH1F *htest = (TH1F*)infData->Get("hppbComb");
  uhist[nbins_cent-1]->hMeas = rebin2(htest,"hMeas_cent0");
  uhist[nbins_cent-1]->hMeas->Print("base");

	
  // Setup jet data branches
  JetData *data[nbins_pthat];
  JetData *dataPP[nbins_pthat];
  JetData *data_pq;
	
  for (int i=0;i<nbins_pthat;i++) data[i] = new JetData(fileName_pthat[i],Form("%sJetAnalyzer/t",algoName[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));
  for (int i=0;i<nbinsPP_pthat;i++) dataPP[i] = new JetData(fileNamePP_pthat[i],Form("%sJetAnalyzer/t",algoNamePP[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));
	
  if (yinglu) data_pq = new JetData(fileName_pthat_pq,Form("%sJetAnalyzer/t",algoName[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));	
	
  // Come back to the output file directory
  ppb_Unfo->cd();
	

  // Get Jet spectra from data file
  cout <<"Reading data..."<<endl;
	
  TCanvas * cInput = new TCanvas("cInput","Input",1200,800);
  cInput->Divide(3,3);

	
  TH1F *hCent = new TH1F("hCent","",nbins_cent,boundaries_cent);
  TH1F *hCentData = new TH1F("hCentData","",40,0,40);
  TH1F *hCentMC = new TH1F("hCentMC","",40,0,40);
	
  TH1F *hVzData = new TH1F("hVzData","",60,-15,15);
  TH1F *hVzMC = new TH1F("hVzMC","",60,-15,15);
  TH1F *hVzPPData = new TH1F("hVzPPData","",60,-15,15);
  TH1F *hVzPPMC = new TH1F("hVzPPMC","",60,-15,15);
  hCent->Sumw2();
  hCentData->Sumw2();
  hCentMC->Sumw2();
  hVzData->Sumw2();
  hVzMC->Sumw2();
  hVzPPData->Sumw2();
  hVzPPMC->Sumw2();
	
  
  TCanvas *c = new TCanvas("c","",600,600);
  TH1F *hPtHatPPb = new TH1F("hPtHatPPb","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatRawPPb = new TH1F("hPtHatRawPPb","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatPP = new TH1F("hPtHatPP","",nbinsPP_pthat,boundariesPP_pthat);
  TH1F *hPtHatRawPP = new TH1F("hPtHatRawPP","",nbinsPP_pthat,boundariesPP_pthat);
  
  RooUnfoldResponse res(uhist[0]->hResMeas,uhist[0]->hResTrue);
  
  cout<<"hi"<<endl;
  for (int i=0;i<nbins_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat,boundaries_pthat);
    data[i]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRawPPb->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
  
  for (int i=0;i<nbinsPP_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
    dataPP[i]->tJet->Project("hPtHatTmp","pthat","abs(vz)<15");
    hPtHatRawPP->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
  
  //apply the JECs to MC from Doga here. 
  TF1 *fgaus_ppb = new TF1("fgaus_ppb","gaus(0)",-20,20);
  fgaus_ppb->SetParameters(1,0,1);
  TF1 *f_ppb = new TF1("f_ppb","[0]/pow(x,[1])",50,300);
  f_ppb->SetParameters(1.052,0.5261);//for akPu3PF in MC for pPb/Pbp
  
  //apply the JECs to MC from Doga here. 
  TF1 *fgaus_pp = new TF1("fgaus_pp","gaus(0)",-20,20);
  fgaus_pp->SetParameters(1,0,1);
  TF1 *f_pp = new TF1("f_pp","[0]/pow(x,[1])",50,300);
  f_pp->SetParameters(0.5688,0.4867);//for ak3PF in MC for pp

  Float_t smeared_pt_ppb = 0,smeared_pt_pp = 0;

  // Fill PPb MC
  if (!useMatrixFromFile) {
    for (int i=0;i<nbins_pthat;i++) {
      if (xsection[i]==0) continue;
      cout <<"Loading pthat"<<boundaries_pthat[i]
	   <<" sample, cross section = "<<xsection[i]
	   << Form(" pthat>%.0f&&pthat<%.0f",boundaries_pthat[i],boundaries_pthat[i+1])<<endl;
      for (Long64_t jentry2=0; jentry2<data[i]->tJet->GetEntries();jentry2++) {
	data[i]->tEvt->GetEntry(jentry2);
	data[i]->tJet->GetEntry(jentry2);
	data[i]->tGenJet->GetEntry(jentry2);
	if(data[i]->pthat<boundaries_pthat[i] || data[i]->pthat>boundaries_pthat[i+1]) continue;
	int pthatBin = hPtHatPPb->FindBin(data[i]->pthat);
	float scale = (xsection[pthatBin-1]-xsection[pthatBin])/hPtHatRawPPb->GetBinContent(pthatBin);
	if(fabs(data[i]->vz)>15) continue;
	int cBin = hCent->FindBin(data[i]->bin)-1;
	double weight_cent=1;
	double weight_pt=1;
	double weight_vz=1;
				
	//weight_cent = fCentralityWeight->Eval(data[i]->bin);
	weight_vz = fVz->Eval(data[i]->vz);
	hCentMC->Fill(data[i]->bin,scale*weight_cent*weight_vz);
	hVzMC->Fill(data[i]->vz,scale*weight_cent*weight_vz);
	if (cBin>=nbins_cent) continue;
	if (cBin==-1) continue;
	hPtHatPPb->Fill(data[i]->pthat,scale*weight_cent*weight_vz);
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
	  if ( data[i]->refpt[k]  < 15. ) continue;
	  if ( data[i]->jteta[k]  > 1.465 || data[i]->jteta[k] < -0.535 ) continue; //eta CM assuming that pPb MC is forward beam
	  //if ( data[i]->chargedMax[k]/data[i]->jtpt[k]<0.01) continue;
	  //if (  ) continue;
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

	    //if(data[i]->jtpt[k]>=50 && data[i]->jtpt[k]<=300){
	    //smeared_pt_ppb = data[i]->jtpt[k]*(1+(f_ppb->Eval(data[i]->jtpt[k]))*fgaus_ppb->GetRandom());
	    //response[cBin]->Fill(smeared_pt_ppb,data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	    //uhist[cBin]-> hMatrix->Fill(data[i]->refpt[k],smeared_pt_ppb,scale*weight_cent*weight_pt*weight_vz);
	    //uhist[cBin]-> hGen->Fill(data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	    //uhist[cBin]-> hRecoMC->Fill(smeared_pt_ppb,scale*weight_cent*weight_pt*weight_vz);	      
	    //}else {
	      response[cBin]->Fill(data[i]->jtpt[k],data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	      uhist[cBin]-> hMatrix->Fill(data[i]->refpt[k],data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	      uhist[cBin]-> hGen->Fill(data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	      uhist[cBin]-> hRecoMC->Fill(data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	      //}
	  }
	  if (isMC&&jentry2>data[i]->tJet->GetEntries()/2. &&!isPyquen) {
	    uhist[cBin]-> hGen->Fill(data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
	    uhist[cBin]-> hMeas->Fill(data[i]->jtpt[k],scale*weight_cent*weight_pt*weight_vz);
	    uhist[cBin]-> hMeasJECSys->Fill(data[i]->jtpt[k]*(1.+0.02/nbins_cent*(nbins_cent-i)),scale*weight_cent*weight_pt*weight_vz);
	    
	    
	  }
	}
	
      }
      
    }
    
		
		
    ////// Pyquen cross check
		
    if(isPyquen && isMC){
			
      for (Long64_t jentry2=0; jentry2<data_pq->tJet->GetEntries();jentry2++) {
	data_pq->tEvt->GetEntry(jentry2);
	data_pq->tJet->GetEntry(jentry2);
	data_pq->tGenJet->GetEntry(jentry2);
	int cBin = hCent->FindBin(data_pq->bin)-1;
	double weight_cent=1;
	double weight_pt=1;
	double weight_vz=1;
				
	//weight_cent = fCentralityWeight->Eval(data_pq->bin);
	weight_vz = fVz->Eval(data_pq->vz);
	if (cBin>=nbins_cent) continue;
	if (cBin==-1) continue;
				
				
	for (int k= 0; k < data_pq->njets; k++) {
					
					
	  uhist[cBin]-> hGen->Fill(data_pq->refpt[k],weight_cent*weight_pt*weight_vz);
	  uhist[cBin]-> hMeas->Fill(data_pq->jtpt[k],weight_cent*weight_pt*weight_vz);
					
	}
      }
      
      
      
    }
    
    // fill pp MC
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
	//if(dataPP[i]->bin<=28) continue;//figure out why this cut is there? ask Yen-Jie 
	int pthatBin = hPtHatPP->FindBin(dataPP[i]->pthat);
	float scale = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/hPtHatRawPP->GetBinContent(pthatBin);
	double weight_cent=1;
	double weight_pt=1;
	double weight_vz=1;
	
	if(fabs(dataPP[i]->vz)>15) continue;

	weight_vz = fVz->Eval(dataPP[i]->vz);
	//if (weight_vz>5||weight_vz<0.5) cout <<dataPP[i]->vz<<" "<<weight_vz<<endl;
	//weight_vz = 1;
	hPtHatPP->Fill(dataPP[i]->pthat,scale*weight_vz);
	int hasLeadingJet = 0;
	hVzPPMC->Fill(dataPP[i]->vz,scale*weight_vz);

	//Float_t vz_temp = dataPP[i]->vz;
	//weight_vz = fVz->Eval(vz_temp);
	

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
	  if ( dataPP[i]->refpt[k]  < 15. ) continue;
	  if ( dataPP[i]->jteta[k]  > 1. || dataPP[i]->jteta[k] < -1. ) continue;

	  //if ( dataPP[i]->chargedMax[k]/dataPP[i]->jtpt[k]<0.01) continue;
	  //if (uhist[nbins_cent]->hMeasMatch!=0) {
	  //   int ptBinNumber = uhist[nbins_cent]->hMeasMatch->FindBin(dataPP[i]->jtpt[k]);
	  //   int ratio = uhist[nbins_cent]->hMeasMatch->GetBinContent(ptBinNumber);
	  //if (ratio!=0) weight_pt = 1./ratio;
	  //}
					
	  if (!isMC||jentry2<dataPP[i]->tJet->GetEntries()/2.) {

	    //if(dataPP[i]->jtpt[k]>=50 && dataPP[i]->jtpt[k]<=300){
	    //smeared_pt_pp = dataPP[i]->jtpt[k]*(1+(f_pp->Eval(dataPP[i]->jtpt[k]))*fgaus_pp->GetRandom());
	    //response[nbins_cent]->Fill(smeared_pt_pp,dataPP[i]->refpt[k],scale*weight_vz);
	    //uhist[nbins_cent]-> hMatrix->Fill(dataPP[i]->refpt[k],smeared_pt_pp,scale*weight_vz);
	    //uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);
	    //uhist[nbins_cent]-> hRecoMC->Fill(smeared_pt_pp,scale*weight_vz);
	    //}else{
	      response[nbins_cent]->Fill(dataPP[i]->jtpt[k],dataPP[i]->refpt[k],scale*weight_vz);
	      uhist[nbins_cent]-> hMatrix->Fill(dataPP[i]->refpt[k],dataPP[i]->jtpt[k],scale*weight_vz);
	      uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);
	      uhist[nbins_cent]-> hRecoMC->Fill(dataPP[i]->jtpt[k],scale*weight_vz);
	      //}
	  }
	  if (isMC&&jentry2>dataPP[i]->tJet->GetEntries()/2.) {
						
	    uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);
	    uhist[nbins_cent]-> hMeas->Fill(dataPP[i]->jtpt[k],scale*weight_vz);
	  }
	}
      }
    }
  }
 
  TCanvas *cMC = new TCanvas("cMC","MC",1000,800);
  cMC->Divide(nbins_cent,nbins_cent);
  TCanvas *cData = new TCanvas("cData","Data",1000,800);
  cData->Divide(nbins_cent,nbins_cent);
  for(int i = 0;i<=nbins_cent;i++){
		
    TLegend *title1 = 0,*title2 = 0;
		
    cMC->cd(i+1);
    cMC->cd(i+1)->SetLogy();
    title1 = myLegend(0.18,0.35,0.48,0.45);//MC
    title1->SetTextSize(0.06);
		
    uhist[i]->hGen->Draw();
    //if(i == 6){
    //  title1->AddEntry(uhist[i]->hGen,"PP MC","");
    //}else{
    //  title1->AddEntry(uhist[i]->hGen,Form("PbPb MC - %2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    //}
    title1->Draw();
		
    cData->cd(i+1);
    cData->cd(i+1)->SetLogy();
    uhist[i]->hMeas->Draw();
		
    title2 = myLegend(0.18,0.7,0.48,0.8);//data
    title2->SetTextSize(0.06);
    //if(i == 6){
    //  title2->AddEntry(uhist[i]->hMeas,"PP Data","");
    //}else {
    //  title2->AddEntry(uhist[i]->hMeas,Form("PbPb Data - %2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    //}
    title2->Draw();
  }
  
  cout<<"checking the input histograms"<<endl;
  for(int i = 0;i<=nbins_cent;i++){
    uhist[i]->hMeas->Print("base");
    uhist[i]->hGen->Print("base");
    cout<<endl<<endl;
  }

  /*
  TCanvas *cCent = new TCanvas("cCent","Centrality",600,600);
  divideBinWidth(hCentData);
  hCentMC->Scale(1./hCentMC->Integral(0,1000));
  hCentData->Scale(1./hCentData->Integral(0,1000));
  hCentMC->SetMarkerColor(2);
  hCentMC->SetLineColor(2);
  hCentData->Draw();
  hCentMC->Draw("same");
  */
  cout <<"Response Matrix..."<<endl;
	
  TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",1200,800);
  cMatrix->Divide(1,2);
  TCanvas* cResponseNorm = new TCanvas("cResponseNorm","Normalized Response Matrix",1200,800);
  cResponseNorm->Divide(1,2);
  for (int i=0;i<=nbins_cent;i++){
    cMatrix->cd(i+1);
    if (!useMatrixFromFile) {
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
    }
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
	
	
  TCanvas *cPtHat = new TCanvas("cPtHat","Pt Hat",600,600);
  cPtHat->SetLogy();
  hPtHatPPb->Draw();
	
  cout << "==================================== TEST =====================================" << endl;
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  char chmet[100];
	
  // ======================= Reconstructed pp and PPb spectra =========================================================
  TCanvas * cPPb = new TCanvas("cPPb","PPb",1200,800);
  cPPb->Divide(3,3);
  cPPb->cd(1);
	
	
  for (int i=0;i<nbins_cent;i++) {
    cPPb->cd(i+1)->SetLogy();
		
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
		
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
		
    bayesianUnfold myUnfolding(uhist[i]->hMatrix,myPrior.hPrior,0);
    myUnfolding.unfold(uhist[i]->hMeas,nBayesianIter); 
    
    //this is where you send the required histogram to be unfolded. 
    //im changing it from uhist[i]->hMeas to uhist[i]->hRecoMC

    bayesianUnfold mcUnfold(uhist[1]->hMatrix,myPrior.hPrior,0);
    mcUnfold.unfold(uhist[1]->hRecoMC,nBayesianIter);
    
    uhist[1]->hUnfoMC = (TH1F*)mcUnfold.hPrior->Clone(Form("Unfolded_Reco_MC_cent%d",1));
    
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
    uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJECSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    //uhist[i]->hRecoBinByBin = (TH1F*) unfold2.Hreco();
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));
		
    /*
    // Do unfolding
    prior myPrior(uhist[i]->hMatrix,uhist[i]->hMeas,0.0);
    myPrior.unfold(uhist[i]->hMeas,1);
    TH1F *hPrior = (TH1F*)uhist[i]->hMatrix->ProjectionX()->Clone(Form("hPrior_cent%d",i));
    hPrior->Scale(uhist[i]->hMeas->Integral(0,1000)/hPrior->Integral(0,1000));
    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;
    // Iteration Systematics
    for (int j=2;j<7;j++)
    {
    bayesianUnfold myUnfoldingSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingSys.unfold(uhist[i]->hMeas,j);
    uhist[i]->hRecoIterSys[j]  = (TH1F*) myUnfoldingSys.hPrior->Clone(Form("hRecoRAA_IterSys%d_cent%d",j,i));
    }
		 
    // Do Bin-by-bin
    TH1F *hBinByBinCorRaw = (TH1F*)uhist[i]->hResponse->ProjectionY();
    TH1F *hMCGen           = (TH1F*)uhist[i]->hResponse->ProjectionX(); // gen
    hBinByBinCorRaw->Divide(hMCGen);
    TF1 *f = new TF1("f","[0]+[1]*x");
    hBinByBinCorRaw->Fit("f","LL ","",100,300);
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone();//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //      TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    delete hBinByBinCorRaw,hMCGen;
    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
    //      uhist[i]->hRecoBinByBin = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
		 
    uhist[i]->hReco         = (TH1F*) uhist[i]->hRecoIterSys[nBayesianIter]->Clone(Form("Unfolded_cent%i",i));
    uhist[i]->hRecoJECSys   = (TH1F*) myUnfoldingJECSys.hPrior->Clone(Form("UnfoldedJeCSys_cent%i",i));
    uhist[i]->hRecoSmearSys   = (TH1F*) myUnfoldingSmearSys.hPrior->Clone(Form("UnfoldedSmearSys_cent%i",i));
    //uhist[i]->hRecoBinByBin = (TH1F*) unfold2.Hreco();
    uhist[i]->hRecoBinByBin->SetName(Form("UnfoldedBinByBin_cent%i",i));
		 
    */
    /*
    if (doToy) {
      TCanvas *cToy = new TCanvas("cToy","toy",600,600);
      int nExp=1000;
      TH1F *hTmp[nbins_truth+1];
      TH1F *hTmp2[nbins_truth+1];
      for (int j=1;j<=nbins_truth;j++) {
	hTmp[j] = new TH1F(Form("hTmp%d",j),"",200,0,10.+uhist[i]->hReco->GetBinContent(j)*2);
	hTmp2[j] = new TH1F(Form("hTmp2%d",j),"",200,0,10.+uhist[i]->hRecoBinByBin->GetBinContent(j)*2);
      }
      for (int exp =0; exp<nExp; exp++) {
	TH1F *hToy = (TH1F*)uhist[i]->hMeas->Clone();
	TH2F *hMatrixToy = (TH2F*)uhist[i]->hMatrix->Clone();
	hToy->SetName("hToy");
	if (exp%100==0) cout <<"Pseudo-experiment "<<exp<<endl;
	for (int j=1;j<=hToy->GetNbinsX();j++) {
	  double value = gRandom->Poisson(uhist[i]->hMeas->GetBinContent(j));
	  hToy->SetBinContent(j,value);
	}
				
	for (int j=1;j<=hMatrixToy->GetNbinsX();j++) {
	  for (int k=1;k<=hMatrixToy->GetNbinsY();k++) {
	    double value = gRandom->Gaus(uhist[i]->hMatrix->GetBinContent(j,k),uhist[i]->hMatrix->GetBinError(j,k));
	    hMatrixToy->SetBinContent(j,k,value);
	  }
	}
	//RooUnfoldBayes unfoldToy(response[i],hToy,2);
	prior myPriorToy(hMatrixToy,hToy,0.0);
	myPriorToy.unfold(hToy,1);
	bayesianUnfold myUnfoldingToy(hMatrixToy,hPrior,0.0);
	myUnfoldingToy.unfold(hToy,nBayesianIter);
	RooUnfoldBinByBin unfoldToy2(response[i],hToy);
	TH1F *hRecoTmp = (TH1F*) myUnfoldingToy.hPrior->Clone();
	TH1F *hRecoTmp2 = (TH1F*) unfoldToy2.Hreco();
				
	for (int j=1;j<=hRecoTmp->GetNbinsX();j++) {
	  hTmp[j]->Fill(hRecoTmp->GetBinContent(j));
	  hTmp2[j]->Fill(hRecoTmp2->GetBinContent(j));
	}
	delete hToy;
	delete hRecoTmp;
	delete hRecoTmp2;
	delete hMatrixToy;
      }
      TF1 *f = new TF1("f","[0]*TMath::Gaus(x,[1],[2])");
      for (int j=1;j<=nbins_truth;j++)
	{
	  f->SetParameters(hTmp[j]->GetMaximum(),hTmp[j]->GetMean(),hTmp[j]->GetRMS());
				
	  if (hTmp[j]->GetMean()>0) {
	    hTmp[j]->Fit("f","LL Q ");
	    hTmp[j]->Fit("f","LL Q ");
	    //	       cToy->SaveAs(Form("toy/cent-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
	    //     	       cout <<j<<" "<<f->GetParameter(2)<<endl;
	    uhist[i]->hReco->SetBinError(j,f->GetParameter(2));
	  }
	  f->SetParameters(hTmp2[j]->GetMaximum(),hTmp2[j]->GetMean(),hTmp2[j]->GetRMS());
	  if (hTmp2[j]->GetMean()>0) {
	    hTmp2[j]->Fit("f","LL Q ");
	    hTmp2[j]->Fit("f","LL Q ");
	    //cToy->SaveAs(Form("toy/cent2-%d-pt-%.0f.gif",i,uhist[i]->hReco->GetBinCenter(j)));
	    //cout <<j<<" "<<f->GetParameter(2)<<endl;
	    uhist[i]->hRecoBinByBin->SetBinError(j,f->GetParameter(2));
	  }
	  delete hTmp[j];
	  delete hTmp2[j];
	}
      cPbPb->cd(i+1);
    }
    */
    //cleanup(uhist[i]->hReco);
    uhist[i]->hMeas->SetMarkerStyle(20);
    uhist[i]->hMeas->SetMarkerColor(kRed);
    uhist[i]->hReco->SetMarkerStyle(25);
    uhist[i]->hReco->Draw("");
    uhist[i]->hReco->SetAxisRange(100,330);
    //      TH1F *hReproduced = (TH1F*)myUnfolding.hReproduced->Clone(Form("hReproduced_cent%d",i));
    //      hReproduced->SetMarkerColor(4);
    //      hReproduced->SetMarkerStyle(24);
    uhist[i]->hMeas->Draw("same");
  }
	
  /*
  cPbPb->cd(8);
  TLegend *pbpblegend = myLegend(0.52,0.65,0.85,0.9);
  uhist[0]->hMeas->SetMarkerStyle(20);
  uhist[0]->hMeas->SetMarkerColor(kRed);
  uhist[0]->hReco->SetMarkerStyle(25);
  pbpblegend->AddEntry(uhist[0]->hMeas,"Measured","pl");
  pbpblegend->AddEntry(uhist[0]->hReco,"Unfolded","pl");
  pbpblegend->Draw();
  */
	
	
	
  
	
	
  //***************         Calculation of RAA        ************
	
  SysData systematics;
	
  /*
	
  // PP histograms
	
  TH1F *hRebinPP         = rebin(uhist[nbins_cent]->hReco, "hRebinPP");
  TH1F *hRebinPP_Npart   = rebin_Npart(uhist[nbins_cent]->hReco, "hRebinPP_Npart");
  TH1F *hRebinBinByBinPP = rebin(uhist[nbins_cent]->hRecoBinByBin, "hRebinBinByBinPP");
  TH1F *hRecoPP          = (TH1F*)uhist[nbins_cent]->hReco->Clone("hRecoPP");
  TH1F *hMeasPP          = (TH1F*)uhist[nbins_cent]->hMeas->Clone("hMeasPP");
  TH1F *hRebinMeasPP     = rebin(uhist[nbins_cent]->hMeas, "hRebinMeasPP");
  TH1F *hRebinGenPP      = rebin(uhist[nbins_cent]->hGen, Form("hRebinGen_cent%d",nbins_cent));
	
	
  //dividebinwidth pp histograms
  divideBinWidth(hRebinPP);
  divideBinWidth(hRebinPP_Npart);
  divideBinWidth(hRebinBinByBinPP);
  divideBinWidth(hRecoPP);
  divideBinWidth(hMeasPP);
  divideBinWidth(hRebinMeasPP);
  divideBinWidth(hRebinGenPP);
	
	
  // Scale PP histograms
	
	
  hRebinPP               ->Scale(1./pplumi/64/1000000);
  hRebinPP_Npart		   ->Scale(1./pplumi/64/1000000);
  hRebinMeasPP           ->Scale(1./pplumi/64/1000000);
  hRecoPP                ->Scale(1./pplumi/64/1000000);
  hRebinBinByBinPP       ->Scale(1./pplumi/64/1000000);
  hMeasPP                ->Scale(1./pplumi/64/1000000);
  hRebinGenPP            ->Scale(1./pplumi/64/1000000);
  */
	
  //***************Scale Factor for trackcut correction ************
	
  float CorFac[6] = {1.0331,1.0331,1.0300,1.0259,1.0217,1.0114};
  //float CorFac[6] = {1.,1.,1.,1.,1.,1.};
	
  // correction factor for PbPb = 0.998, 0.998, 0.995, 0.991,0.987,0.977, pp = 0.966
	
	
	
	
	
	
	
	
  //************************     Making Canvas       ***********************
	
	
	
  TCanvas * cIterSys = new TCanvas("cIterSys","Iteration Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cIterSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	
  TCanvas * cJECSys = new TCanvas("cJECSys","JEC Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cJECSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
  TCanvas * cSmearSys = new TCanvas("cSmearSys","Smear Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSmearSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	
  TCanvas * cRpA = new TCanvas("cRpA","RpA",1200,800);
  makeMultiPanelCanvasWithGap(cRpA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
	
  TCanvas * cSys = new TCanvas("cSys","Total Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
  TCanvas * cIterSysPP = new TCanvas("cIterSysPP","Iteration Systematics for PP",800,600);
	
  TCanvas * cPPMCclosure = new TCanvas("cPPMCclosure","cPPMCclosure",600,450);
	
	
	
  //************   MC Closure Test For PP        *************
	
  /*	
  TLine *line = new TLine(100,1,300,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
	
  if(isMC)
    {
		
      cout <<"  MC Closure Test For PP   " << endl;
		
      cPPMCclosure->cd();
		
      hRebinMeasPP->Divide(hRebinGenPP);
      hRebinBinByBinPP->Divide(hRebinGenPP);
      hRebinPP->Divide(hRebinGenPP);
		
		
      hRebinPP->SetAxisRange(100,300,"X");
      hRebinPP->SetAxisRange(0,2,"Y");
      hRebinPP ->SetLineColor(kBlack);
      hRebinPP ->SetMarkerColor(kBlack);
      hRebinMeasPP->SetMarkerStyle(24);
      hRebinMeasPP->SetLineColor(4);
      hRebinMeasPP->SetMarkerColor(4);
      hRebinBinByBinPP->SetMarkerStyle(33);
      hRebinBinByBinPP->SetLineColor(kRed);
      hRebinBinByBinPP->SetMarkerColor(kRed);
		
      makeHistTitle(hRebinPP,"","Jet p_{T} (GeV/c)","Reco / Truth");
      hRebinPP->GetYaxis()->SetTitleOffset(1.5);
      hRebinPP->GetXaxis()->SetTitleOffset(1.3);
      hRebinPP->Draw("");
      hRebinBinByBinPP->Draw("same");
      hRebinMeasPP->Draw("same");
      line->Draw();
		
      TLegend *legpp = myLegend(0.52,0.65,0.85,0.9);
      legpp->AddEntry(hRebinPP,"pp Bayesian","pl");
      legpp->AddEntry(hRebinBinByBinPP,"pp Bin-by-bin","pl");
      legpp->AddEntry(hRebinMeasPP,"pp no unfolding","pl");
      legpp->Draw();
      putCMSPrel(0.2,0.83,0.06);
      drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,20);
      drawText("PYTHIA",0.6,0.4,22);
      drawText("| #eta | <2 ",0.6,0.31,22);
		
		
    }
	
  */	
  //************    Iteration Sys PP  ak3PF        *************
  /*	
  cout <<" Plotting  Iteration Sys  for PP  akPu3PF   " << endl;
  cIterSysPP->cd();
	
  // Iteration systematics
  TH1F *hRecoIterSysPP[10];
  TH1F *hRebinPP_tmp         = rebin(uhist[nbins_cent]->hReco, "hRebinPP_tmp");
  TLegend *legBayesianIterPP = myLegend(0.6,0.7,0.9,0.9);
	
  for (int j=2;j<7;j++) {
    hRecoIterSysPP[j] = rebin(uhist[nbins_cent]->hRecoIterSys[j],Form("hRecoIterSysPP_IterSys%d",j));
    hRecoIterSysPP[j]->SetLineColor(colorCode[j-2]);
    hRecoIterSysPP[j]->SetMarkerColor(colorCode[j-2]);
    hRecoIterSysPP[j]->Divide(hRebinPP_tmp);
    if (j==2){
      makeHistTitle(hRecoIterSysPP[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
      hRecoIterSysPP[j]->SetAxisRange(0,2,"Y");
      hRecoIterSysPP[j]->Draw();
    } else {
      hRecoIterSysPP[j]->Draw("same");
    }
		
    checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoIterSysPP[j],0,1.1);
    legBayesianIterPP->AddEntry(hRecoIterSysPP[j],Form("Iteration %d",j),"pl");
  }
  legBayesianIterPP->Draw();
  line->Draw();
  drawEnvelope(systematics.hSysIter[nbins_cent],"hist same");
	
	
  */
	
	
	
  // calculate the maximum from all cent bins
  for (int i=0;i<nbins_cent;i++) {
		
    TH1F *hRecoRpAIterSys[10];
    TH1F *hRebinRpA = rebin(uhist[i]->hReco, Form("hRebinRpA_cent%d",i));
		
    for (int j=2;j<7;j++) {
      hRecoRpAIterSys[j] = rebin(uhist[i]->hRecoIterSys[j],Form("hRecoRpA_IterSys%d_cent%d",j,i));
      hRecoRpAIterSys[j]->SetLineColor(colorCode[j-2]);
      hRecoRpAIterSys[j]->SetMarkerColor(colorCode[j-2]);
      hRecoRpAIterSys[j]->Divide(hRebinRpA);
      checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoRpAIterSys[j],0,1.1);
    }
		
		
  }
  for (int i=0;i<nbins_cent;i++) {
    //cRpA->cd(nbins_cent-i);
    
    TH1F *hRebinRpA = rebin(uhist[i]->hReco, Form("hRebinRpA_cent%d",i));

    /*		
    TH1F *hRebinRAA = rebin(uhist[i]->hReco, Form("hRebinRAA_cent%d",i));
    TH1F *hRebinRAA_Npart = rebin_Npart(uhist[i]->hReco, Form("hRebinRAA_Npart_cent%d",i));
    TH1F *hRebinBinByBinRAA = rebin(uhist[i]->hRecoBinByBin, Form("hRebinBinByBinRAA_cent%d",i));
    TH1F *hRebinMeasRAA = rebin(uhist[i]->hMeas, Form("hRebinMeasRAA_cent%d",i));
    TH1F *hRecoRAA         = (TH1F*)uhist[i]->hReco->Clone(Form("hRecoRAA_cent%d",i));
    TH1F *hMeasRAA         = (TH1F*)uhist[i]->hMeas->Clone(Form("hMeasRAA_cent%d",i));
    TH1F *hRecoRAAJECSys   = (TH1F*)uhist[i]->hRecoJECSys->Clone(Form("hRecoRAAJECSys_cent%d",i));
    TH1F *hRecoRAASmearSys   = (TH1F*)uhist[i]->hRecoSmearSys->Clone(Form("hRecoRAASmearSys_cent%d",i));
    TH1F *hRebinRAAJECSys = rebin(hRecoRAAJECSys, Form("hRebinRAAJECSys_cent%d",i));
    TH1F *hRebinRAASmearSys = rebin(hRecoRAASmearSys, Form("hRebinRAASmearSys_cent%d",i));
		
    // dividebin width
    divideBinWidth(hRebinRAA);
    divideBinWidth(hRebinRAA_Npart);
    divideBinWidth(hRebinBinByBinRAA);
    divideBinWidth(hRebinMeasRAA);
    divideBinWidth(hRebinRAAJECSys);
    divideBinWidth(hRebinRAASmearSys);
    divideBinWidth(hRecoRAA);
    divideBinWidth(hMeasRAA);
    divideBinWidth(hRecoRAAJECSys);
    divideBinWidth(hRecoRAASmearSys);
    
    ///// Set Histograms of Genjet to Save
    uhist[i]-> hGen->SetName(Form("hGen_cent%i",i));
    
    TLine *l = new TLine(100,1,300,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    TLegend *title=0;
    if(!isMC) title = myLegend(0.18,0.7,0.48,0.8);//data
    if (isMC) title = myLegend(0.18,0.35,0.48,0.45);//MC
    title->AddEntry(hRecoRAA,Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    title->SetTextSize(0.06);
    */
    
    // Iteration systematics
    TH1F *hRecoRpAIterSys[10];
    TLegend *legBayesianIter = myLegend(0.6,0.6,0.9,0.9);
    legBayesianIter->SetTextSize(0.042);
    
    for (int j=2;j<7;j++) {
      
      
      //************    Iteration Sys   akPu3PF        *************
      
      cout <<" Plotting  Iteration Sys   akPu3PF   " << endl;
      cIterSys->cd(nbins_cent-i);
      
      hRecoRpAIterSys[j] = rebin(uhist[i]->hRecoIterSys[j],Form("hRecoRpA_IterSys%d_cent%d",j,i));
      hRecoRpAIterSys[j]->SetLineColor(colorCode[j-2]);
      hRecoRpAIterSys[j]->SetMarkerColor(colorCode[j-2]);
      hRecoRpAIterSys[j]->Divide(hRebinRpA);
      if (j==2){
	makeHistTitle(hRecoRpAIterSys[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
	hRecoRpAIterSys[j]->SetAxisRange(0,2,"Y");
	hRecoRpAIterSys[j]->Draw();
      } else {
	hRecoRpAIterSys[j]->Draw("same");
      }
      
      checkMaximumSys(systematics.hSysIter[i],hRecoRpAIterSys[j],0,1.1);
      legBayesianIter->AddEntry(hRecoRpAIterSys[j],Form("Iteration %d",j),"pl");
    }
    if (i==nbins_cent-1) legBayesianIter->Draw();
    //title->Draw();
    //l->Draw();
    
    if (i==nbins_cent-1) {
      putCMSPrel(0.2,0.83,0.06);
      drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,21);
    }
    
    if (i==nbins_cent-2)
      {   drawText("PbPb         #sqrt{s_{NN}} = 2.76 TeV",0.2,0.83,19);
	drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
      }
    
    if (i==nbins_cent-4)
      {	drawText("Bayesian",0.68,0.83,19);
	drawText("| #eta | < 2 ",0.28,0.83,24);
      }
		
    DrawPanelLabel(i);
    
    
    if (useFixedIterativeSys) {
      // always use central systematics
      systematics.hSysIter[i] = (TH1F*) systematics.hSysIter[nbins_cent]->Clone();
    }
    drawEnvelope(systematics.hSysIter[i],"hist same");
    
    
    //calculate the RpA here: 
    // histograms needed: hRebinRAA and 
    
    // do it in another macro - for sake of speed. 


    /*	
		
		
		
    cRAA->cd(nbins_cent-i);
    makeHistTitle(hRebinRAA,"","Jet p_{T} (GeV/c)","Jet R_{AA}");
    if (!isMC) {
			
      // Scale PbPb Hisotograms  and Calculate RAA
			
      //hRebinRAA            ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRebinRAA_Npart      ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRebinMeasRAA        ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRebinBinByBinRAA    ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRecoRAA             ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hMeasRAA             ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRecoRAAJECSys       ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRebinRAASmearSys     ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRebinRAAJECSys       ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      //hRecoRAASmearSys     ->Scale(1./CorFac[i]/lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
			
      hRebinRAA            ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRebinRAA_Npart      ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRebinMeasRAA        ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRebinBinByBinRAA    ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRecoRAA             ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hMeasRAA             ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRecoRAAJECSys       ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRebinRAASmearSys    ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRebinRAAJECSys      ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
      hRecoRAASmearSys     ->Scale(1./lumi/7.65/1000000/0.025/(boundaries_cent[i+1]-boundaries_cent[i])/ncoll[i]);
			
      hRebinMeasRAA->Divide(hRebinMeasPP);
      hRebinRAA->Divide(hRebinPP);
      hRebinRAA_Npart->Divide(hRebinPP_Npart);
      hRebinRAAJECSys->Divide(hRebinPP);
      hRebinRAASmearSys->Divide(hRebinPP);
      hRecoRAA->Divide(hRecoPP);
      hRebinBinByBinRAA->Divide(hRebinBinByBinPP);
      hMeasRAA->Divide(hMeasPP);
      hRecoRAAJECSys->Divide(hRecoPP);
      hRecoRAASmearSys->Divide(hRecoPP);
      hRebinRAA->SetYTitle("Jet R_{AA}");
			
      TH1F *hRebinRAASmearSys = rebin(hRecoRAASmearSys, Form("hRebinRAASmearSys_cent%d",i));
    } else {
			
      //************   MC Closure Test For PbPb        *************
			
			
      cout <<"  MC Closure Test For PbPb   " << endl;
			
			
      TH1F *hRebinGen        = rebin(uhist[i]->hGen, Form("hRebinGen_cent%d",i));
      hRebinRAA->SetYTitle("Reco / Truth");
      hRebinMeasRAA->Divide(hRebinGen);
      hRebinRAA->Divide(hRebinGen);
      hRebinRAAJECSys->Divide(hRebinGen);
      hRebinRAASmearSys->Divide(hRebinGen);
      hRecoRAA->Divide(uhist[i]->hGen);
      hRebinBinByBinRAA->Divide(hRebinGen);
      hMeasRAA->Divide(uhist[i]->hGen);
      hRecoRAAJECSys->Divide(uhist[i]->hGen);
			
    }
    hRebinRAA->SetAxisRange(100,330,"X");
    hRebinRAA->SetAxisRange(0,2,"Y");
    hRebinMeasRAA ->SetMarkerStyle(24);
    hMeasRAA ->SetLineColor(5);
    hMeasRAA ->SetMarkerColor(5);
    hRebinMeasRAA ->SetLineColor(kBlack);
    hRebinMeasRAA ->SetMarkerColor(kBlack);
    hRecoRAA->SetMarkerStyle(24);
    hRecoRAA->SetLineColor(4);
    hRecoRAA->SetMarkerColor(4);
    hRebinBinByBinRAA->SetMarkerStyle(33);
    hRebinBinByBinRAA->SetLineColor(kRed);
    hRebinBinByBinRAA->SetMarkerColor(kRed);
		
    //************   JEC Sys  akPu3PF        *************
		
		
    hRebinRAAJECSys->Divide(hRebinRAAJECSys,hRebinRAA,1,1,"B");
    hRebinRAAJECSys->SetAxisRange(0.61,1.39,"Y");
    hRebinRAAJECSys->Draw("p");
    makeHistTitle(hRebinRAAJECSys,"","Jet p_{T} (GeV/c)","Ratio",2);
    hRebinRAAJECSys->SetAxisRange(100,330,"X");
    TF1 *fPol = new TF1("fPol","[0]+[1]*x");
    hRebinRAAJECSys->Fit("fPol","");
		
		
		
		
    cout <<" Plotting JEC Sys  akPu3PF  " << endl;
		
    cJECSys->cd(nbins_cent-i);
    hRebinRAAJECSys->Draw("p");
    checkMaximumSys(systematics.hSysJEC[i],functionHist(fPol,systematics.hSysJEC[i],Form("hist_sysJEC_cent%d",i)));
    l->Draw();
		
		
    //************    SmearSys   akPu3PF        *************
		
    cout <<" Plotting SmearSys akPu3PF " << endl;
		
    cSmearSys->cd(nbins_cent-i);
		
    hRebinRAASmearSys->Divide(hRebinRAASmearSys,hRebinRAA,1,1,"B");
    hRebinRAASmearSys->SetAxisRange(0.61,1.39,"Y");
    hRebinRAASmearSys->Draw("p");
    makeHistTitle(hRebinRAASmearSys,"","Jet p_{T} (GeV/c)","Ratio",2);
    hRebinRAASmearSys->SetAxisRange(100,330,"X");
    hRebinRAASmearSys->Fit("fPol","");
    hRebinRAASmearSys->Draw("p");
    checkMaximumSys(systematics.hSysSmear[i],functionHist(fPol,systematics.hSysSmear[i],Form("hist_sysSmear_cent%d",i)));
    l->Draw();
		
		
    title->Draw();
		
    */	
		
		
    //************    Different methods RAA in akPu3PF        *************
    /*	
		
    cout <<" Plotting different methods RAA in akPu3PF " << endl;
		
    cRAA->cd(nbins_cent-i);
    hRebinRAA->Draw();
    //if(!isMC){
    //  systematics.calcTotalSys(i);
    //  systematics.Draw(hRebinRAA,i);
    //}
		
    DrawPanelLabel(i);
		
    ////// Put Correlated and Uncorrelated Uncertainties in the Bayesian Unfolding histogram
    TH1F *hRebinRAA_corr = new TH1F(*hRebinRAA);
    for (Int_t j = 1; j < hRebinRAA_corr->GetNbinsX() + 1; j++) {
      Float_t x = hRebinRAA_corr->GetXaxis()->GetBinCenter(j);
      Float_t s = 0;
			
      switch (algo) {
      case 2:
	break;
      case 3:
	switch (i) {
	case 0:
	  s = 2.044234807705282 + 0.0028322876264653307*x + 3.2827115118940787e-6*pow(x,2);
	  break;
	case 1:
	  s = 0.9941448362655844 + 0.011261508479420763*x - 0.000010142050572618295*pow(x,2);
	  break;
	case 2:
	  s = 0.6746697572381065 + 0.02023907788830483*x - 0.00004519177768958705*pow(x,2);
	  break;
	case 3:
	  s = 0.622822048623053 + 0.01784370943830341*x - 0.00003605323648754686*pow(x,2);
	  break;
	case 4:
	  s = 0.7997724959427492 + 0.013842026536881468*x - 0.000020916084265574625*pow(x,2);
	  break;
	case 5:
	  s = 0.5639310618959498 + 0.015922843743525164*x - 0.000030872114735102704*pow(x,2);
	  break;
	}
	break;
      case 4:
	break;
      }
      hRebinRAA_corr->SetBinError(j, hRebinRAA_corr->GetBinError(j) * s);
    }
    hRebinRAA_corr->SetLineColor(kPink-4);
    hRebinRAA_corr->SetLineWidth(5);
    //hRebinRAA_corr->Draw("e0x0same");
    //hRebinRAA_corr->Draw("same");
		
		
		
		
    title->Draw();
		
    TGraphErrors *smearing, *smearing_sys;
    TGraphErrors *gsvd;
    if (!isMC) {
      smearing = new TGraphErrors(Form("data/smearing/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      smearing->SetName("smearing");
      smearing->SetMarkerStyle(21);
      smearing->SetMarkerColor(kBlue);
      smearing->SetLineColor(kBlue);
      smearing->Draw("p same");
			
			
			
      // ****************Pawan's Editting******************************
			
      smearing_sys = new TGraphErrors(Form("data/smearing/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      smearing_sys->SetName("smearing_sys");
      smearing_sys->SetMarkerStyle(21);
      smearing_sys->SetMarkerColor(kBlue);
      smearing_sys->SetLineColor(kBlue);
			
			
      double fx=0,fy=0;
      for(int ibin=0;ibin<smearing->GetN();ibin++){
	smearing->GetPoint(ibin, fx, fy);
	double erSy = smearing->GetErrorX(ibin);  //! systematics
	double erSt = smearing->GetErrorY(ibin);  //! statistics
	smearing->SetPointError(ibin,0,erSt);     //! only statistical
	smearing_sys->SetPointError(ibin,0,erSy); //! only systematic
				
	TBox *b = new TBox(hRebinRAA_corr->GetBinLowEdge(ibin+1),fy-erSy,hRebinRAA_corr->GetBinLowEdge(ibin+2),fy+erSy);
	b->SetFillColor(kViolet+6);
	b->SetFillStyle(3006);
	b->SetLineColor(kViolet+6);
	//b->Draw();
      }
			
			
			
      // ****************   Pawan's Editting End    ******************************
    }
		
    hRebinRAA->SetMarkerStyle(20);
    hRebinRAA->SetMarkerColor(kBlack);
    hRebinRAA->Draw("same");
    hRebinMeasRAA->Draw("same");
		
		
    if (!isMC) {
      gsvd = new TGraphErrors(Form("data/gsvd/%.0f-%.0f.dat",boundaries_cent[i]*2.5,boundaries_cent[i+1]*2.5));
      gsvd->SetMarkerStyle(34);
      gsvd->SetMarkerColor(kGreen+3);
      gsvd->SetLineColor(kGreen+3);
      gsvd->Draw("p same");
    }
		
    hRebinBinByBinRAA->Draw("same");
		
		
    title->Draw();
    l->Draw();
		
    if (i==nbins_cent-1) {
      TLegend *leg = myLegend(0.6,0.65,0.95,0.9);
      leg->SetTextSize(0.05);
      leg->AddEntry(hRebinRAA,"Bayesian","pl");
      leg->AddEntry(hRebinBinByBinRAA,"Bin-by-bin","pl");
      if (!isMC){
	leg->AddEntry(smearing,"Smearing","pl");
	leg->AddEntry(gsvd,"GSVD","pl");
      }
      leg->AddEntry(hRebinMeasRAA,"No unfolding","pl");
      leg->Draw();
      putCMSPrel(0.2,0.83,0.06);
      drawText("Anti-k_{T} Particle Flow Jets   R = 0.3",0.2,0.23,20);
    }
    //uhist[i]->hMeas->Write();
		
    if (i==nbins_cent-2 &&!isMC)
      {   drawText("PbPb        #sqrt{s_{NN}} = 2.76 TeV",0.2,0.83,19);
	drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
      }
    if (i==nbins_cent-2 && isMC)
			
      {   drawText("PYTHIA+HYDJET",0.5,0.83,22);
	drawText("| #eta | <2 ",0.5,0.75,22);
			
      }
    if (i==nbins_cent-4 &&!isMC)
      drawText("| #eta | < 2 ",0.28,0.83,24);
		
    */	
		
    /*		
    //************     Total Sys       *************
		
		
    TLegend *title_ = myLegend(0.18,0.8,0.48,0.9);
    title_->AddEntry(hRecoRAA,Form("%2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    title_->SetTextSize(0.06);
		
		
    cout <<" Plotting Total Systematics " << endl;
		
    cSys->cd(nbins_cent-i);
    systematics.DrawComponent(i);
		
    DrawPanelLabel(i);
		
		
    if (i==nbins_cent-1) {
      putCMSPrel(0.5,0.2,0.06);
			
    }
		
    if (i==nbins_cent-2)
      {   drawText("PbPb  #sqrt{s_{NN}} = 2.76 TeV",0.5,0.83,19);
	drawText("#int L dt = 150 #mub^{-1}",0.5,0.72,19);
	drawText("Anti-k_{T} PF Jets   R = 0.3",0.4,0.21,19);
      }
		
    title_->Draw();
    l->Draw();
    */
  }   // centrality for loop end
	
	
	
	
  //************   Test  Plotting ends      ***************
	
	
	
  if (isMC) {
    cRpA->Update();
    cRpA->SaveAs("MCClosureTest/MCClosureTest.gif");
    cRpA->SaveAs("MCClosureTest/MCClosureTest.C");
    cRpA->SaveAs("MCClosureTest/MCClosureTest.pdf");
    cIterSys->Update();
    cIterSys->SaveAs("MCClosureTest/IterSys.gif");
    cIterSys->SaveAs("MCClosureTest/IterSys.C");
    cIterSys->SaveAs("MCClosureTest/IterSys.pdf");
    cPPMCclosure->Update();
    cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.gif");
    cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.C");
    cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.pdf");
		
		
  } else {
		
    /*
   
    cPPb->Update();
    cPPb->SaveAs(Form("result-2013-ppb-%s-cent-%d/data_vs_unfo.gif",algoName[algo],nbins_cent),"RECREATE");
    cPPb->SaveAs(Form("result-2013-ppb-%s-cent-%d/data_vs_unfo.C",algoName[algo],nbins_cent),"RECREATE");
    cPPb->SaveAs(Form("result-2013-ppb-%s-cent-%d/data_vs_unfo.pdf",algoName[algo],nbins_cent),"RECREATE");
    cIterSys->Update();
    cIterSys->SaveAs(Form("result-2013-ppb-%s-cent-%d/IterSys.gif",algoName[algo],nbins_cent),"RECREATE");
    cIterSys->SaveAs(Form("result-2013-ppb-%s-cent-%d/IterSys.C",algoName[algo],nbins_cent),"RECREATE");
    cIterSys->SaveAs(Form("result-2013-ppb-%s-cent-%d/IterSys.pdf",algoName[algo],nbins_cent),"RECREATE");
    */
    cMC->Update();
    cMC->SaveAs(Form("result-2013-ppb-%s-cent-%d/MC.gif",algoName[algo],nbins_cent),"RECREATE");
    cMC->SaveAs(Form("result-2013-ppb-%s-cent-%d/MC.C",algoName[algo],nbins_cent),"RECREATE");
    cMC->SaveAs(Form("result-2013-ppb-%s-cent-%d/MC.pdf",algoName[algo],nbins_cent),"RECREATE");
    cData->Update();
    cData->SaveAs(Form("result-2013-ppb-%s-cent-%d/Data.gif",algoName[algo],nbins_cent),"RECREATE");
    cData->SaveAs(Form("result-2013-ppb-%s-cent-%d/Data.C",algoName[algo],nbins_cent),"RECREATE");
    cData->SaveAs(Form("result-2013-ppb-%s-cent-%d/Data.pdf",algoName[algo],nbins_cent),"RECREATE");
    
    
    cMatrix->Update();
   
    cMatrix->SaveAs(Form("result-2013-ppb-%s-cent-%d/ResponseMatrix.jpg",algoName[algo],nbins_cent),"RECREATE");
    cMatrix->SaveAs(Form("result-2013-ppb-%s-cent-%d/ResponseMatrix.C",algoName[algo],nbins_cent),"RECREATE");
    cMatrix->SaveAs(Form("result-2013-ppb-%s-cent-%d/ResponseMatrix.pdf",algoName[algo],nbins_cent),"RECREATE");
    cResponseNorm->Update();
    cResponseNorm->SaveAs(Form("result-2013-ppb-%s-cent-%d/NormalizedResponseMatrix.jpg",algoName[algo],nbins_cent),"RECREATE");
    cResponseNorm->SaveAs(Form("result-2013-ppb-%s-cent-%d/NormalizedResponseMatrix.C",algoName[algo],nbins_cent),"RECREATE");
    cResponseNorm->SaveAs(Form("result-2013-ppb-%s-cent-%d/NormalizedResponseMatrix.pdf",algoName[algo],nbins_cent),"RECREATE");
    
		
  }
  divideBinWidth(hCent);
	

	
  


		
  ppb_Unfo->Write();
  //pbpb_Unfo->Close();
	
  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;
	
}






