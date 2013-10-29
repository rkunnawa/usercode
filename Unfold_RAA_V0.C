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
#include "headers/bayesianUnfold.h"
#include "headers/prior.h"

using namespace std;


//==============================================================================
// Unfolding Ying Lu 08 07 11
// Update Yen-Jie Lee 06.22.12
//==============================================================================

void Unfold_RAA_V0(int method = 1,int algo = 3,bool useSpectraFromFile = 0, bool useMatrixFromFile = 0, int doToy = 0, int isMC = 0,char *spectraFileName = "pbpbSpectra.root", int doJECSys = 0) // algo 2 =akpu2 ; 3 =akpu3 ; 4 =akpu4 ;1 = icpu5
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
	
  int nBayesianIter = 6;// 4 or 6 
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
	
  // ================ Bin Size ======================================================================
	
	
  // ================ PbPb PtBin ======================================================================
	
	
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
	
	
  ////// New MC samples
	
  if (yinglu) {
    boundaries_pthat[0]=30;
    fileName_pthat[0]="/hadoop/store/user/belt/hiForest2/Pythia30_HydjetDrum_mix01_HiForest2_v19.root";   
    xsection[0]= 1.079e-02;
	
    boundaries_pthat[1]=50;
    fileName_pthat[1]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet50_HydjetDrum_v27_mergedV1.root";   
    xsection[1]= 1.021e-03;
	
    boundaries_pthat[2]=80;
    fileName_pthat[2]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet80_HydjetDrum_v27_mergedV1.root";   
    xsection[2]= 9.913e-05;
	
    boundaries_pthat[3]=100;
    fileName_pthat[3]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet100_HydjetDrum_v27_mergedV1.root ";   
    xsection[3]= 3.069e-05 ;
	
    boundaries_pthat[4]=120;
    fileName_pthat[4]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet120_HydjetDrum_v27_mergedV1.root";   
    xsection[4]= 1.128e-05;
	
    boundaries_pthat[5]=170;
    fileName_pthat[5]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet170_HydjetDrum_v27_mergedV1.root";   
    xsection[5]= 1.470e-06;
	
    boundaries_pthat[6]=200;
    fileName_pthat[6]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet200_HydjetDrum_v28_mergedV1.root";   
    xsection[6]= 5.310e-07;
	
    boundaries_pthat[7]=250;
    fileName_pthat[7]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet250_HydjetDrum_v28_mergedV1.root";   
    xsection[7]= 1.192e-7;
	
    boundaries_pthat[8]=300;
    fileName_pthat[8]="/hadoop/store/user/belt/hiForest2/v27_v28/Dijet300_HydjetDrum_v28_mergedV1.root";   
    xsection[8]= 3.176e-08;
	
    xsection[9] = 0;
    boundaries_pthat[9]=1000;
  } else { 
	
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

    // yen-jie's file location
    /*
      boundaries_pthat[0]=30;
      fileName_pthat[0]="/d102/yjlee/hiForest2MC/Pythia30_HydjetDrum_mix01_HiForest2_v19.root";   
      xsection[0]= 1.079e-02;
	
      boundaries_pthat[1]=50;
      fileName_pthat[1]="/d102/yjlee/hiForest2MC/Pythia50_HydjetDrum_mix01_HiForest2_v19.root";   
      xsection[1]= 1.021e-03;
	
      boundaries_pthat[2]=80;
      fileName_pthat[2]="/d102/yjlee/hiForest2MC/Pythia80_HydjetDrum_mix01_HiForest2_v20.root";   
      xsection[2]= 9.913e-05;
	
      boundaries_pthat[3]=100;
      fileName_pthat[3]="/d102/yjlee/hiForest2MC/v27/mergedFile100v8.root";   
      xsection[3]= 3.069e-05 ;

      boundaries_pthat[4]=120;
      fileName_pthat[4]="/d102/yjlee/hiForest2MC/Pythia120_HydjetDrum_mix01_HiForest2_v21_ivan.root";   
      xsection[4]= 1.128e-05;
	
      boundaries_pthat[5]=170;
      fileName_pthat[5]="/d102/yjlee/hiForest2MC/Pythia170_HydjetDrum_mix01_HiForest2_v19.root";   
      xsection[5]= 1.470e-06;
	
      boundaries_pthat[6]=200;
      fileName_pthat[6]="/d102/yjlee/hiForest2MC/Pythia200_HydjetDrum_mix01_HiForest2_v21_ivan.root";   
      xsection[6]= 5.310e-07;
	
      boundaries_pthat[7]=250;
      fileName_pthat[7]="/d102/yjlee/hiForest2MC/Pythia250_HydjetDrum_mix01_HiForest2_v21_ivan.root";   
      xsection[7]= 1.192e-7;
	
      boundaries_pthat[8]=300;
      fileName_pthat[8]="/d102/yjlee/hiForest2MC/Pythia300_HydjetDrum_mix01_HiForest2_v21_ivan.root";   
      xsection[8]= 3.176e-08;
	
      xsection[9] = 0;
      boundaries_pthat[9]=1000;	
    */
  }
      	
	
  // ================ pp PtBin ======================================================================
  const int nbinsPP_pthat = 8;
  Double_t boundariesPP_pthat[nbinsPP_pthat+1];
  char *fileNamePP_pthat[nbinsPP_pthat+1];
  Double_t xsectionPP[nbinsPP_pthat+1];
	
  if (yinglu) {	
    boundariesPP_pthat[0]=30;
    fileNamePP_pthat[0]="/hadoop/store/user/belt/hiForest2/pp276Dijet30_merged.root";   
    xsectionPP[0]= 1.079e-02;
	
    boundariesPP_pthat[1]=50;
    fileNamePP_pthat[1]="/hadoop/store/user/belt/hiForest2/pp276Dijet50_merged.root";   
    xsectionPP[1]= 1.021e-03;
	
    boundariesPP_pthat[2]=80;
    fileNamePP_pthat[2]="/hadoop/store/user/belt/hiForest2/pp276Dijet80_merged.root";   
    xsectionPP[2]= 9.913e-05;
	
    boundariesPP_pthat[3]=120;
    fileNamePP_pthat[3]="/hadoop/store/user/belt/hiForest2/pp276Dijet120_merged.root";   
    xsectionPP[3]= 1.128e-05;
	
    boundariesPP_pthat[4]=170;
    fileNamePP_pthat[4]="/hadoop/store/user/belt/hiForest2/pp276Dijet170_merged.root";   
    xsectionPP[4]= 1.470e-06;
	
    boundariesPP_pthat[5]=200;
    fileNamePP_pthat[5]="/hadoop/store/user/belt/hiForest2/pp276Dijet200_merged.root";   
    xsectionPP[5]= 5.310e-07;
	
    boundariesPP_pthat[6]=250;
    fileNamePP_pthat[6]="/hadoop/store/user/belt/hiForest2/pp276Dijet250_merged.root";   
    xsectionPP[6]= 1.192e-7;
	
    boundariesPP_pthat[7]=300;
    fileNamePP_pthat[7]="/hadoop/store/user/belt/hiForest2/pp276Dijet300_merged.root";   
    xsectionPP[7]= 3.176e-08;
	
    xsectionPP[8] = 0;
    boundariesPP_pthat[8]=1000;
  } else {

	
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
    // yen-jie's file location
    boundariesPP_pthat[0]=30;
    fileNamePP_pthat[0]="/d102/yjlee/hiForest2MC/pp276Dijet30_merged.root";   
    xsectionPP[0]= 1.079e-02;
	
    boundariesPP_pthat[1]=50;
    fileNamePP_pthat[1]="/d102/yjlee/hiForest2MC/pp276Dijet50_merged.root";   
    xsectionPP[1]= 1.021e-03;
	
    boundariesPP_pthat[2]=80;
    fileNamePP_pthat[2]="/d102/yjlee/hiForest2MC/pp276Dijet80_merged.root";   
    xsectionPP[2]= 9.913e-05;
	
    boundariesPP_pthat[3]=120;
    fileNamePP_pthat[3]="/d102/yjlee/hiForest2MC/pp276Dijet120_merged.root";   
    xsectionPP[3]= 1.128e-05;
	
    boundariesPP_pthat[4]=170;
    fileNamePP_pthat[4]="/d102/yjlee/hiForest2MC/pp276Dijet170_merged.root";   
    xsectionPP[4]= 1.470e-06;
	
    boundariesPP_pthat[5]=200;
    fileNamePP_pthat[5]="/d102/yjlee/hiForest2MC/pp276Dijet200_merged.root";   
    xsectionPP[5]= 5.310e-07;
	
    boundariesPP_pthat[6]=250;
    fileNamePP_pthat[6]="/d102/yjlee/hiForest2MC/pp276Dijet250_merged.root";   
    xsectionPP[6]= 1.192e-7;
	
    boundariesPP_pthat[7]=300;
    fileNamePP_pthat[7]="/d102/yjlee/hiForest2MC/pp276Dijet300_merged.root";   
    xsectionPP[7]= 3.176e-08;
	
    xsectionPP[8] = 0;
    boundariesPP_pthat[8]=1000;
    */
  }
	
  //*******************lumi number for the sample***************************
  float lumi=150.;
  float pplumi=5300.;
  //float lumi=129.;
  //float pplumi=212.;
  //*************************************************************************
	
  // Output file
  TFile *pbpb_Unfo = new TFile(Form("pbpb_pp_Unfo_%s_cent_%d.root",algoName[algo],nbins_cent),"RECREATE");
	
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
  TF1 *fVz;
  TF1* fCentralityWeight;
  TCut dataSelection;
  TCut dataSelectionPbPb;
  TCut dataSelectionPP;
  TCut TriggerSelectionPP;
  TCut TriggerSelectionPbPb;

  if (isMC) {
    // MC closure test, no reweighting
    fVz = new TF1("fVz","1");
    fCentralityWeight =  new TF1("fCentralityWeight","1");
    dataSelection = "abs(vz)<15&&trackMax/jtpt>0.01&&abs(jteta)<2";
  } else {
    // Reweight to describe data
    fVz = new TF1("fVz","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
    fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05);
    fCentralityWeight = new TF1("fCentralityWeight","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
    fCentralityWeight->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04);
    //dataSelectionPbPb = "abs(vz)<15&&pcollisionEventSelection&&pHBHENoiseFilter&&abs(jteta)<2";
    //dataSelectionPP = "abs(vz)<15&&pPAcollisionEventSelectionPA&&pHBHENoiseFilter&&abs(jteta)<2";
    //TriggerSelectionPP = "HLT_PAJet80_NoJetID_v1";
    //TriggerSelectionPbPb ="HLT_HIJet80_v1";

    dataSelectionPbPb = "abs(eta)<2";
    dataSelectionPP = "abs(eta)<2";
		
    TriggerSelectionPP = "jet80";
    TriggerSelectionPbPb = "jet80";

  }
	
  // Vertex reweighting for pp
  TF1 *fVzPP = new TF1("fVzPP","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x");
  fVzPP->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05);
	
   	
  // Read data file
  TFile *infData;
	
  if (isMC) {
    if (yinglu) {
      infData = new TFile("/hadoop/store/user/belt/hiForest2/Pythia80_HydjetDrum_mix01_HiForest2_v20.root");
    } else {
		
    }
    cout << "This is a MC closure test"<<endl;
  } else {
    if (yinglu) {
      infData = new TFile("/hadoop/store/user/belt/hiForest2/promptskim-hihighpt-hltjet80-pt90-v20.root");
    } else {
      //infData = new TFile("/d102/yjlee/hiForest2/promptskim-hihighpt-hltjet80-pt90-v20.root");
      //infData = new TFile("/mnt/hadoo/cms/store/user/ymao/hiForest/PbPb2011/HiTracking");
      infData = new TFile("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/ntuple_2011_pbpbJet80.root");
    }
    cout << "This is a data analysis"<<endl;
  }

  TTree *tDataJet = (TTree*)infData->Get("ntjet");

  /*
    TTree *tDataEvt = (TTree*)infData->Get("hiEvtAnalyzer/HiTree");
    TTree *tDataSkim = (TTree*)infData->Get("skimanalysis/HltTree");
    TTree *tDataHlt = (TTree*)infData->Get("hltanalysis/HltTree");
    TTree *tDataJet  = (TTree*)infData->Get(Form("%sJetAnalyzer/t",algoName[algo]));
    tDataJet->AddFriend(tDataEvt);
    tDataJet->AddFriend(tDataSkim);
    tDataJet->AddFriend(tDataHlt);

    int dataNJet;
    tDataJet->SetBranchAddress("nref",&dataNJet);
  */

  Float_t dataNJet;
  tDataJet->SetBranchAddress("nrefe",&dataNJet);
	
  TFile *infPP;
	
  if (isMC) {
    if (yinglu) {
      infPP = new TFile("/hadoop/store/user/belt/hiForest2/pp276Dijet80_merged.root");
    } else {
		
    }
  } else {
    //infPP = new TFile("/hadoop/store/user/belt/hiForest2/HiForest-ppskim-hihighpt-pt90-v1_v3_part.root");
    if (yinglu) {
      infPP = new TFile("/hadoop/store/user/belt/hiForest2/pp_merged_full.root");
    } else {
      //infPP = new TFile("/d102/yjlee/hiForest2PP/pp_merged_full.root");
      infPP = new TFile("/mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/ntuple_2013_ppJet80.root");
    }   
  }

  TTree *tPPJet = (TTree*)infPP->Get("ntjet");

  /*
    TTree *tPPEvt = (TTree*)infPP->Get("hiEvtAnalyzer/HiTree");
    TTree *tPPSkim = (TTree*)infPP->Get("skimanalysis/HltTree");
    TTree *tPPHlt = (TTree*)infPP->Get("hltanalysis/HltTree");
    TTree *tPPJet  = (TTree*)infPP->Get(Form("%sJetAnalyzer/t",algoNamePP[algo]));
    tPPJet->AddFriend(tPPEvt);
    tPPJet->AddFriend(tPPSkim);
    tPPJet->AddFriend(tPPHlt);
  */
	
  // Setup jet data branches
  JetData *data[nbins_pthat]; 
  JetData *dataPP[nbins_pthat];
  JetData *data_pq; 
	
  for (int i=0;i<nbins_pthat;i++) data[i] = new JetData(fileName_pthat[i],Form("%sJetAnalyzer/t",algoName[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));	
  for (int i=0;i<nbinsPP_pthat;i++) dataPP[i] = new JetData(fileNamePP_pthat[i],Form("%sJetAnalyzer/t",algoNamePP[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));	
	
  if (yinglu) data_pq = new JetData(fileName_pthat_pq,Form("%sJetAnalyzer/t",algoName[algo]),Form("%sJetAnalyzer/t",algoNameGen[algo]));	
	
  TFile *fSpectra(0);
  TFile *akPu2PFResult(0);
  TFile *akPu3PFResult(0);
  TFile *akPu4PFResult(0);
  TFile *akPu2CaloResult(0);
  TFile *akPu3CaloResult(0);
  TFile *akPu4CaloResult(0);
	
	
  if (useSpectraFromFile||useMatrixFromFile){
    fSpectra = new TFile(spectraFileName,"read");
    akPu2PFResult = new TFile("result/pbpb_Unfo_akPu2PF.root");
    akPu3PFResult = new TFile("result/pbpb_Unfo_akPu3PF.root");
    akPu4PFResult = new TFile("result/pbpb_Unfo_akPu4PF.root");
    akPu2CaloResult = new TFile("result/pbpb_Unfo_akPu2Calo.root");
    akPu3CaloResult = new TFile("result/pbpb_Unfo_akPu3Calo.root");
    akPu4CaloResult = new TFile("result/pbpb_Unfo_akPu4Calo.root");
  }
	
  // Come back to the output file dir
  pbpb_Unfo->cd();
	
	
  Float_t rndGaus[1000];
	
  TTree *tRandom = new TTree("tRandom","");
  tRandom->Branch("nrefe",&dataNJet,"nrefe/F");
  tRandom->Branch("rndGaus",rndGaus,"rndGaus[nrefe]/F");
	
  TRandom rnd;
  /*
    if (!useSpectraFromFile) {
    for (int i=0;i<tDataJet->GetEntries();i++)
    {
    tDataJet->GetEntry(i);
    if (i%10000==0)  cout <<i<<" / "<<tDataJet->GetEntries()<<endl;
    for (int j=0;j<dataNJet;j++)
    {
    //rndGaus[j]=rnd.Gaus(0,1);
    rndGaus[j] = 1; //extremely stupid but just go with it for now. 
    }   
    tRandom->Fill();
    }
    }
  */
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
	
	
  if (!useSpectraFromFile) {
    tDataJet->AddFriend(tRandom);
    tDataJet->Project("hCentData","bin",dataSelectionPbPb);
    tDataJet->Project("hVzData","vz",dataSelectionPbPb);
  }
	
  for (int i=0;i<=nbins_cent;i++){
    //TCut centCut = Form("hiBin<%.0f&&hiBin>=%.0f",boundaries_cent[i+1],boundaries_cent[i]);
    TCut centCut = Form("bin<%.0f&&bin>=%.0f",boundaries_cent[i+1],boundaries_cent[i]);
    if (useSpectraFromFile) {
      uhist[i]->hMeas = (TH1F*)fSpectra->Get(Form("hMeas_cent%d",i));
      uhist[i]->hMeasJECSys = (TH1F*)fSpectra->Get(Form("hMeasJECSys_cent%d",i));
      uhist[i]->hMeasSmearSys = (TH1F*)fSpectra->Get(Form("hMeasSmearSys_cent%d",i));
    } else {
      if (!isMC) {
	tDataJet->Project(Form("hMeas_cent%d",i),"pt", dataSelectionPbPb&&centCut&&TriggerSelectionPbPb);
	tDataJet->Project(Form("hMeasJECSys_cent%d",i),Form("pt*%f",1.+0.02/nbins_cent*(nbins_cent-i)), dataSelectionPbPb&&centCut&&TriggerSelectionPbPb);
	//tDataJet->Project(Form("hMeasSmearSys_cent%d",i),"pt+rndGaus", dataSelectionPbPb&&centCut&&TriggerSelectionPbPb);
	tDataJet->Project(Form("hMeasSmearSys_cent%d",i),"pt", dataSelectionPbPb&&centCut&&TriggerSelectionPbPb);
	//tDataJet->Project(Form("hGen_cent%d",i),"refpt", dataSelection&&centCut);
	//removeError(uhist[i]->hGen);
      }   
    }
		
    if (useMatrixFromFile) {
      hCentMC = (TH1F*) fSpectra->Get("hCentMC");
      uhist[i]->hMatrix = (TH2F*) fSpectra->Get(Form("hMatrix_cent%d",i));
      uhist[i]->hMeasMatch = (TH1F*)((TH2F*) fSpectra->Get(Form("hMatrix_cent%d",i)))->ProjectionY();
      uhist[i]->hMeasMatch->Divide(uhist[i]->hMeas);
    } else {
      uhist[i]->hMeasMatch = 0;
    }
    uhist[i]->hMeas->Draw();
  }
	
  if (!isMC) 
    tPPJet->Project(Form("hMeas_cent%d",nbins_cent),"pt",dataSelectionPP&&TriggerSelectionPP);
  tPPJet->Project("hVzPPData","vz",dataSelectionPP&&TriggerSelectionPP);
	
	
  TCanvas *c = new TCanvas("c","",600,600);
  TH1F *hPtHat = new TH1F("hPtHat","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatRaw = new TH1F("hPtHatRaw","",nbins_pthat,boundaries_pthat);
  TH1F *hPtHatPP = new TH1F("hPtHatPP","",nbinsPP_pthat,boundariesPP_pthat);
  TH1F *hPtHatRawPP = new TH1F("hPtHatRawPP","",nbinsPP_pthat,boundariesPP_pthat);
	
  RooUnfoldResponse res(uhist[0]->hResMeas,uhist[0]->hResTrue);
	
  for (int i=0;i<nbins_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbins_pthat,boundaries_pthat);
    data[i]->tJet->Project("hPtHatTmp","pthat");
    hPtHatRaw->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
	
  for (int i=0;i<nbinsPP_pthat;i++) {
    TH1F *hPtHatTmp = new TH1F("hPtHatTmp","",nbinsPP_pthat,boundariesPP_pthat);
    dataPP[i]->tJet->Project("hPtHatTmp","pthat","abs(vz)<15");
    hPtHatRawPP->Add(hPtHatTmp);
    delete hPtHatTmp;
  }
	
  // Fill PbPb MC   
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
	    uhist[cBin]-> hGen->Fill(data[i]->refpt[k],scale*weight_cent*weight_pt*weight_vz);
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
		
	weight_cent = fCentralityWeight->Eval(data_pq->bin);
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
	if(dataPP[i]->bin<=28) continue;
	int pthatBin = hPtHatPP->FindBin(dataPP[i]->pthat);
	float scale = (xsectionPP[pthatBin-1]-xsectionPP[pthatBin])/hPtHatRawPP->GetBinContent(pthatBin);
	if(fabs(dataPP[i]->vz)>15) continue;
	double weight_cent=1;
	double weight_pt=1;
	double weight_vz=1;
				
	weight_vz = fVzPP->Eval(dataPP[i]->vz);
	if (weight_vz>5||weight_vz<0.5) cout <<dataPP[i]->vz<<" "<<weight_vz<<endl;
	//weight_vz = 1;
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
	    uhist[nbins_cent]-> hGen->Fill(dataPP[i]->refpt[k],scale*weight_vz);   

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
  cMC->Divide(3,3);
  TCanvas *cData = new TCanvas("cData","Data",1000,800);
  cData->Divide(3,3);
  for(int i = 0;i<=nbins_cent;i++){

    TLegend *title1 = 0,*title2 = 0;

    cMC->cd(i+1);
    cMC->cd(i+1)->SetLogy();
    title1 = myLegend(0.18,0.35,0.48,0.45);//MC
    title1->SetTextSize(0.06);

    uhist[i]->hGen->Draw();
    if(i == 6){
      title1->AddEntry(uhist[i]->hGen,"PP MC","");
    }else{
      title1->AddEntry(uhist[i]->hGen,Form("PbPb MC - %2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    }
    title1->Draw();

    cData->cd(i+1);
    cData->cd(i+1)->SetLogy();
    uhist[i]->hMeas->Draw();

    title2 = myLegend(0.18,0.7,0.48,0.8);//data
    title2->SetTextSize(0.06);
    if(i == 6){
      title2->AddEntry(uhist[i]->hMeas,"PP Data","");
    }else {
      title2->AddEntry(uhist[i]->hMeas,Form("PbPb Data - %2.0f-%2.0f%%",2.5*boundaries_cent[i],2.5*boundaries_cent[i+1]),"");
    }
    title2->Draw();
  }

	
  TCanvas *cCent = new TCanvas("cCent","Centrality",600,600);
  divideBinWidth(hCentData);
  hCentMC->Scale(1./hCentMC->Integral(0,1000));
  hCentData->Scale(1./hCentData->Integral(0,1000));
  hCentMC->SetMarkerColor(2);
  hCentMC->SetLineColor(2);
  hCentData->Draw();
  hCentMC->Draw("same");
	
  cout <<"Response Matrix..."<<endl;
	
  TCanvas * cMatrix = new TCanvas("cMatrix","Matrix",1200,800);
  cMatrix->Divide(3,3);
  TCanvas* cResponseNorm = new TCanvas("cResponseNorm","Normalized Response Matrix",1200,800);
  cResponseNorm->Divide(3,3);
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
  hPtHat->Draw();
	
  cout << "==================================== TEST =====================================" << endl;
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  char chmet[100]; 
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
  TCanvas * cPbPb = new TCanvas("cPbPb","PbPb",1200,800);
  cPbPb->Divide(3,3); 
  cPbPb->cd(1);
	
	
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

    bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrix,hPrior,0);
    myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);

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
	

  cPbPb->cd(8);
  TLegend *pbpblegend = myLegend(0.52,0.65,0.85,0.9);
  uhist[0]->hMeas->SetMarkerStyle(20);
  uhist[0]->hMeas->SetMarkerColor(kRed);
  uhist[0]->hReco->SetMarkerStyle(25);
  pbpblegend->AddEntry(uhist[0]->hMeas,"Measured","pl");
  pbpblegend->AddEntry(uhist[0]->hReco,"Unfolded","pl");
  pbpblegend->Draw();
	
	
	
	
	
	
	
  //***************         Calculation of RAA        ************
	
  SysData systematics;
	
	
	
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
	

  TCanvas * cRAA = new TCanvas("cRAA","RAA",1200,800);
  makeMultiPanelCanvasWithGap(cRAA,3,2,0.01,0.01,0.16,0.2,0.04,0.04);

	
  TCanvas * cSys = new TCanvas("cSys","Total Systematics",1200,800);
  makeMultiPanelCanvasWithGap(cSys,3,2,0.01,0.01,0.16,0.2,0.04,0.04);
	
  TCanvas * cIterSysPP = new TCanvas("cIterSysPP","Iteration Systematics for PP",800,600);
	
  TCanvas * cPPMCclosure = new TCanvas("cPPMCclosure","cPPMCclosure",600,450);
 
	
	
  //************   MC Closure Test For PP        *************
	
	
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
	
	
  //************    Iteration Sys PP  ak3PF        *************
	 
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
	 
	 
	 
	 
	
	
  // calculate the maximum from all cent bins
  for (int i=0;i<nbins_cent;i++) {
		
    TH1F *hRecoRAAIterSys[10];
    TH1F *hRebinRAA = rebin(uhist[i]->hReco, Form("hRebinRAA_cent%d",i));
		
    for (int j=2;j<7;j++) {
      hRecoRAAIterSys[j] = rebin(uhist[i]->hRecoIterSys[j],Form("hRecoRAA_IterSys%d_cent%d",j,i));
      hRecoRAAIterSys[j]->SetLineColor(colorCode[j-2]);
      hRecoRAAIterSys[j]->SetMarkerColor(colorCode[j-2]);
      hRecoRAAIterSys[j]->Divide(hRebinRAA);
      checkMaximumSys(systematics.hSysIter[nbins_cent],hRecoRAAIterSys[j],0,1.1);
    }
		
		
  }
  for (int i=0;i<nbins_cent;i++) {
    cRAA->cd(nbins_cent-i);
		
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
		
    // Iteration systematics
    TH1F *hRecoRAAIterSys[10];
    TLegend *legBayesianIter = myLegend(0.6,0.6,0.9,0.9);
    legBayesianIter->SetTextSize(0.042);
		
    for (int j=2;j<7;j++) {
			
			
      //************    Iteration Sys   akPu3PF        *************
			
      cout <<" Plotting  Iteration Sys   akPu3PF   " << endl;
      cIterSys->cd(nbins_cent-i);
			
      hRecoRAAIterSys[j] = rebin(uhist[i]->hRecoIterSys[j],Form("hRecoRAA_IterSys%d_cent%d",j,i));
      hRecoRAAIterSys[j]->SetLineColor(colorCode[j-2]);
      hRecoRAAIterSys[j]->SetMarkerColor(colorCode[j-2]);
      hRecoRAAIterSys[j]->Divide(hRebinRAA);
      if (j==2){
	makeHistTitle(hRecoRAAIterSys[j],"","Jet p_{T} (GeV/c)","Ratio (Unfolded / Nominal)");
	hRecoRAAIterSys[j]->SetAxisRange(0,2,"Y");
	hRecoRAAIterSys[j]->Draw(); 
      } else {
	hRecoRAAIterSys[j]->Draw("same");
      }
			
      checkMaximumSys(systematics.hSysIter[i],hRecoRAAIterSys[j],0,1.1);
      legBayesianIter->AddEntry(hRecoRAAIterSys[j],Form("Iteration %d",j),"pl"); 	  
    }
    if (i==nbins_cent-1) legBayesianIter->Draw();
    title->Draw();
    l->Draw();
		
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
		
	
		
		
    //************    Different methods RAA in akPu3PF        *************
		
		
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
		
  }   // centrality for loop end
	
	
	
	
  //************   Test  Plotting ends      ***************
	
	
	
  if (isMC) {
    cRAA->Update();
    cRAA->SaveAs("MCClosureTest/MCClosureTest.gif"); 
    cRAA->SaveAs("MCClosureTest/MCClosureTest.C"); 
    cRAA->SaveAs("MCClosureTest/MCClosureTest.pdf"); 
    cIterSys->Update();
    cIterSys->SaveAs("MCClosureTest/IterSys.gif");
    cIterSys->SaveAs("MCClosureTest/IterSys.C");
    cIterSys->SaveAs("MCClosureTest/IterSys.pdf");
    cPPMCclosure->Update();
    cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.gif");
    cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.C");
    cPPMCclosure->SaveAs("MCClosureTest/MCClosureTestPP.pdf");
	
		
  } else {
	
	
    cRAA->Update();
    cRAA->SaveAs(Form("result-%s-cent-%d/result.gif",algoName[algo],nbins_cent),"RECREATE"); 
    cRAA->SaveAs(Form("result-%s-cent-%d/result.C",algoName[algo],nbins_cent),"RECREATE"); 
    cRAA->SaveAs(Form("result-%s-cent-%d/result.pdf",algoName[algo],nbins_cent),"RECREATE"); 
    cPbPb->Update();
    cPbPb->SaveAs(Form("result-%s-cent-%d/data_vs_unfo.gif",algoName[algo],nbins_cent),"RECREATE");
    cPbPb->SaveAs(Form("result-%s-cent-%d/data_vs_unfo.C",algoName[algo],nbins_cent),"RECREATE");
    cPbPb->SaveAs(Form("result-%s-cent-%d/data_vs_unfo.pdf",algoName[algo],nbins_cent),"RECREATE");
    cIterSys->Update();
    cIterSys->SaveAs(Form("result-%s-cent-%d/IterSys.gif",algoName[algo],nbins_cent),"RECREATE");
    cIterSys->SaveAs(Form("result-%s-cent-%d/IterSys.C",algoName[algo],nbins_cent),"RECREATE");
    cIterSys->SaveAs(Form("result-%s-cent-%d/IterSys.pdf",algoName[algo],nbins_cent),"RECREATE");	
    cIterSysPP->Update();
    cIterSysPP->SaveAs(Form("result-%s-cent-%d/IterSysPP.gif",algoName[algo],nbins_cent),"RECREATE");
    cIterSysPP->SaveAs(Form("result-%s-cent-%d/IterSysPP.C",algoName[algo],nbins_cent),"RECREATE");
    cIterSysPP->SaveAs(Form("result-%s-cent-%d/IterSysPP.pdf",algoName[algo],nbins_cent),"RECREATE");
    cSys->Update();
    cSys->SaveAs(Form("result-%s-cent-%d/TotalSys.gif",algoName[algo],nbins_cent),"RECREATE");
    cSys->SaveAs(Form("result-%s-cent-%d/TotalSys.C",algoName[algo],nbins_cent),"RECREATE");
    cSys->SaveAs(Form("result-%s-cent-%d/TotalSys.pdf",algoName[algo],nbins_cent),"RECREATE");
    cSmearSys->Update();
    cSmearSys->SaveAs(Form("result-%s-cent-%d/SmearSys.gif",algoName[algo],nbins_cent),"RECREATE"); 
    cSmearSys->SaveAs(Form("result-%s-cent-%d/SmearSys.C",algoName[algo],nbins_cent),"RECREATE"); 
    cSmearSys->SaveAs(Form("result-%s-cent-%d/SmearSys.pdf",algoName[algo],nbins_cent),"RECREATE"); 
    cJECSys->Update();
    cJECSys->SaveAs(Form("result-%s-cent-%d/JECSys.gif",algoName[algo],nbins_cent),"RECREATE");
    cJECSys->SaveAs(Form("result-%s-cent-%d/JECSys.C",algoName[algo],nbins_cent),"RECREATE");
    cJECSys->SaveAs(Form("result-%s-cent-%d/JECSys.pdf",algoName[algo],nbins_cent),"RECREATE");
    cMC->Update();
    cMC->SaveAs(Form("result-%s-cent-%d/MC.gif",algoName[algo],nbins_cent),"RECREATE");
    cMC->SaveAs(Form("result-%s-cent-%d/MC.C",algoName[algo],nbins_cent),"RECREATE");
    cMC->SaveAs(Form("result-%s-cent-%d/MC.pdf",algoName[algo],nbins_cent),"RECREATE");		
    cData->Update();
    cData->SaveAs(Form("result-%s-cent-%d/Data.gif",algoName[algo],nbins_cent),"RECREATE");
    cData->SaveAs(Form("result-%s-cent-%d/Data.C",algoName[algo],nbins_cent),"RECREATE");
    cData->SaveAs(Form("result-%s-cent-%d/Data.pdf",algoName[algo],nbins_cent),"RECREATE");
    cMatrix->Update();
    cMatrix->SaveAs(Form("result-%s-cent-%d/ResponseMatrix.gif",algoName[algo],nbins_cent),"RECREATE");
    cMatrix->SaveAs(Form("result-%s-cent-%d/ResponseMatrix.C",algoName[algo],nbins_cent),"RECREATE");
    cMatrix->SaveAs(Form("result-%s-cent-%d/ResponseMatrix.pdf",algoName[algo],nbins_cent),"RECREATE");
    cResponseNorm->Update();
    cResponseNorm->SaveAs(Form("result-%s-cent-%d/NormalizedResponseMatrix.gif",algoName[algo],nbins_cent),"RECREATE");
    cResponseNorm->SaveAs(Form("result-%s-cent-%d/NormalizedResponseMatrix.C",algoName[algo],nbins_cent),"RECREATE");
    cResponseNorm->SaveAs(Form("result-%s-cent-%d/NormalizedResponseMatrix.pdf",algoName[algo],nbins_cent),"RECREATE");
		
		
  }
  divideBinWidth(hCent);
	
	
	
  pbpb_Unfo->Write();
  //pbpb_Unfo->Close();

  timer.Stop();
  cout<<" Total time taken CPU = "<<timer.CpuTime()<<endl;
  cout<<" Total time taken Real = "<<timer.RealTime()<<endl;

}






