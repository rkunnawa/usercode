#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TChain.h"
#include "TF1.h"
#include "TMath.h"

//These includes cause some complications in CMSSW_5_3_8_HI_patchX. Commented out for pp.
//If you want to recalculate the JECs on the fly again, just uncomment everything in the updateJEC blocks

//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

// ******* GLOBAL DECLARATIONS **********
//const int QCDpthatBins = 8;
const int QCDpthatBins = 9;
//const int QCDpthatBins = 1;
const int dataFiles = 6;

const TString algo = "ak3PF" ;
const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
const Double_t jetPtBin[]={3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;
const int pthatbin[10] = {15,30,50,80,120,170,220,280,370, 9999};
const double wght[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};

//***************************************


//**********************************************************
// Count the MC events to appropriately weight the pthat bins
//**********************************************************

int *countMCevents(std::string infile, bool usePUsub, int isMC){

  TChain *ch = NULL;
//  if(usePUsub) 
    ch = new TChain(Form("%sJetAnalyzer/t", algo.Data()));
//  else ch = new TChain("ak3PFJetAnalyzer/t");
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  int *MCentries = new int[QCDpthatBins];
  for(int ifile=0; ifile<QCDpthatBins; ifile++){
    instr >> filename;
    ch->Add(filename.c_str());
    if(ifile==0) MCentries[ifile] = ch->GetEntries();
    else {
      int sum = 0 ;
      for(int i = 0 ; i <ifile ; i++) sum+=MCentries[i];
     MCentries[ifile] = ch->GetEntries()-sum; 
     }
 //     MCentries[ifile] = ch->GetEntries(Form("pthat<%d",pthatbin[ifile+1]));
 //   MCentries[ifile] = ch->GetEntries();
  }
//  MCentries[0] = ch->GetEntries("pthat<15");
/*  MCentries[0] = ch->GetEntries("pthat>=15 && pthat<30");
  MCentries[1] = ch->GetEntries("pthat>=30 && pthat<50");
  MCentries[2] = ch->GetEntries("pthat>=50 && pthat<80");
  MCentries[3] = ch->GetEntries("pthat>=80 && pthat<120");
  MCentries[4] = ch->GetEntries("pthat>=120 && pthat<170");
  MCentries[5] = ch->GetEntries("pthat>=170 && pthat<220");
  MCentries[6] = ch->GetEntries("pthat>=220 && pthat<280");

  MCentries[0] = ch->GetEntries("pthat<30");
  MCentries[1] = ch->GetEntries("pthat<50");
  MCentries[2] = ch->GetEntries("pthat<80");
  MCentries[3] = ch->GetEntries("pthat<120");
  MCentries[4] = ch->GetEntries("pthat<170");
  MCentries[5] = ch->GetEntries("pthat<220");
  MCentries[6] = ch->GetEntries("pthat<280");
  if(QCDpthatBins==8)MCentries[7] = ch->GetEntries("pthat>=280");
  else {
   MCentries[7] = ch->GetEntries("pthat>=280 && pthat<370");
//   MCentries[7] = ch->GetEntries("pthat<370");
   MCentries[8] = ch->GetEntries("pthat>=370");
  }
*/  for(int i=0; i<QCDpthatBins; i++){
    cout << "QCD MCentries[" << i << "]: " << MCentries[i] << endl;
  }


  
  return MCentries;
}

//**********************************************************
// Trigger-Combine the data in order to unfold properly later
// need to adjust here to use the leading jet binning combination
//**********************************************************

//[0] = MB, [1] = Jet20, [2] = Jet40, [3] = Jet60, [4] = Jet80, [5] = Jet100
double trigComb(bool *triggerDecision, double *pscl){
  double weight=1;
  // if(triggerDecision[0] && !triggerDecision[1] && !triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[0]); //Removing finnicky Jet20 sample
  if(triggerDecision[1] && !triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[1]);
  if(triggerDecision[2] && !triggerDecision[3]) weight = 1./(1./pscl[1] + 1./pscl[2] - (1./(pscl[1]*pscl[2])));
  if(triggerDecision[3]) weight = 1.;
  return weight;
}

//**********************************************************
// "get" the trigger prescales by counting trigger overlap
//**********************************************************

double* getPscls(std::string infile, int nFiles, bool usePUsub){
  
  TChain *dataCH = NULL;
//  if(usePUsub){
    dataCH = new TChain(Form("%sJetAnalyzer/t", algo.Data()));
//  }
//  else dataCH = new TChain("ak3PFJetAnalyzer/t");
  TChain *dataCH2 = new TChain("hltanalysis/HltTree");
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  for(int ifile=0; ifile<nFiles; ifile++){
    instr >> filename;
    dataCH->Add(filename.c_str());
    dataCH2->Add(filename.c_str());
  }
  dataCH->AddFriend(dataCH2, "hltanalysis/HltTree");
//    TFile *fin=NULL;
//    for(int ifile=0; ifile<nFiles; ifile++){
//        fin = TFile::Open(filename.c_str(), "readonly");
//        Long64_t nentries = dataCH->GetEntries();
//        for (Long64_t i=0; i<nentries;i++) {
//            dataCH->GetEntry(i);
//
//        }
//    }
//    bool jet40 = false ;
//    bool jet60 = false ;
//    bool jet75 = false ;
//    bool jet95 = false ;
//    bool jet120 = false ;
//    bool combined = false ;
    
  //Set up trigger combination prescales for data
  double ov1, ov2, ov3, ov4;
  ov1 = dataCH->GetEntries("HLT_PAJet20_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov2 = dataCH->GetEntries("HLT_PAJet40_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov3 = dataCH->GetEntries("HLT_PAJet60_NoJetID_v1 && HLT_PAJet80_NoJetID_v1");
  ov4 = dataCH->GetEntries("HLT_PAJet80_NoJetID_v1");
  double *pscls = new double[4];
  pscls[0] = ov4/ov1;
  pscls[1] = ov4/ov2;
  pscls[2] = ov4/ov3;
  pscls[3] = 1.;
  return pscls;
}

//**********************************************************
// ~~~ MAIN PROGRAM ~~~
//**********************************************************

void analyzeTreesPA(int isRecopp=1, int ppPbPb=0, int isMC=1, int doNtuples=1, int doJets=1, int doTracks=1, int updateJEC=0, int cbin=-1, bool ExpandedTree=false, bool usePUsub=0)
{
  // isMC=0 --> Real data, ==1 --> QCD, ==2 --> bJet, ==3 --> cJet
  Float_t minJetPt=0.;
  
  Float_t maxJetEta=3.;
  Float_t minMuPt=5;
  
  // cbin = -1 --> 0-100%
  // cbin = 0 --> 0-20%
  // cbin = 1 --> 20-50%
  // cbin =2 --> 50-100%
  if(!ppPbPb) cbin=-1;
  int useWeight=1;

  //Get vz weight histograms for MC if needed
//  TH1D *hMCvz[2], hDatavz;
//  TFile *fMCvz, *fDatavz, *fBvz;
//  if(isMC){
      TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
      fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
    
//    fMCvz = new TFile("MCvzDistr.root");
//    hMCvz[0] = (TH1D*)(fMCvz->Get("hvz"))->Clone("hMCvz_0");
//    fDatavz = new TFile("DatavzDistro.root");
//    hDatavz = (TH1D*)(fMCvz->Get("hvzData"))->Clone("hDatavz");
//    fBvz = new TFile("BvzDistr.root");
//    hMCvz[1] = (TH1D*)(fBvz->Get("hvzB"))->Clone("hMCvz_1");
//    
//    hMCvz[0]->Scale(1./hMCvz[0]->Integral());
//    hMCvz[1]->Scale(1./hMCvz[1]->Integral());
//    hDatavz->Scale(1./hDatavz->Integral());
//  }

  cout << "Analyzing Trees! Assuming " << QCDpthatBins  << endl;

  int pthatbin[10] = {15,30,50,80,120,170,220,280,370, 9999};
  double wght[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};
/*  if(QCDpthatBins==8) {
    pthatbin[]= {15,30,50,80,120,170,220,280, 9999};
    wght[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 0};
  }
 else {
   pthatbin[]= {15,30,50,80,120,170,220,280,370, 9999};
    wght[]={5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05,1.018E-05,2.477E-06,6.160E-07, 1.088E-07, 0};
  } 
*/ //  int pthatbin[QCDpthatBins+1] = {80, 9999};
  double w = 1.;
//  double wght[QCDpthatBins+1]={4.412E-04, 0};

  TFile *fin=NULL;
  std::string infile;
  int *MCentr = NULL;
  double *pscls = NULL;

  //File load
    if(!isMC){
      infile = "pPbDataForest.txt";
    }
    else {
    //  infile = "pPbForestList.txt";
      infile = "PbpForestList.txt";
    //  infile = "pPbForestListProd25.txt";
    //  infile = "ppForestList.txt";
    //  infile = "test.txt";
    }

  if(!isMC){
   // pscls = getPscls(infile,QCDpthatBins,usePUsub);
  }
  
  int dupRuns[6] = {181912,181913,181938,181950,181985,182124};
  
  std::vector<int> usedEvents[6];
  int nDup=0;

  //Declaration of leaves types
  Int_t evt;
  Int_t run;
  Int_t bin;
  Int_t lumi;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
  Float_t hf;
  Float_t vz;
  Int_t nref;
  Float_t rawpt[1000];
  Float_t jtpt[1000];
  Float_t jteta[1000];
  Float_t jty[1000];
  Float_t jtphi[1000];
  Float_t jtpu[1000];
  Float_t discr_ssvHighEff[1000];
  Float_t discr_ssvHighPur[1000];
  //Float_t discr_csvMva[1000];
  Float_t discr_csvSimple[1000];
  // Float_t discr_muByIp3[1000];
  Float_t discr_muByPt[1000];
  Float_t discr_prob[1000];
  Float_t discr_probb[1000];
  Float_t discr_tcHighEff[1000];
  Float_t discr_tcHighPur[1000];
  Int_t nsvtx[1000];
  Int_t svtxntrk[1000];
  Float_t svtxdl[1000];
  Float_t svtxdls[1000];
  Float_t svtxm[1000];
  Float_t svtxpt[1000];
  Int_t nIPtrk[1000];
  Int_t nselIPtrk[1000];
  Int_t nIP;
  Int_t ipJetIndex[10000];
  Float_t ipPt[10000];
  Float_t ipProb0[10000];
  Float_t ipProb1[10000];
  Float_t ip2d[10000];
  Float_t ip2dSig[10000];
  Float_t ip3d[10000];
  Float_t ip3dSig[10000];
  Float_t ipDist2Jet[10000];
  Float_t ipDist2JetSig[10000];
  Float_t ipClosest2Jet[10000];
  //Float_t mue[1000];
  Float_t mupt[1000];
  Float_t muptPF[1000];
  Float_t mueta[1000];
  Float_t muphi[1000];
  //Float_t mudr[1000];
  Float_t muptrel[1000];
  //Int_t muchg[1000];
  Float_t pthat;
  Int_t beamId1;
  Int_t beamId2;
  Int_t subid[1000];
  Float_t refpt[1000];
  Float_t refeta[1000];
  Float_t refy[1000];
  Float_t refphi[1000];
  Float_t refdphijt[1000];
  Float_t refdrjt[1000];
  Float_t refparton_pt[1000];
  Int_t refparton_flavor[1000];
  Int_t refparton_flavorForB[1000];
  Bool_t refparton_isGSP[1000];
  Int_t          HLT_PAZeroBiasPixel_SingleTrack_v1;
  Int_t          HLT_PAJet20_NoJetID_v1;
  Int_t          HLT_PAJet40_NoJetID_v1;
  Int_t          HLT_PAJet60_NoJetID_v1;
  Int_t          HLT_PAJet80_NoJetID_v1;
  Int_t          HLT_PAJet100_NoJetID_v1;
  Int_t pVertexFilterCutGplusUpsPP;
  Int_t pPAcollisionEventSelectionPA;
  Int_t pHBHENoiseFilter;
  Int_t pprimaryvertexFilter;

  /*
Int_t ngen;
Int_t genmatchindex[1000];
Float_t genpt[1000];
Float_t geneta[1000];
Float_t geny[1000];
Float_t genphi[1000];
Float_t gendphijt[1000];
Float_t gendrjt[1000];
*/

  float chargedMax[1000];
  float photonMax[1000];
  float neutralMax[1000];
  float chargedSum[1000];
  float photonSum[1000];
  float neutralSum[1000];
  float muSum[1000];
  float eSum[1000];

  Double_t weight, xSecWeight, centWeight, vzWeight;
 
  centWeight = 1. ; 
  int nHLTBit;
  bool hltBit[12];

  Int_t pvSel;
  Int_t hbheNoiseSel;
  Int_t spikeSel;
  Int_t collSell;

  TFile *fout=NULL;

   if ( isRecopp ) { // pp reco, jet triggered
      if(isMC) fout=new TFile(Form("histos/PbpMCProd24_ppReco_%s_QCDjetTrig_JetPt%.fnoIPupperCut.root",algo.Data(), minJetPt),"recreate");
      else fout=new TFile(Form("histos/ppdata_ppReco_%s_jetTrig_noIPupperCut.root", algo.Data()),"recreate");
    }
    else if( isRecopp && !usePUsub){
      if(isMC) fout=new TFile("histos/ppMC_ppReco_ak3PF_gsp2_QCDjetTrig_noIPupperCut.root","recreate");
      else fout=new TFile("histos/ppdata_ppReco_ak3PF_gsp2_jetTrig_noIPupperCut.root","recreate");
    }
    else if (!isRecopp) { // hi reco, jet triggered
      if(isMC)fout=new TFile("histos/ppMC_hiReco_jetTrig_addGSP_up.root","recreate");
      else fout=new TFile("histos/ppdata_hiReco_jetTrig_regPFforJets.root","recreate");
    }


  TH1D *hbin = new TH1D("hbin","hbin",100,0.5,100.5);
  TH1D *hbinw = new TH1D("hbinw","hbinw",100,0.5,100.5);
  hbin->Sumw2(); hbinw->Sumw2();

  TH1D *hvz = new TH1D("hvz","hvz",120,-15.,15.);
  TH1D *hvzw = new TH1D("hvzw","hvzw",120,-15.,15.);
  hvz->Sumw2(); hvzw->Sumw2();
  TH1F * hJetPtEtaBin[netabin];
  for(Int_t ieta = 0 ; ieta < netabin ; ieta++){
     hJetPtEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), nJetPtBin, jetPtBin);
     hJetPtEtaBin[ieta]->Sumw2();
   }

//  TH1D *hjtpt = new TH1D("hjtpt","hjtpt",1000,0,1000);
  TH1D *hjtpt = new TH1D("hjtpt","hjtpt", nJetPtBin, jetPtBin);

  TH1D *hjtptB = new TH1D("hjtptB","hjtptB",68,80,330);
  TH1D *hjtptC = new TH1D("hjtptC","hjtptC",68,80,330);
  TH1D *hjtptL = new TH1D("hjtptL","hjtptL",68,80,330);
  TH1D *hjtptU = new TH1D("hjtptU","hjtptU",68,80,330);
  hjtpt->Sumw2(); hjtptB->Sumw2(); hjtptC->Sumw2(); hjtptL->Sumw2(); hjtptU->Sumw2();

  TH1D *hrawpt = new TH1D("hrawpt","hrawpt",500,0,500);
  TH1D *hrawptB = new TH1D("hrawptB","hrawptB",68,80,330);
  TH1D *hrawptC = new TH1D("hrawptC","hrawptC",68,80,330);
  TH1D *hrawptL = new TH1D("hrawptL","hrawptL",68,80,330);
  hrawpt->Sumw2(); hrawptB->Sumw2(); hrawptC->Sumw2(); hrawptL->Sumw2();

  TH1D *hjteta = new TH1D("hjteta","hjteta",40,-3,3);
  TH1D *hjtetaB = new TH1D("hjtetaB","hjtetaB",40,-2,2);
  TH1D *hjtetaC = new TH1D("hjtetaC","hjtetaC",40,-2,2);
  TH1D *hjtetaL = new TH1D("hjtetaL","hjtetaL",40,-2,2);
  hjteta->Sumw2(); hjtetaB->Sumw2(); hjtetaC->Sumw2(); hjtetaL->Sumw2();

  TH1D *hjtphi = new TH1D("hjtphi","hjtphi",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiB = new TH1D("hjtphiB","hjtphiB",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiC = new TH1D("hjtphiC","hjtphiC",40,-1.*acos(-1.),acos(-1.));
  TH1D *hjtphiL = new TH1D("hjtphiL","hjtphiL",40,-1.*acos(-1.),acos(-1.));
  hjtphi->Sumw2(); hjtphiB->Sumw2(); hjtphiC->Sumw2(); hjtphiL->Sumw2();

  TH1D *hdiscr_csvSimple = new TH1D("hdiscr_csvSimple","hdiscr_csvSimple",25,0,1);
  TH1D *hdiscr_csvSimpleB = new TH1D("hdiscr_csvSimpleB","hdiscr_csvSimpleB",25,0,1);
  TH1D *hdiscr_csvSimpleC = new TH1D("hdiscr_csvSimpleC","hdiscr_csvSimpleC",25,0,1);
  TH1D *hdiscr_csvSimpleL = new TH1D("hdiscr_csvSimpleL","hdiscr_csvSimpleL",25,0,1);
  hdiscr_csvSimple->Sumw2(); hdiscr_csvSimpleB->Sumw2(); hdiscr_csvSimpleC->Sumw2(); hdiscr_csvSimpleL->Sumw2();
  
  TH1D *hdiscr_prob = new TH1D("hdiscr_prob","hdiscr_prob",25,0,2.5);
  TH1D *hdiscr_probB = new TH1D("hdiscr_probB","hdiscr_probB",25,0,2.5);
  TH1D *hdiscr_probC = new TH1D("hdiscr_probC","hdiscr_probC",25,0,2.5);
  TH1D *hdiscr_probL = new TH1D("hdiscr_probL","hdiscr_probL",25,0,2.5);
  hdiscr_prob->Sumw2(); hdiscr_probB->Sumw2(); hdiscr_probC->Sumw2(); hdiscr_probL->Sumw2();
  
  TH1D *hdiscr_ssvHighEff = new TH1D("hdiscr_ssvHighEff","hdiscr_ssvHighEff",25,1,6);
  TH1D *hdiscr_ssvHighEffB = new TH1D("hdiscr_ssvHighEffB","hdiscr_ssvHighEffB",25,1,6);
  TH1D *hdiscr_ssvHighEffC = new TH1D("hdiscr_ssvHighEffC","hdiscr_ssvHighEffC",25,1,6);
  TH1D *hdiscr_ssvHighEffL = new TH1D("hdiscr_ssvHighEffL","hdiscr_ssvHighEffL",25,1,6);
  hdiscr_ssvHighEff->Sumw2(); hdiscr_ssvHighEffB->Sumw2(); hdiscr_ssvHighEffC->Sumw2(); hdiscr_ssvHighEffL->Sumw2();
  
  TH1D *hdiscr_ssvHighPur = new TH1D("hdiscr_ssvHighPur","hdiscr_ssvHighPur",25,1,6);
  TH1D *hdiscr_ssvHighPurB = new TH1D("hdiscr_ssvHighPurB","hdiscr_ssvHighPurB",25,1,6);
  TH1D *hdiscr_ssvHighPurC = new TH1D("hdiscr_ssvHighPurC","hdiscr_ssvHighPurC",25,1,6);
  TH1D *hdiscr_ssvHighPurL = new TH1D("hdiscr_ssvHighPurL","hdiscr_ssvHighPurL",25,1,6);
  hdiscr_ssvHighPur->Sumw2(); hdiscr_ssvHighPurB->Sumw2(); hdiscr_ssvHighPurC->Sumw2(); hdiscr_ssvHighPurL->Sumw2();

  TH1D *hdiscr_tcHighEff = new TH1D("hdiscr_tcHighEff","hdiscr_tcHighEff",25,1,6);
  TH1D *hdiscr_tcHighEffB = new TH1D("hdiscr_tcHighEffB","hdiscr_tcHighEffB",25,1,6);
  TH1D *hdiscr_tcHighEffC = new TH1D("hdiscr_tcHighEffC","hdiscr_tcHighEffC",25,1,6);
  TH1D *hdiscr_tcHighEffL = new TH1D("hdiscr_tcHighEffL","hdiscr_tcHighEffL",25,1,6);
  hdiscr_tcHighEff->Sumw2(); hdiscr_tcHighEffB->Sumw2(); hdiscr_tcHighEffC->Sumw2(); hdiscr_tcHighEffL->Sumw2();
  
  TH1D *hdiscr_tcHighPur = new TH1D("hdiscr_tcHighPur","hdiscr_tcHighPur",25,1,6);
  TH1D *hdiscr_tcHighPurB = new TH1D("hdiscr_tcHighPurB","hdiscr_tcHighPurB",25,1,6);
  TH1D *hdiscr_tcHighPurC = new TH1D("hdiscr_tcHighPurC","hdiscr_tcHighPurC",25,1,6);
  TH1D *hdiscr_tcHighPurL = new TH1D("hdiscr_tcHighPurL","hdiscr_tcHighPurL",25,1,6);
  hdiscr_tcHighPur->Sumw2(); hdiscr_tcHighPurB->Sumw2(); hdiscr_tcHighPurC->Sumw2(); hdiscr_tcHighPurL->Sumw2();
  
  TH1D *hnsvtx = new TH1D("hnsvtx","hnsvtx",6,-0.5,5.5);
  TH1D *hnsvtxB = new TH1D("hnsvtxB","hnsvtxB",6,-0.5,5.5);
  TH1D *hnsvtxC = new TH1D("hnsvtxC","hnsvtxC",6,-0.5,5.5);
  TH1D *hnsvtxL = new TH1D("hnsvtxL","hnsvtxL",6,-0.5,5.5);
  hnsvtx->Sumw2(); hnsvtxB->Sumw2(); hnsvtxC->Sumw2(); hnsvtxL->Sumw2();
  
  TH1D *hsvtxntrk = new TH1D("hsvtxntrk","hsvtxntrk",12,-0.5,11.5);
  TH1D *hsvtxntrkB = new TH1D("hsvtxntrkB","hsvtxntrkB",12,-0.5,11.5);
  TH1D *hsvtxntrkC = new TH1D("hsvtxntrkC","hsvtxntrkC",12,-0.5,11.5);
  TH1D *hsvtxntrkL = new TH1D("hsvtxntrkL","hsvtxntrkL",12,-0.5,11.5);
  hsvtxntrk->Sumw2(); hsvtxntrkB->Sumw2(); hsvtxntrkC->Sumw2(); hsvtxntrkL->Sumw2();

  TH1D *hsvtxdl = new TH1D("hsvtxdl","hsvtxdl",20,0,10);
  TH1D *hsvtxdlB = new TH1D("hsvtxdlB","hsvtxdlB",20,0,10);
  TH1D *hsvtxdlC = new TH1D("hsvtxdlC","hsvtxdlC",20,0,10);
  TH1D *hsvtxdlL = new TH1D("hsvtxdlL","hsvtxdlL",20,0,10);
  hsvtxdl->Sumw2(); hsvtxdlB->Sumw2(); hsvtxdlC->Sumw2(); hsvtxdlL->Sumw2();

  TH1D *hsvtxdls = new TH1D("hsvtxdls","hsvtxdls",40,0,80);
  TH1D *hsvtxdlsB = new TH1D("hsvtxdlsB","hsvtxdlsB",40,0,80);
  TH1D *hsvtxdlsC = new TH1D("hsvtxdlsC","hsvtxdlsC",40,0,80);
  TH1D *hsvtxdlsL = new TH1D("hsvtxdlsL","hsvtxdlsL",40,0,80);
  hsvtxdls->Sumw2(); hsvtxdlsB->Sumw2(); hsvtxdlsC->Sumw2(); hsvtxdlsL->Sumw2();
  
  TH1D *hsvtxm = new TH1D("hsvtxm","hsvtxm",32,0,8);
  TH1D *hsvtxmB = new TH1D("hsvtxmB","hsvtxmB",32,0,8);
  TH1D *hsvtxmC = new TH1D("hsvtxmC","hsvtxmC",32,0,8);
  TH1D *hsvtxmL = new TH1D("hsvtxmL","hsvtxmL",32,0,8);
  hsvtxm->Sumw2(); hsvtxmB->Sumw2(); hsvtxmC->Sumw2(); hsvtxmL->Sumw2();
  
  TH1D *hsvtxmSV3 = new TH1D("hsvtxmSV3","hsvtxmSV3",32,0,8);
  TH1D *hsvtxmSV3B = new TH1D("hsvtxmSV3B","hsvtxmSV3B",32,0,8);
  TH1D *hsvtxmSV3C = new TH1D("hsvtxmSV3C","hsvtxmSV3C",32,0,8);
  TH1D *hsvtxmSV3L = new TH1D("hsvtxmSV3L","hsvtxmSV3L",32,0,8);
  hsvtxmSV3->Sumw2(); hsvtxmSV3B->Sumw2(); hsvtxmSV3C->Sumw2(); hsvtxmSV3L->Sumw2();
  
  TH1D *hsvtxpt = new TH1D("hsvtxpt","hsvtxpt",20,0,100);
  TH1D *hsvtxptB = new TH1D("hsvtxptB","hsvtxptB",20,0,100);
  TH1D *hsvtxptC = new TH1D("hsvtxptC","hsvtxptC",20,0,100);
  TH1D *hsvtxptL = new TH1D("hsvtxptL","hsvtxptL",20,0,100);
  hsvtxpt->Sumw2(); hsvtxptB->Sumw2(); hsvtxptC->Sumw2(); hsvtxptL->Sumw2();
  
  TH1D *hsvtxptSV3 = new TH1D("hsvtxptSV3","hsvtxptSV3",20,0,100);
  TH1D *hsvtxptSV3B = new TH1D("hsvtxptSV3B","hsvtxptSV3B",20,0,100);
  TH1D *hsvtxptSV3C = new TH1D("hsvtxptSV3C","hsvtxptSV3C",20,0,100);
  TH1D *hsvtxptSV3L = new TH1D("hsvtxptSV3L","hsvtxptSV3L",20,0,100);
  hsvtxptSV3->Sumw2(); hsvtxptSV3B->Sumw2(); hsvtxptSV3C->Sumw2(); hsvtxptSV3L->Sumw2();
  
  TH1D *hnIPtrk = new TH1D("hnIPtrk","hnIPtrk",100,0,100);
  TH1D *hnIPtrkB = new TH1D("hnIPtrkB","hnIPtrkB",100,0,100);
  TH1D *hnIPtrkC = new TH1D("hnIPtrkC","hnIPtrkC",100,0,100);
  TH1D *hnIPtrkL = new TH1D("hnIPtrkL","hnIPtrkL",100,0,100);
  hnIPtrk->Sumw2(); hnIPtrkB->Sumw2(); hnIPtrkC->Sumw2(); hnIPtrkL->Sumw2();
  
  TH1D *hnselIPtrk = new TH1D("hnselIPtrk","hnselIPtrk",100,0,100);
  TH1D *hnselIPtrkB = new TH1D("hnselIPtrkB","hnselIPtrkB",100,0,100);
  TH1D *hnselIPtrkC = new TH1D("hnselIPtrkC","hnselIPtrkC",100,0,100);
  TH1D *hnselIPtrkL = new TH1D("hnselIPtrkL","hnselIPtrkL",100,0,100);
  hnselIPtrk->Sumw2(); hnselIPtrkB->Sumw2(); hnselIPtrkC->Sumw2(); hnselIPtrkL->Sumw2();
  
  TH1D *hmuptrel = new TH1D("hmuptrel","hmuptrel",40,0,4);
  TH1D *hmuptrelB = new TH1D("hmuptrelB","hmuptrelB",40,0,4);
  TH1D *hmuptrelC = new TH1D("hmuptrelC","hmuptrelC",40,0,4);
  TH1D *hmuptrelL = new TH1D("hmuptrelL","hmuptrelL",40,0,4);
  hmuptrel->Sumw2(); hmuptrelB->Sumw2(); hmuptrelC->Sumw2(); hmuptrelL->Sumw2();
  
  TH1D *hmuptrelSV2 = new TH1D("hmuptrelSV2","hmuptrelSV2",40,0,4);
  TH1D *hmuptrelSV2B = new TH1D("hmuptrelSV2B","hmuptrelSV2B",40,0,4);
  TH1D *hmuptrelSV2C = new TH1D("hmuptrelSV2C","hmuptrelSV2C",40,0,4);
  TH1D *hmuptrelSV2L = new TH1D("hmuptrelSV2L","hmuptrelSV2L",40,0,4);
  hmuptrelSV2->Sumw2(); hmuptrelSV2B->Sumw2(); hmuptrelSV2C->Sumw2(); hmuptrelSV2L->Sumw2();
  
  TH1D *hmuptrelSV3 = new TH1D("hmuptrelSV3","hmuptrelSV3",40,0,4);
  TH1D *hmuptrelSV3B = new TH1D("hmuptrelSV3B","hmuptrelSV3B",40,0,4);
  TH1D *hmuptrelSV3C = new TH1D("hmuptrelSV3C","hmuptrelSV3C",40,0,4);
  TH1D *hmuptrelSV3L = new TH1D("hmuptrelSV3L","hmuptrelSV3L",40,0,4);
  hmuptrelSV3->Sumw2(); hmuptrelSV3B->Sumw2(); hmuptrelSV3C->Sumw2(); hmuptrelSV3L->Sumw2();
  
  TH1D *hipPt = new TH1D("hipPt","hipPt",40,0,40);
  TH1D *hipPtB = new TH1D("hipPtB","hipPtB",40,0,40);
  TH1D *hipPtC = new TH1D("hipPtC","hipPtC",40,0,40);
  TH1D *hipPtL = new TH1D("hipPtL","hipPtL",40,0,40);
  hipPt->Sumw2(); hipPtB->Sumw2(); hipPtC->Sumw2(); hipPtL->Sumw2();
  
  TH1D *hipProb0 = new TH1D("hipProb0","hipProb0",40,-1,1);
  TH1D *hipProb0B = new TH1D("hipProb0B","hipProb0B",40,-1,1);
  TH1D *hipProb0C = new TH1D("hipProb0C","hipProb0C",40,-1,1);
  TH1D *hipProb0L = new TH1D("hipProb0L","hipProb0L",40,-1,1);
  hipProb0->Sumw2(); hipProb0B->Sumw2(); hipProb0C->Sumw2(); hipProb0L->Sumw2();
  
  TH1D *hipProb1 = new TH1D("hipProb1","hipProb1",40,-1,1);
  TH1D *hipProb1B = new TH1D("hipProb1B","hipProb1B",40,-1,1);
  TH1D *hipProb1C = new TH1D("hipProb1C","hipProb1C",40,-1,1);
  TH1D *hipProb1L = new TH1D("hipProb1L","hipProb1L",40,-1,1);
  hipProb1->Sumw2(); hipProb1B->Sumw2(); hipProb1C->Sumw2(); hipProb1L->Sumw2();
  
  TH1D *hip2d = new TH1D("hip2d","hip2d",40,-0.1,0.1);
  TH1D *hip2dB = new TH1D("hip2dB","hip2dB",40,-0.1,0.1);
  TH1D *hip2dC = new TH1D("hip2dC","hip2dC",40,-0.1,0.1);
  TH1D *hip2dL = new TH1D("hip2dL","hip2dL",40,-0.1,0.1);
  hip2d->Sumw2(); hip2dB->Sumw2(); hip2dC->Sumw2(); hip2dL->Sumw2();
  
  TH1D *hip2dSig = new TH1D("hip2dSig","hip2dSig",70,-35,35);
  TH1D *hip2dSigB = new TH1D("hip2dSigB","hip2dSigB",70,-35,35);
  TH1D *hip2dSigC = new TH1D("hip2dSigC","hip2dSigC",70,-35,35);
  TH1D *hip2dSigL = new TH1D("hip2dSigL","hip2dSigL",70,-35,35);
  hip2dSig->Sumw2(); hip2dSigB->Sumw2(); hip2dSigC->Sumw2(); hip2dSigL->Sumw2();

  TH1D *hip2d1 = new TH1D("hip2d1","hip2d1",40,-0.1,0.1);
  TH1D *hip2d1B = new TH1D("hip2d1B","hip2d1B",40,-0.1,0.1);
  TH1D *hip2d1C = new TH1D("hip2d1C","hip2d1C",40,-0.1,0.1);
  TH1D *hip2d1L = new TH1D("hip2d1L","hip2d1L",40,-0.1,0.1);
  hip2d1->Sumw2(); hip2d1B->Sumw2(); hip2d1C->Sumw2(); hip2d1L->Sumw2();

  TH1D *hip2dSig1 = new TH1D("hip2dSig1","hip2dSig1",70,-35,35);
  TH1D *hip2dSig1B = new TH1D("hip2dSig1B","hip2dSig1B",70,-35,35);
  TH1D *hip2dSig1C = new TH1D("hip2dSig1C","hip2dSig1C",70,-35,35);
  TH1D *hip2dSig1L = new TH1D("hip2dSig1L","hip2dSig1L",70,-35,35);
  hip2dSig1->Sumw2(); hip2dSig1B->Sumw2(); hip2dSig1C->Sumw2(); hip2dSig1L->Sumw2();

  TH1D *hip2d2 = new TH1D("hip2d2","hip2d2",40,-0.1,0.1);
  TH1D *hip2d2B = new TH1D("hip2d2B","hip2d2B",40,-0.1,0.1);
  TH1D *hip2d2C = new TH1D("hip2d2C","hip2d2C",40,-0.1,0.1);
  TH1D *hip2d2L = new TH1D("hip2d2L","hip2d2L",40,-0.1,0.1);
  hip2d2->Sumw2(); hip2d2B->Sumw2(); hip2d2C->Sumw2(); hip2d2L->Sumw2();

  TH1D *hip2dSig2 = new TH1D("hip2dSig2","hip2dSig2",70,-35,35);
  TH1D *hip2dSig2B = new TH1D("hip2dSig2B","hip2dSig2B",70,-35,35);
  TH1D *hip2dSig2C = new TH1D("hip2dSig2C","hip2dSig2C",70,-35,35);
  TH1D *hip2dSig2L = new TH1D("hip2dSig2L","hip2dSig2L",70,-35,35);
  hip2dSig2->Sumw2(); hip2dSig2B->Sumw2(); hip2dSig2C->Sumw2(); hip2dSig2L->Sumw2();

  TH1D *hip2d3 = new TH1D("hip2d3","hip2d3",40,-0.1,0.1);
  TH1D *hip2d3B = new TH1D("hip2d3B","hip2d3B",40,-0.1,0.1);
  TH1D *hip2d3C = new TH1D("hip2d3C","hip2d3C",40,-0.1,0.1);
  TH1D *hip2d3L = new TH1D("hip2d3L","hip2d3L",40,-0.1,0.1);
  hip2d3->Sumw2(); hip2d3B->Sumw2(); hip2d3C->Sumw2(); hip2d3L->Sumw2();

  TH1D *hip2dSig3 = new TH1D("hip2dSig3","hip2dSig3",70,-35,35);
  TH1D *hip2dSig3B = new TH1D("hip2dSig3B","hip2dSig3B",70,-35,35);
  TH1D *hip2dSig3C = new TH1D("hip2dSig3C","hip2dSig3C",70,-35,35);
  TH1D *hip2dSig3L = new TH1D("hip2dSig3L","hip2dSig3L",70,-35,35);
  hip2dSig3->Sumw2(); hip2dSig3B->Sumw2(); hip2dSig3C->Sumw2(); hip2dSig3L->Sumw2();
  
  TH1D *hip3d = new TH1D("hip3d","hip3d",40,-0.1,0.1);
  TH1D *hip3dB = new TH1D("hip3dB","hip3dB",40,-0.1,0.1);
  TH1D *hip3dC = new TH1D("hip3dC","hip3dC",40,-0.1,0.1);
  TH1D *hip3dL = new TH1D("hip3dL","hip3dL",40,-0.1,0.1);
  hip3d->Sumw2(); hip3dB->Sumw2(); hip3dC->Sumw2(); hip3dL->Sumw2();
  
  TH1D *hip3dSig = new TH1D("hip3dSig","hip3dSig",70,-35,35);
  TH1D *hip3dSigB = new TH1D("hip3dSigB","hip3dSigB",70,-35,35);
  TH1D *hip3dSigC = new TH1D("hip3dSigC","hip3dSigC",70,-35,35);
  TH1D *hip3dSigL = new TH1D("hip3dSigL","hip3dSigL",70,-35,35);
  hip3dSig->Sumw2(); hip3dSigB->Sumw2(); hip3dSigC->Sumw2(); hip3dSigL->Sumw2();
  
  TH1D *hip3d1 = new TH1D("hip3d1","hip3d1",40,-0.1,0.1);
  TH1D *hip3d1B = new TH1D("hip3d1B","hip3d1B",40,-0.1,0.1);
  TH1D *hip3d1C = new TH1D("hip3d1C","hip3d1C",40,-0.1,0.1);
  TH1D *hip3d1L = new TH1D("hip3d1L","hip3d1L",40,-0.1,0.1);
  hip3d1->Sumw2(); hip3d1B->Sumw2(); hip3d1C->Sumw2(); hip3d1L->Sumw2();

  TH1D *hip3dSig1 = new TH1D("hip3dSig1","hip3dSig1",70,-35,35);
  TH1D *hip3dSig1B = new TH1D("hip3dSig1B","hip3dSig1B",70,-35,35);
  TH1D *hip3dSig1C = new TH1D("hip3dSig1C","hip3dSig1C",70,-35,35);
  TH1D *hip3dSig1L = new TH1D("hip3dSig1L","hip3dSig1L",70,-35,35);
  hip3dSig1->Sumw2(); hip3dSig1B->Sumw2(); hip3dSig1C->Sumw2(); hip3dSig1L->Sumw2();

  TH1D *hip3d2 = new TH1D("hip3d2","hip3d2",40,-0.1,0.1);
  TH1D *hip3d2B = new TH1D("hip3d2B","hip3d2B",40,-0.1,0.1);
  TH1D *hip3d2C = new TH1D("hip3d2C","hip3d2C",40,-0.1,0.1);
  TH1D *hip3d2L = new TH1D("hip3d2L","hip3d2L",40,-0.1,0.1);
  hip3d2->Sumw2(); hip3d2B->Sumw2(); hip3d2C->Sumw2(); hip3d2L->Sumw2();

  TH1D *hip3dSig2 = new TH1D("hip3dSig2","hip3dSig2",70,-35,35);
  TH1D *hip3dSig2B = new TH1D("hip3dSig2B","hip3dSig2B",70,-35,35);
  TH1D *hip3dSig2C = new TH1D("hip3dSig2C","hip3dSig2C",70,-35,35);
  TH1D *hip3dSig2L = new TH1D("hip3dSig2L","hip3dSig2L",70,-35,35);
  hip3dSig2->Sumw2(); hip3dSig2B->Sumw2(); hip3dSig2C->Sumw2(); hip3dSig2L->Sumw2();

  TH1D *hip3d3 = new TH1D("hip3d3","hip3d3",40,-0.1,0.1);
  TH1D *hip3d3B = new TH1D("hip3d3B","hip3d3B",40,-0.1,0.1);
  TH1D *hip3d3C = new TH1D("hip3d3C","hip3d3C",40,-0.1,0.1);
  TH1D *hip3d3L = new TH1D("hip3d3L","hip3d3L",40,-0.1,0.1);
  hip3d3->Sumw2(); hip3d3B->Sumw2(); hip3d3C->Sumw2(); hip3d3L->Sumw2();

  TH1D *hip3dSig3 = new TH1D("hip3dSig3","hip3dSig3",70,-35,35);
  TH1D *hip3dSig3B = new TH1D("hip3dSig3B","hip3dSig3B",70,-35,35);
  TH1D *hip3dSig3C = new TH1D("hip3dSig3C","hip3dSig3C",70,-35,35);
  TH1D *hip3dSig3L = new TH1D("hip3dSig3L","hip3dSig3L",70,-35,35);
  hip3dSig3->Sumw2(); hip3dSig3B->Sumw2(); hip3dSig3C->Sumw2(); hip3dSig3L->Sumw2();

  TH1D *hipDist2Jet = new TH1D("hipDist2Jet","hipDist2Jet",40,-0.1,0);
  TH1D *hipDist2JetB = new TH1D("hipDist2JetB","hipDist2JetB",40,-0.1,0);
  TH1D *hipDist2JetC = new TH1D("hipDist2JetC","hipDist2JetC",40,-0.1,0);
  TH1D *hipDist2JetL = new TH1D("hipDist2JetL","hipDist2JetL",40,-0.1,0);
  hipDist2Jet->Sumw2(); hipDist2JetB->Sumw2(); hipDist2JetC->Sumw2(); hipDist2JetL->Sumw2();
  
  TH1D *hipDist2JetSig = new TH1D("hipDist2JetSig","hipDist2JetSig",40,-0.1,0.1);
  TH1D *hipDist2JetSigB = new TH1D("hipDist2JetSigB","hipDist2JetSigB",40,-0.1,0.1);
  TH1D *hipDist2JetSigC = new TH1D("hipDist2JetSigC","hipDist2JetSigC",40,-0.1,0.1);
  TH1D *hipDist2JetSigL = new TH1D("hipDist2JetSigL","hipDist2JetSigL",40,-0.1,0.1);
  hipDist2JetSig->Sumw2(); hipDist2JetSigB->Sumw2(); hipDist2JetSigC->Sumw2(); hipDist2JetSigL->Sumw2();
  
  TH1D *hipClosest2Jet = new TH1D("hipClosest2Jet","hipClosest2Jet",40,0,1);
  TH1D *hipClosest2JetB = new TH1D("hipClosest2JetB","hipClosest2JetB",40,0,1);
  TH1D *hipClosest2JetC = new TH1D("hipClosest2JetC","hipClosest2JetC",40,0,1);
  TH1D *hipClosest2JetL = new TH1D("hipClosest2JetL","hipClosest2JetL",40,0,1);
  hipClosest2Jet->Sumw2(); hipClosest2JetB->Sumw2(); hipClosest2JetC->Sumw2(); hipClosest2JetL->Sumw2();

  Double_t t_jtpt, t_jteta, t_jtphi, t_rawpt, t_refpt, t_subid, t_discr_prob, t_discr_ssvHighEff, t_discr_ssvHighPur, t_discr_csvSimple, t_svtxm;
  Double_t t_pthat, t_weight;
  Int_t t_refparton_flavorForB;
  Int_t trigIndex, t_bin;
  Double_t         t_hiHFplusEta4;
  Double_t         t_hiHFminusEta4;
  Int_t t_nIP;
  Double_t t_ipPt[100], t_ipProb0[100];
  Int_t t_ipJetIndex[100];

  TTree *nt = new TTree("nt","");
  nt->Branch("jtpt",&t_jtpt,"jtpt/D");
  nt->Branch("jteta",&t_jteta,"jteta/D");
  nt->Branch("jtphi",&t_jtphi,"jtphi/D");
  nt->Branch("subid",&t_subid,"subid/I");
  nt->Branch("rawpt",&t_rawpt,"rawpt/D");
  nt->Branch("refpt",&t_refpt,"refpt/D");
  nt->Branch("refparton_flavorForB",&t_refparton_flavorForB,"refparton_flavorForB/I");
  nt->Branch("discr_prob",&t_discr_prob,"discr_prob/D");
  nt->Branch("discr_ssvHighEff",&t_discr_ssvHighEff,"discr_ssvHighEff/D");
  nt->Branch("discr_ssvHighPur",&t_discr_ssvHighPur,"discr_ssvHighPur/D");
  nt->Branch("discr_csvSimple",&t_discr_csvSimple,"discr_csvSimple/D");
  nt->Branch("svtxm",&t_svtxm,"svtxm/D");
  nt->Branch("bin",&t_bin,"bin/I");
  nt->Branch("hiHFplusEta4",&t_hiHFplusEta4,"hiHFplusEta4/D");
  nt->Branch("hiHFminusEta4",&t_hiHFminusEta4,"hiHFminusEta4/D");
  if(ppPbPb){
    nt->Branch("trigIndex",&trigIndex,"trigIndex/I");
  }
  if(ExpandedTree){
    nt->Branch("nIP",&t_nIP);
    nt->Branch("ipPt",t_ipPt,"ipPt[nIP]/D");
    nt->Branch("ipProb0",t_ipProb0,"ipProb0[nIP]/D");
    nt->Branch("ipJetIndex",t_ipJetIndex,"ipJetIndex[nIP]/I");
  }
  if(!ppPbPb){
      nt->Branch("HLT_PAZeroBiasPixel_SingleTrack_v1 ",&HLT_PAZeroBiasPixel_SingleTrack_v1,"HLT_PAZeroBiasPixel_SingleTrack_v1/I");
    nt->Branch("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_NoJetID_v1,"HLT_PAJet20_noJetID_v1/I");
    nt->Branch("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_NoJetID_v1,"HLT_PAJet40_noJetID_v1/I");
    nt->Branch("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_NoJetID_v1,"HLT_PAJet60_noJetID_v1/I");
    nt->Branch("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_NoJetID_v1,"HLT_PAJet80_noJetID_v1/I");
    nt->Branch("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_NoJetID_v1,"HLT_PAJet100_noJetID_v1/I");
    nt->Branch("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP,"pVertexFilterCutGplusUpsPP/I");
  }

  if(isMC) nt->Branch("pthat",&t_pthat,"pthat/D");
  nt->Branch("weight",&t_weight,"weight/D");

  // grab the JEC's
   
  //JetCorrectorParameters* parHI442x_l2, * parHI442x_l3;
  // vector<JetCorrectorParameters> vpar_HI442x;
  // FactorizedJetCorrector *_JEC_HI442x=NULL;
   
  if(updateJEC){
     
    //cout<<" updating the JECs, USING REGPF "<<endl;

    //string L2Name = "JEC/JEC_regPF_L2Relative_AK3PF.txt";
    //string L3Name = "JEC/JEC_regPF_L3Absolute_AK3PF.txt";
    string L2Name = "JEC/JEC_dijet_L2Relative_AK3PF.txt";
    string L3Name = "JEC/JEC_dijet_L3Absolute_AK3PF.txt";
     
    // parHI442x_l2 = new JetCorrectorParameters(L2Name.c_str());
    //parHI442x_l3 = new JetCorrectorParameters(L3Name.c_str());
          
    // vpar_HI442x.push_back(*parHI442x_l2);
    // vpar_HI442x.push_back(*parHI442x_l3);
    // _JEC_HI442x = new FactorizedJetCorrector(vpar_HI442x);
  }
   
  std::ifstream instr(infile.c_str(), std::ifstream::in);
  std::string filename;
  int nFiles=0;
  if(ppPbPb) nFiles=1;
  else if(isMC){
    nFiles=QCDpthatBins;
  }
  else{
    nFiles=dataFiles;
  }
  for(int ifile=0; ifile<nFiles; ifile++){
    
    //Add b/c statistics to the HF statistics
    if(!ppPbPb){
      if((isMC && ifile<QCDpthatBins) || !isMC){
        instr >> filename;
      }
      std::cout << "File: " << filename << std::endl;
      fin = TFile::Open(filename.c_str());
    }
    TTree *t;
  //  if(usePUsub) 
       t = (TTree*) fin->Get(Form("%sJetAnalyzer/t", algo.Data()));
  //  else t = (TTree*) fin->Get("ak3PFJetAnalyzer/t");
    TTree *tSkim = (TTree*) fin->Get("skimanalysis/HltTree");
    TTree *tEvt = NULL;
    if(!ppPbPb) tEvt = (TTree*) fin->Get("hiEvtAnalyzer/HiTree");
    TTree *tHlt = NULL;
    if(!ppPbPb) tHlt = (TTree*) fin->Get("hltanalysis/HltTree");
    if(!t || !tSkim || (!tEvt&&!ppPbPb) || (!tHlt&&!ppPbPb)){ cout << "Error! Can't find one of the trees!" << endl; exit(0);}
     
    if(tEvt) t->AddFriend("hiEvtAnalyzer/HiTree");
    if(tHlt) t->AddFriend("hltanalysis/HltTree");
    if(tSkim) t->AddFriend("skimanalysis/HltTree");
     
    Long64_t nentries = t->GetEntries();

    t->SetBranchAddress("evt",&evt);
    t->SetBranchAddress("lumi",&lumi);
    if(cbin != -1 || ppPbPb) t->SetBranchAddress("bin",&bin);
    if(!isMC) t->SetBranchAddress("run",&run);
    if(ppPbPb) t->SetBranchAddress("hf",&hf);
    t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("nref",&nref);
    t->SetBranchAddress("rawpt",rawpt);
    t->SetBranchAddress("jtpt",jtpt);
    t->SetBranchAddress("jteta",jteta);
    t->SetBranchAddress("jtphi",jtphi);
    if(!ppPbPb){
      t->SetBranchAddress("jtpu",jtpu);
      t->SetBranchAddress("jty",jty);
    }
    t->SetBranchAddress("discr_ssvHighEff",discr_ssvHighEff);
    t->SetBranchAddress("discr_ssvHighPur",discr_ssvHighPur);
    //t->SetBranchAddress("discr_csvMva",discr_csvMva);
    t->SetBranchAddress("discr_csvSimple",discr_csvSimple);
    //t->SetBranchAddress("discr_muByIp3",discr_muByIp3);
    t->SetBranchAddress("discr_muByPt",discr_muByPt);
    t->SetBranchAddress("discr_prob",discr_prob);

    t->SetBranchAddress("discr_probb",discr_probb);
    t->SetBranchAddress("discr_tcHighEff",discr_tcHighEff);
    t->SetBranchAddress("discr_tcHighPur",discr_tcHighPur);
    t->SetBranchAddress("nsvtx",nsvtx);
    t->SetBranchAddress("svtxntrk",svtxntrk);
    t->SetBranchAddress("svtxdl",svtxdl);
    t->SetBranchAddress("svtxdls",svtxdls);
    t->SetBranchAddress("svtxpt",svtxpt);
    
    t->SetBranchAddress("svtxm",svtxm);

    t->SetBranchAddress("nIPtrk",nIPtrk);
    t->SetBranchAddress("nselIPtrk",nselIPtrk);
    t->SetBranchAddress("nIP",&nIP);

    if(doTracks){
      t->SetBranchAddress("ipJetIndex",ipJetIndex);
      t->SetBranchAddress("ipPt",ipPt);
      t->SetBranchAddress("ipProb0",ipProb0);
      //t->SetBranchAddress("ipProb1",ipProb1);
      t->SetBranchAddress("ip2d",ip2d);
      t->SetBranchAddress("ip2dSig",ip2dSig);
      t->SetBranchAddress("ip3d",ip3d);
      t->SetBranchAddress("ip3dSig",ip3dSig);
      t->SetBranchAddress("ipDist2Jet",ipDist2Jet);
      //t->SetBranchAddress("ipDist2JetSig",ipDist2JetSig);
      t->SetBranchAddress("ipClosest2Jet",ipClosest2Jet);
    }

    t->SetBranchAddress("mupt",mupt);
    if(ppPbPb) t->SetBranchAddress("muptPF",muptPF);
    if(!ppPbPb) t->SetBranchAddress("pVertexFilterCutGplusUpsPP",&pVertexFilterCutGplusUpsPP);

    /*
t->SetBranchAddress("mue",mue);
t->SetBranchAddress("mueta",mueta);
t->SetBranchAddress("muphi",muphi);
t->SetBranchAddress("mudr",mudr);
t->SetBranchAddress("muptrel",muptrel);
t->SetBranchAddress("muchg",muchg);
*/
    if(isMC){
      t->SetBranchAddress("pthat",&pthat);
      t->SetBranchAddress("beamId1",&beamId1);
      t->SetBranchAddress("beamId2",&beamId2);
      t->SetBranchAddress("refpt",refpt);
      t->SetBranchAddress("subid",subid);
      t->SetBranchAddress("refeta",refeta);
      t->SetBranchAddress("refy",refy);
      t->SetBranchAddress("refphi",refphi);
      t->SetBranchAddress("refdphijt",refdphijt);
      t->SetBranchAddress("refdrjt",refdrjt);
      t->SetBranchAddress("refparton_pt",refparton_pt);
      t->SetBranchAddress("refparton_flavor",refparton_flavor);
      t->SetBranchAddress("refparton_flavorForB",refparton_flavorForB);
      t->SetBranchAddress("refparton_isGSP",refparton_isGSP);
      

      TBranch* tweight;
    //  if(isMC){
        tweight = t->GetBranch("weight");
        if(!tweight){
         if(ifile==0){
         cout << "Weight not found in Tree. Calculating..." << endl;
         useWeight=0;
         }
        }
        if(!ppPbPb && !useWeight && ifile==0){
         MCentr = countMCevents(infile, usePUsub, isMC);
        // MCentr = (int) nentries;
         // if(isMC>1){
         // for(int lm=HFpthatBins+2; lm<QCDpthatBins+1; lm++){
         // MCentr[HFpthatBins] += MCentr[lm]; //hack because we go to pthat bin 540 in QCD jet and only pthat bin 170 in b/c jet MC
         // }
         // }
         for(int i=0; i<QCDpthatBins; i++){
         cout << "MCentr["<<i<<"]: " << *(MCentr+i) << endl;
         }
        }
  //    }
    }
    t->SetBranchAddress("chargedMax",chargedMax);
    t->SetBranchAddress("photonMax",photonMax);
    t->SetBranchAddress("neutralMax",neutralMax);
    t->SetBranchAddress("chargedSum",chargedSum);
    t->SetBranchAddress("photonSum",photonSum);
    t->SetBranchAddress("neutralSum",neutralSum);
    t->SetBranchAddress("muSum",muSum);
    t->SetBranchAddress("eSum",eSum);
    
    if(isMC&&useWeight){
      t->SetBranchAddress("weight",&weight);
      t->SetBranchAddress("xSecWeight",&xSecWeight);
      if(ppPbPb)t->SetBranchAddress("centWeight",&centWeight);
      t->SetBranchAddress("vzWeight",&vzWeight);
    }

    if(ppPbPb){
      t->SetBranchAddress("nHLTBit",&nHLTBit);
      t->SetBranchAddress("hltBit",hltBit);
      
      tSkim->SetBranchAddress("pvSel",&pvSel);
      tSkim->SetBranchAddress("hbheNoiseSel",&hbheNoiseSel);
      tSkim->SetBranchAddress("spikeSel",&spikeSel);
      tSkim->SetBranchAddress("collSell",&collSell);
    }
    if(!ppPbPb){
      t->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
      t->SetBranchAddress("HLT_PAJet20_NoJetID_v1",&HLT_PAJet20_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&HLT_PAJet40_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&HLT_PAJet60_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_PAJet80_NoJetID_v1);
      t->SetBranchAddress("HLT_PAJet100_NoJetID_v1",&HLT_PAJet100_NoJetID_v1);
      t->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
      t->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
      t->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
      t->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4);
      t->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4);
    }


    int gspCounter=0;
    //nentries=10;
    for (Long64_t i=0; i<nentries;i++) {
       
      if (i%100000==0) cout<<" i = "<<i<<" out of "<<nentries<<" ("<<(int)(100*(float)i/(float)nentries)<<"%)"<<endl;
      
      tSkim->GetEntry(i);
      t->GetEntry(i);
      tEvt->GetEntry(i);
      if(ppPbPb && isMC){
        // temporarily remove cuts from MC
        if(!pvSel||!spikeSel) continue; //hbheNoise doesn't work in mixed events
      }
      if(ppPbPb){
        //if(!pvSel||!hbheNoiseSel||!spikeSel) continue;
        // turn off spike and on coll Sel
        if(!pvSel||!hbheNoiseSel||!collSell){
         //cout<<" selection failed, pvSel="<<pvSel<<", hbheNoiseSel="<<hbheNoiseSel<<" , collSell="<<collSell<<endl;
         continue;
        }
      }

      //Cut to remove events that correspond to the twiki "good events" but not the golden lumi filter
      if(!isMC){
        if(((int)run>211256 ) || ((int)run<210676 )) continue;
      }

      if(!ppPbPb){
        if(!isMC){
          if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA) continue;
    //  if(!HLT_PAZeroBiasPixel_SingleTrack_v1 && !HLT_PAJet40_NoJetID_v1 && !HLT_PAJet60_NoJetID_v1 && !HLT_PAJet80_NoJetID_v1) continue;
        }
      /*  else{
        //  if(!pHBHENoiseFilter || !pPAcollisionEventSelectionPA) continue;
          if(!pHBHENoiseFilter ) continue;
        }
     */  
       //! for 0-90% centrality bins selection using HF sum energy
       t_hiHFplusEta4 = hiHFplusEta4 ;
       t_hiHFminusEta4 = hiHFminusEta4 ;
     // if((t_hiHFplusEta4+t_hiHFminusEta4)<2.87) continue ;
      }

      
      if(ppPbPb){
        if(cbin==-1){
         // do nothing
         t_bin=bin;
        }
        else if(cbin==0){
         if(bin>=8) continue;
         else t_bin=bin;
        }
        else if(cbin==1) {
         if(bin<8||bin>=20) continue;
         else t_bin=bin;
        }
        else if(cbin==2){
         if(bin<20) continue;
         else t_bin=bin;
        }
        else {
         cout<<" bin not defined "<<endl;
         return;
        }
      }
      else t_bin=100;
      
//      if(isMC&&!ppPbPb){
//        if(beamId1==2112 || beamId2==2112) continue;
//      }
      

      if(fabs(vz)>15.) continue;
      
      // pileup rejection
      if(ppPbPb && hf>150000.){
        cout<<" rejecting pileup, "<<" hf "<<hf<<" bin "<<bin<<endl;
        for(int ij=0;ij<nref;ij++) if(jtpt[ij]>20.&&fabs(jteta[ij])<3.)cout<<" # associated tracks = "<<nselIPtrk[ij]<<endl;
        continue;
      }
      bool isNoise=false;
      if(ppPbPb){
        for(int ij=0; ij<nref; ij++){        
         if(jtpt[ij]>minJetPt&&fabs(jteta[ij])<3){
         if(neutralMax[ij]/(neutralMax[ij]+chargedMax[ij]+photonMax[ij])>0.975){
         cout<<" cleaning event with jet of "<<jtpt[ij]<<", eta "<<jteta[ij]<<" noise = "<<neutralMax[ij]/(neutralMax[ij]+chargedMax[ij]+photonMax[ij])<<endl;
         isNoise=true;
         }
         if(chargedSum[ij]+photonSum[ij]+neutralSum[ij]+muSum[ij]+eSum[ij]<0.5*rawpt[ij]){
         cout<<" cleaning event with jet of "<<jtpt[ij]<<", eta "<<jteta[ij]<<" sum PF pt = "<<chargedSum[ij]+photonSum[ij]+neutralSum[ij]+muSum[ij]+eSum[ij]<<endl;
         isNoise=true;
         }
         }
        }
      }
      if(isNoise) continue;
      
        
      
  /*    if(updateJEC){
        
        //for(int ij=0; ij<nref; ij++){        
        // _JEC_HI442x->setJetEta(jteta[ij]);
        // _JEC_HI442x->setJetPt(rawpt[ij]);
        // jtpt[ij] = rawpt[ij]*_JEC_HI442x->getCorrection();
        //}        
      }
    */  
      if(useWeight){
        if(isMC)w=weight;
      }
      //trigger weighting in pp data
      if(!ppPbPb && !isMC){
        bool trgDec[6] = {(bool)HLT_PAZeroBiasPixel_SingleTrack_v1, (bool)HLT_PAJet20_NoJetID_v1, (bool)HLT_PAJet40_NoJetID_v1, (bool)HLT_PAJet60_NoJetID_v1, (bool)HLT_PAJet80_NoJetID_v1, (bool)HLT_PAJet100_NoJetID_v1};
      //  w = trigComb(trgDec, pscls);
      }

      
      //Do the weighting = x-sec / Nentries, where Nentries is weighted differently for B/C jets and QCD jets
      if(isMC){
        t_pthat=pthat;
        int j=0;
        while(pthat>pthatbin[j] && j<QCDpthatBins) {j++;
      //  while(j<QCDpthatBins) {j++;
     //  cout <<"j==" << j <<endl ; 
      // for(j=0 ; j <QCDpthatBins ; j++){
         if(j==QCDpthatBins) 
           w = ((wght[j-1])/MCentr[j-1]); //wght[0] = pthat>15, MCentr[0] = pthat<15. I know it's dumb - bear with me.
         else 
          w = ((wght[j-1]-wght[j])/MCentr[j-1]);
      //  }
       }
      if(ifile<QCDpthatBins-1 && pthat>pthatbin[ifile+1]) continue ;
      }
      if(isMC){
        bool isFiltered=0;
        double vzWeight=1;
        int vzbin = (int) TMath::Ceil(vz+15.+0.4); // 0.4 is the pixel detector shift
      //  if(vzbin>0&&vzbin<=30)
       //     vzWeight = hDatavz->GetBinContent(vzbin)/hMCvz[isFiltered]->GetBinContent(vzbin);
            vzWeight = fVz->Eval(vz);
        t_weight=w*vzWeight;
        xSecWeight=w ;
      }
      else t_weight=w;

      int useEvent=0;
      
      int trackPosition =0;

      for(int ij=0;ij<nref;ij++){
      
        trackPosition+=nselIPtrk[ij];
      
//        if(useGSP==2){
//         if(refparton_isGSP[ij]==1){
//         gspCounter++;
//         if(gspCounter%2==0) continue;
//         }        
//        }
//        if(useGSP==3){
//         if(refparton_isGSP[ij]==0){
//         gspCounter++;
//         if(gspCounter%2==0) continue;
//         }        
//        }
          if(rawpt[ij]<30.&&fabs(jteta[ij])>maxJetEta) continue ; 
        //  if(isMC && subid[ij]!=0) continue ;
          if(isMC && jtpt[ij]>4*pthat) continue ;
        if(jtpt[ij]>minJetPt && fabs(jteta[ij])<maxJetEta){
         if(doNtuples){
        
         t_jtpt=jtpt[ij];
         t_jteta=jteta[ij];
         t_jtphi=jtphi[ij];
         t_rawpt=rawpt[ij];
         t_refpt=refpt[ij];
         t_subid=subid[ij];
         t_refparton_flavorForB=refparton_flavorForB[ij];
         t_discr_prob=discr_prob[ij];
         t_discr_ssvHighEff=discr_ssvHighEff[ij];
         t_discr_ssvHighPur=discr_ssvHighPur[ij];
         t_discr_csvSimple=discr_csvSimple[ij];
         t_svtxm=svtxm[ij];
        
         //Find jet tracks that correspond to the jet & apply proximity cuts
         if(ExpandedTree){
         t_nIP=nselIPtrk[ij];
         int counter=0;
         for(int itrk=0; itrk<nIP; itrk++){
                if(ipJetIndex[itrk] == ij){
                 t_ipProb0[counter] = ipProb0[itrk];
                 t_ipPt[counter] = ipPt[itrk];
                 t_ipJetIndex[counter] = ij;
                 counter++;
                }
         }
         }
         nt->Fill();


         }

        
         if(!doJets) continue;
        
        
         useEvent=1;
        
       //  hjtpt->Fill(jtpt[ij],w);
        if(TMath::Abs((jteta[ij]+0.465))<=1.)hjtpt->Fill(jtpt[ij],w);
          Int_t dEtaBin = -1 ;
          for(Int_t ieta = 0 ; ieta <netabin; ieta++){
           if((jteta[ij]+0.465)>deta[ieta]&&(jteta[ij]+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
          }
         if(dEtaBin!=-1) hJetPtEtaBin[dEtaBin]->Fill(jtpt[ij],w);

         hrawpt->Fill(rawpt[ij],w);

         hjteta->Fill(jteta[ij],w);

         hjtphi->Fill(jtphi[ij],w);

         //*
         hdiscr_csvSimple->Fill(discr_csvSimple[ij],w);

        
         hdiscr_prob->Fill(discr_prob[ij],w);

         hdiscr_ssvHighEff->Fill(discr_ssvHighEff[ij],w);
        
         hdiscr_ssvHighPur->Fill(discr_ssvHighPur[ij],w);

         hdiscr_tcHighEff->Fill(discr_tcHighEff[ij],w);
        
         hdiscr_tcHighPur->Fill(discr_tcHighPur[ij],w);

         //*
         hnsvtx->Fill(nsvtx[ij],w);
        
            if(nsvtx[ij]>0){
                
                hsvtxntrk->Fill(svtxntrk[ij],w);
                
                
                // require at least 1 tracks as in btagging @ 7 TeV note
                if(svtxntrk[ij]>1){
                    
                    hsvtxdl->Fill(svtxdl[ij],w);
                    
                    hsvtxdls->Fill(svtxdls[ij],w);
                    
                    
                    hsvtxm->Fill(svtxm[ij],w);
                    
                    
                    hsvtxpt->Fill(svtxpt[ij],w);
                    
                    
                    if(svtxntrk[ij]>=3) {
                        
                        hsvtxmSV3->Fill(svtxm[ij],w);
                        
                        
                        hsvtxptSV3->Fill(svtxpt[ij],w);
                        
                    }
                }
            }
        
         hnIPtrk->Fill(nIPtrk[ij],w);

        
        
         hnselIPtrk->Fill(nselIPtrk[ij],w);
        

        }  //! jet pt cuts
      }  //! end of jet loop
      if(useEvent){
        if(isMC){
         hbinw->Fill(bin,w);
         hbin->Fill(bin,xSecWeight*vzWeight);
        }
        else hbin->Fill(bin);
        
        if(isMC){
         hvzw->Fill(vz,w);
         if(ppPbPb)hvz->Fill(vz,xSecWeight*centWeight);
         else hvz->Fill(vz,xSecWeight);
        }
        else hvz->Fill(vz);
      }
      
    }  //! events loop
  } //! file loop
  fout->cd();

  hbin->Write(); hbinw->Write(); hvz->Write(); hvzw->Write();
  
  hjtpt->Write();
  for(Int_t ieta = 0 ; ieta <netabin; ieta++){
      hJetPtEtaBin[ieta]->Write();
    } 
  hrawpt->Write();
  
  hjteta->Write();
  
  hjtphi->Write();
  
  hdiscr_csvSimple->Write();
  
  hdiscr_prob->Write();
  
  hdiscr_ssvHighEff->Write();

  hdiscr_ssvHighPur->Write();

  hdiscr_tcHighEff->Write();

  hdiscr_tcHighPur->Write();

  hnsvtx->Write();

  hsvtxntrk->Write();
  
  hsvtxdl->Write();

  hsvtxdls->Write();

  hsvtxm->Write();

  hsvtxmSV3->Write();

  hsvtxpt->Write();

  hsvtxptSV3->Write();

  hnIPtrk->Write();

  hnselIPtrk->Write();

  hipPt->Write();

  hipProb0->Write();

  hipProb1->Write();

  hip2d->Write();

  hip2dSig->Write();

  hip2d1->Write();

  hip2dSig1->Write();

  hip2d2->Write();

  hip2dSig2->Write();

  hip2d3->Write();

  hip2dSig3->Write();

  hip3d->Write();

  hip3dSig->Write();

  hip3d1->Write();

  hip3dSig1->Write();

  hip3d2->Write();

  hip3dSig2->Write();

  hip3d3->Write();

  hip3dSig3->Write();


  hipDist2Jet->Write();

  hipDist2JetSig->Write();

  hipClosest2Jet->Write();

  nt->Write();
  
  fout->Close();

}



