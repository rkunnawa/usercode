// test script
	
 

  // Unfolding starts: 

	
  cout << "==================================== TEST =====================================" << endl;
	
  cout << "==================================== UNFOLD ===================================" << endl;
	
  char chmet[100];
	
  // ======================= Reconstructed pp and PbPb spectra =========================================================
  TCanvas * cPbPb = new TCanvas("cPbPb","PbPb",1200,800);
  cPbPb->Divide(2,1);
	
  for (int i=0;i<=nbins_cent;i++) {
    cPbPb->cd(i+1)->SetLogy();
	  
    // Do unfolding
    prior myPrior(uhist[i]->hMatrix,uhist[i]->hMeas,0.0);
    myPrior.unfold(uhist[i]->hMeas,1);
    TH1F *hPrior = (TH1F*)uhist[i]->hMatrix->ProjectionX()->Clone(Form("hPrior_cent%d",i));
    hPrior->Scale(uhist[i]->hMeas->Integral(0,1000)/hPrior->Integral(0,1000));
    //bayesianUnfold myUnfoldingJECSys(uhist[i]->hMatrix,hPrior,0);
    //myUnfoldingJECSys.unfold(uhist[i]->hMeasJECSys,nBayesianIter);
    //bayesianUnfold myUnfoldingSmearSys(uhist[i]->hMatrix,hPrior,0);
    //myUnfoldingSmearSys.unfold(uhist[i]->hMeasSmearSys,nBayesianIter);
    cout <<"Unfolding bin "<<i<<endl;
    // Iteration Systematics
    for (int j=2;j<7;j++){
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
    TH1F* hBinByBinCor = (TH1F*)hBinByBinCorRaw->Clone("hBinByBinCor");//functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    //      TH1F* hBinByBinCor = (TH1F*)functionHist(f,hBinByBinCorRaw,Form("hBinByBinCor_cent%d",i));
    delete hBinByBinCorRaw,hMCGen;
    uhist[i]->hRecoBinByBin = (TH1F*) uhist[i]->hMeas->Clone(Form("hRecoBinByBin_cent%d",i));
    uhist[i]->hRecoBinByBin->Divide(hBinByBinCor);
    //      uhist[i]->hRecoBinByBin = (TH1F*) hMCReco->Clone(Form("hRecoBinByBin_cent%d",i));
		
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
  hRecoBinByBinRAA->Scale(1./362.24/1147500000);
  hRecoRAA->Scale(1./362.24/1147500000);
  hMeasRAA->Scale(1./362.24/1147500000); // ncoll * total min bias events. 

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
  
  TFile f("RAA_PbPb2011_PP2013_new_vz_HLT_80_dup_removed.root","RECREATE");
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

  f.Close();
