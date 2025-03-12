#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

#include "TRandom3.h"

#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TSystemFile.h"

#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

void GetFiles(char const *input, vector<string> &files,
              string filenamePattern=".root") {
  TSystemDirectory dir(input, input);
  TList *list = dir.GetListOfFiles();

  if (list) {
    TSystemFile *file;
    string fname;
    std::cout << "Getting files from " << input << std::endl;
    TIter next(list);
    while ((file = (TSystemFile *)next())) {
      fname = file->GetName();
      std::cout << "File: " << fname << std::endl;
      if (file->IsDirectory() && (fname.find(".") == string::npos)) {
        string newDir = string(input) + fname + "/";
        GetFiles(newDir.c_str(), files);
      } else if ((fname.find(filenamePattern.c_str()) != string::npos)) {
        // continue if the filename is .root_hist
        if (fname.find("_hist") != string::npos) continue;
        files.push_back(string(input) + fname);
        cout << files.back() << endl;
      }
    }
  }

  return;
}

void FillChain(TChain &chain, vector<string> &files) {
  for (auto file : files) {
    chain.Add(file.c_str());
  }
}

// convert histogram to different format such that the range of values goes from 0 to nBins and each bin has entries of the old histogram with the same number
TH1D* convertHistToUnfoldingFormat(TH1D* inHist){
    int nBins = inHist->GetNbinsX();
    double xMin = inHist->GetXaxis()->GetXmin();
    double xMax = inHist->GetXaxis()->GetXmax();
    TH1D* outHist = new TH1D(Form("%s_unfolded", inHist->GetName()), Form("%s_unfolded", inHist->GetTitle()), nBins, 0, nBins);
    for(int iBin = 1; iBin <= nBins; iBin++){
        double binContent = inHist->GetBinContent(iBin);
        outHist->SetBinContent(iBin, binContent);
    }
    return outHist;
}

void createUnfoldingHistograms(char const *infnameMCDir){
    // pick up the matched samples
    vector<string> files;
    GetFiles(infnameMCDir, files, "_Matched.root");

    TChain pairChain("PairTree");
    FillChain(pairChain, files);
    TTreeReader pairReader(&pairChain);
    TTreeReaderValue<int> nPairs(pairReader, "NPair");
    TTreeReaderArray<double> GenE(pairReader, "GenE1");
    TTreeReaderArray<double> GenE2(pairReader, "GenE2");
    TTreeReaderArray<double> GenX(pairReader, "GenX1");
    TTreeReaderArray<double> GenY(pairReader, "GenY1");
    TTreeReaderArray<double> GenZ(pairReader, "GenZ1");
    TTreeReaderArray<double> RecoE(pairReader, "RecoE1");
    TTreeReaderArray<double> RecoE2(pairReader, "RecoE2");
    TTreeReaderArray<double> RecoX(pairReader, "RecoX1");
    TTreeReaderArray<double> RecoY(pairReader, "RecoY1");
    TTreeReaderArray<double> RecoZ(pairReader, "RecoZ1");
    TTreeReaderArray<double> DistanceGen(pairReader, "DistanceGen");
    TTreeReaderArray<double> DistanceReco(pairReader, "DistanceReco");
    TTreeReaderArray<double> Distance1(pairReader, "Distance1");
    TTreeReaderArray<double> Distance2(pairReader, "Distance2");
    TTreeReaderArray<double> E1E2Gen(pairReader, "E1E2Gen");
    TTreeReaderArray<double> E1E2Reco(pairReader, "E1E2Reco");

    TChain particleChain("MatchedTree");
    FillChain(particleChain, files);
    TTreeReader particleReader(&particleChain);
    TTreeReaderValue<int> nParticles(particleReader, "NParticle");
    TTreeReaderArray<double> GenEParticle(particleReader, "GenE");
    TTreeReaderArray<double> GenXParticle(particleReader, "GenX");
    TTreeReaderArray<double> GenYParticle(particleReader, "GenY");
    TTreeReaderArray<double> GenZParticle(particleReader, "GenZ");
    TTreeReaderArray<double> RecoEParticle(particleReader, "RecoE");
    TTreeReaderArray<double> RecoXParticle(particleReader, "RecoX");
    TTreeReaderArray<double> RecoYParticle(particleReader, "RecoY");
    TTreeReaderArray<double> RecoZParticle(particleReader, "RecoZ");
    TTreeReaderArray<double> DistanceParticle(particleReader, "Distance");



    int nBins = 100;
    std::vector<double> trackPtBinsReco = {0.2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40};
    std::vector<double> trackPtBinsGen = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50};

    std::vector<double> deltaRBinsReco = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9};
    std::vector<double> deltaRBinsGen = {0.0, 0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1};

    // create energy weighting bins, from 0 to 1000 GeV in 50 GeV bins
    std::vector<double> e1e2BinsReco = {1.0,10.0, 20.0, 50, 100,  200, 300, 500};
    std::vector<double> e1e2BinsGen = {0.1,1.0,10.0, 20.0, 50, 100, 200, 300, 500, 700};


    TH2D *hRecoGenDeltaR = new TH2D("hRecoGenDeltaR", "Reco vs Gen Delta R", nBins, 0, 3.2, nBins, 0, 3.2);
    TH2D *hRecoGenTrackPt = new TH2D("hRecoGenTrackPt", "Reco vs Gen Track Pt", nBins, 01., 50, nBins, 0, 50);
    TH2D *hRecoGenTrackE = new TH2D("hRecoGenTrackE", "Reco vs Gen Track E", nBins, 0.1, 50, nBins, 0, 50);
    TH2D *hRecoGenEnergyWeighting = new TH2D("hRecoGenEnergyWeighting", "Reco vs Gen Energy Weighting", nBins, 0.1, 1000, nBins, 0.1, 1000);
    TH1D *hDeltaDeltaR = new TH1D("hDeltaDeltaR", "Delta Delta R", nBins, -3.2, 3.2);
    TH1D *hDelta1 = new TH1D("hDelta1", "Delta 1", nBins, -3.2, 3.2);
    TH1D *hDelta2 = new TH1D("hDelta2", "Delta 2", nBins, -3.2, 3.2);
    TH1D *hDelta3 = new TH1D("hDelta3", "Delta 3", nBins, -3.2, 3.2);
  

    TH2D* singleTrackGenERecoE = new TH2D("singleTrackGenERecoE", "Single Track Gen E vs Reco E", 100, 0, 50, 100, 0, 50);
    TH1D* singleTrackMatchingDistance = new TH1D("singleTrackMatchingDistance", "Single Track Matching Distance", 100, 0, 3.2);
    TH2D* singleTrackMatchingDistanceVsEta = new TH2D("singleTrackMatchingDistanceVsEta", "Single Track Matching Distance vs Eta", 100, 0, 3.2, 100, -3.2, 3.2);

    // track pT unfolding hists
    TH1D *hTrackPtSmeared = new TH1D("hTrackPtSmeared", "Track Pt Smeared", trackPtBinsReco.size()-1, trackPtBinsReco.data());
    TH1D *hTrackPtGen = new TH1D("hTrackPtGen", "Track Pt Gen", trackPtBinsGen.size()-1, trackPtBinsGen.data());
    TH2D *hTrackPtResp = new TH2D("hTrackPtResp", "Track Pt Response", trackPtBinsReco.size()-1, trackPtBinsReco.data(), trackPtBinsGen.size()-1, trackPtBinsGen.data());
    TH1D *hTrackPtSmeared_SplitMC = new TH1D("hTrackPtSmeared_SplitMC", "Track Pt Smeared", trackPtBinsReco.size()-1, trackPtBinsReco.data());
    TH1D *hTrackPtGen_SplitMC = new TH1D("hTrackPtGen_SplitMC", "Track Pt Gen", trackPtBinsGen.size()-1, trackPtBinsGen.data());
    TH2D *hTrackPtResp_SplitMC = new TH2D("hTrackPtResp_SplitMC", "Track Pt Response", trackPtBinsReco.size()-1, trackPtBinsReco.data(), trackPtBinsGen.size()-1, trackPtBinsGen.data());

    // detlat R unfolding hists
    TH1D *hDeltaRSmeared = new TH1D("hDeltaRSmeared", "Delta R Smeared", deltaRBinsReco.size()-1, deltaRBinsReco.data());
    TH1D *hDeltaRGen = new TH1D("hDeltaRGen", "Delta R Gen", deltaRBinsGen.size()-1, deltaRBinsGen.data());
    TH2D *hDeltaRResp = new TH2D("hDeltaRResp", "Delta R Response", deltaRBinsReco.size()-1, deltaRBinsReco.data(), deltaRBinsGen.size()-1, deltaRBinsGen.data());
    TH1D *hDeltaRSmeared_SplitMC = new TH1D("hDeltaRSmeared_SplitMC", "Delta R Smeared", deltaRBinsReco.size()-1, deltaRBinsReco.data());
    TH1D *hDeltaRGen_SplitMC = new TH1D("hDeltaRGen_SplitMC", "Delta R Gen", deltaRBinsGen.size()-1, deltaRBinsGen.data());
    TH2D *hDeltaRResp_SplitMC = new TH2D("hDeltaRResp_SplitMC", "Delta R Response", deltaRBinsReco.size()-1, deltaRBinsReco.data(), deltaRBinsGen.size()-1, deltaRBinsGen.data());

    // energy weighting unfolding hists
    TH1D *hE1E2Smeared = new TH1D("hE1E2Smeared", "E1E2 Smeared", e1e2BinsReco.size()-1, e1e2BinsReco.data());
    TH1D *hE1E2Gen = new TH1D("hE1E2Gen", "E1E2 Gen", e1e2BinsGen.size()-1, e1e2BinsGen.data());
    TH2D *hE1E2Resp = new TH2D("hE1E2Resp", "E1E2 Response", e1e2BinsReco.size()-1, e1e2BinsReco.data(), e1e2BinsGen.size()-1, e1e2BinsGen.data());
    TH1D *hE1E2Smeared_SplitMC = new TH1D("hE1E2Smeared_SplitMC", "E1E2 Smeared", e1e2BinsReco.size()-1, e1e2BinsReco.data());
    TH1D *hE1E2Gen_SplitMC = new TH1D("hE1E2Gen_SplitMC", "E1E2 Gen", e1e2BinsGen.size()-1, e1e2BinsGen.data());
    TH2D *hE1E2Resp_SplitMC = new TH2D("hE1E2Resp_SplitMC", "E1E2 Response", e1e2BinsReco.size()-1, e1e2BinsReco.data(), e1e2BinsGen.size()-1, e1e2BinsGen.data());


    TRandom3* rand = new TRandom3();

    Long64_t totalEvents = pairReader.GetEntries(true);
    for (Long64_t iEvent = 0; iEvent < totalEvents; iEvent++) {
        pairReader.Next();
        for (int iPair = 0; iPair < *nPairs; iPair++) {
            if(RecoE[iPair] < 0 || RecoE2[iPair] < 0) continue; // remove non-matched pairs
            hRecoGenDeltaR->Fill(DistanceReco[iPair], DistanceGen[iPair]);

            // if(DistanceReco[iPair] < 0.5 && DistanceGen[iPair] > 2){
            //     std::cout << "------ Found a mismatch ------" << std::endl;
            //     std::cout << "Reco Delta R: " << DistanceReco[iPair] << " Gen Delta R: " << DistanceGen[iPair] << std::endl;
            //     std::cout << "RecoE: " << RecoE[iPair] << " GenE: " << GenE[iPair] << std::endl;
            //     std::cout << "Distance Reco Gen Particle 1: " << Distance1[iPair] << " Distance Reco-Gen Particle 2: " << Distance2[iPair] << std::endl;
            //     std::cout << " ----------------- " << std::endl;
            // } 
            // get the pT of the reco and gen tracks
            double RecoPt = sqrt(RecoX[iPair]*RecoX[iPair] + RecoY[iPair]*RecoY[iPair] + RecoZ[iPair]*RecoZ[iPair]);
            double GenPt = sqrt(GenX[iPair]*GenX[iPair] + GenY[iPair]*GenY[iPair] + GenZ[iPair]*GenZ[iPair]);
            hRecoGenTrackPt->Fill(RecoPt, GenPt);
            hTrackPtSmeared->Fill(RecoPt);
            hTrackPtGen->Fill(GenPt);
            hDeltaRSmeared->Fill(DistanceReco[iPair]);
            hDeltaRGen->Fill(DistanceGen[iPair]);
            hDeltaRResp->Fill(DistanceReco[iPair], DistanceGen[iPair]);
            hE1E2Smeared->Fill(E1E2Reco[iPair]);
            hE1E2Gen->Fill(E1E2Gen[iPair]);
            hE1E2Resp->Fill(E1E2Reco[iPair], E1E2Gen[iPair]);
            int RecoPtBin = hTrackPtResp->GetXaxis()->FindBin(RecoPt);
            int GenPtBin = hTrackPtResp->GetYaxis()->FindBin(GenPt);
            hTrackPtResp->Fill(RecoPt, GenPt);
            double split = rand->Rndm();
            if(split < 0.5){
                hTrackPtSmeared_SplitMC->Fill(RecoPt);
                hTrackPtGen_SplitMC->Fill(GenPt);
                hDeltaRSmeared_SplitMC->Fill(DistanceReco[iPair]);
                hDeltaRGen_SplitMC->Fill(DistanceGen[iPair]);
                hE1E2Smeared_SplitMC->Fill(E1E2Reco[iPair]);
                hE1E2Gen_SplitMC->Fill(E1E2Gen[iPair]);
            }
            else{
                hTrackPtResp_SplitMC->Fill(RecoPt, GenPt);
                hDeltaRResp_SplitMC->Fill(DistanceReco[iPair], DistanceGen[iPair]);
                hE1E2Resp_SplitMC->Fill(E1E2Reco[iPair], E1E2Gen[iPair]);
            }
            hRecoGenTrackE->Fill(RecoE[iPair], GenE[iPair]);
            hRecoGenEnergyWeighting->Fill(E1E2Reco[iPair], E1E2Gen[iPair]);
            hDeltaDeltaR->Fill(DistanceReco[iPair] - DistanceGen[iPair]);
            if(DistanceGen[iPair] < 0.5) hDelta1->Fill(DistanceReco[iPair] - DistanceGen[iPair]);
            if(DistanceGen[iPair] > 0.5 && DistanceGen[iPair] < 1.0) hDelta2->Fill(DistanceReco[iPair] - DistanceGen[iPair]);
            if(DistanceGen[iPair] > 1.0) hDelta3->Fill(DistanceReco[iPair] - DistanceGen[iPair]);
        }
    } // end loop over the number of events

    // loop over the particles
    Long64_t totalParticles = particleReader.GetEntries(true);
    for (Long64_t iEvent = 0; iEvent < totalParticles; iEvent++) {
        particleReader.Next();
        for (int iParticle = 0; iParticle < *nParticles; iParticle++) {
            singleTrackGenERecoE->Fill(RecoEParticle[iParticle], GenEParticle[iParticle]);
            singleTrackMatchingDistance->Fill(DistanceParticle[iParticle]);
        }
    } // end loop over the number of events

    // now draw the histograms 
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFile* unfoldingRootFile = new TFile(Form("%s/unfoldingHistograms.root",infnameMCDir), "RECREATE");
    unfoldingRootFile->cd();
    hTrackPtSmeared->Write();
    hTrackPtGen->Write();
    hTrackPtResp->Write();
    hTrackPtSmeared_SplitMC->Write();
    hTrackPtGen_SplitMC->Write();
    hTrackPtResp_SplitMC->Write();
    hDeltaRSmeared->Write();
    hDeltaRGen->Write();
    hDeltaRResp->Write();
    hDeltaRSmeared_SplitMC->Write();
    hDeltaRGen_SplitMC->Write();
    hDeltaRResp_SplitMC->Write();
    hE1E2Smeared->Write();
    hE1E2Gen->Write();
    hE1E2Resp->Write();
    hE1E2Smeared_SplitMC->Write();
    hE1E2Gen_SplitMC->Write();
    hE1E2Resp_SplitMC->Write();
    // convert to the unfolding format and write that
    TH1D* hTrackPtSmearedUnfolded = convertHistToUnfoldingFormat(hTrackPtSmeared);
    TH1D* hTrackPtGenUnfolded = convertHistToUnfoldingFormat(hTrackPtGen);
    hTrackPtSmearedUnfolded->Write();
    hTrackPtGenUnfolded->Write();
    unfoldingRootFile->Close();

    // make plot directory
    system(Form("mkdir -p %s/plot", infnameMCDir));

    std::string tag = "MetricVar1";

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetLogz();
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.05);
    hRecoGenDeltaR->SetXTitle("#Delta#it{R}_{reco}");
    hRecoGenDeltaR->SetYTitle("#Delta#it{R}_{gen}");
    hRecoGenDeltaR->Draw("colz");
    c1->SaveAs(Form("%s/plot/RecoGenDeltaR_%s.pdf", infnameMCDir, tag.c_str()));

    c1->SetLogz(false);
    hRecoGenDeltaR->Draw("colz");
    c1->SaveAs(Form("%s/plot/RecoGenDeltaR_%s_linear.pdf", infnameMCDir, tag.c_str()));

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    c2->SetLogz();
    c2->SetTopMargin(0.05);
    hRecoGenTrackPt->SetXTitle("#it{p}_{T, reco} (GeV/#it{c})");
    hRecoGenTrackPt->SetYTitle("#it{p}_{T, gen} (GeV/#it{c})");
    hRecoGenTrackPt->Draw("colz");
    c2->SetRightMargin(0.15);
    c2->SaveAs(Form("%s/plot/RecoGenTrackPt_%s.pdf", infnameMCDir, tag.c_str()));

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
    c3->SetLogz();
    c3->SetRightMargin(0.15);
    c3->SetTopMargin(0.05);
    hRecoGenEnergyWeighting->SetXTitle("#it{E}_{1}#it{E}_{2} reco");
    hRecoGenEnergyWeighting->SetYTitle("#it{E}_{1}#it{E}_{2} gen");
    hRecoGenEnergyWeighting->Draw("colz");
    c3->SaveAs(Form("%s/plot/RecoGenEnergyWeighting_%s.pdf", infnameMCDir, tag.c_str()));

    TCanvas *c4 = new TCanvas("c4", "c4", 800, 800);
    c4->SetLogz();
    c4->SetRightMargin(0.15);
    c4->SetTopMargin(0.05);
    c4->SetLeftMargin(0.15);
    hRecoGenTrackE->SetXTitle("#it{E}_{reco} (GeV)");
    hRecoGenTrackE->SetYTitle("#it{E}_{gen} (GeV)");
    hRecoGenTrackE->Draw("colz");
    c4->SaveAs(Form("%s/plot/RecoGenTrackE_%s.pdf", infnameMCDir, tag.c_str()));

    TCanvas *c5 = new TCanvas("c5", "c5", 800, 800);
    c5->SetRightMargin(0.05);
    c5->SetTopMargin(0.05);
    c5->SetLeftMargin(0.15);
    c5->SetLogy();
    hDeltaDeltaR->SetLineColor(kBlack);
    hDeltaDeltaR->SetMarkerStyle(20);
    hDeltaDeltaR->SetMarkerColor(kBlack);
    hDeltaDeltaR->Scale(1./hDeltaDeltaR->Integral(), "width");
    hDeltaDeltaR->SetXTitle("#Delta#it{R}_{reco} - #Delta#it{R}_{gen}");
    hDeltaDeltaR->SetYTitle("Probability");
    hDeltaDeltaR->GetYaxis()->SetRangeUser(1e-5,20);
    hDeltaDeltaR->Draw();

    hDelta1->SetLineColor(kRed);
    hDelta1->SetMarkerStyle(20);
    hDelta1->SetMarkerColor(kRed);
    hDelta1->Scale(1./hDelta1->Integral(), "width");
    hDelta1->Draw("same");

    hDelta2->SetLineColor(kBlue);
    hDelta2->SetMarkerStyle(20);
    hDelta2->SetMarkerColor(kBlue);
    hDelta2->Scale(1./hDelta2->Integral(), "width");
    hDelta2->Draw("same");

    hDelta3->SetLineColor(kGreen);
    hDelta3->SetMarkerStyle(20);
    hDelta3->SetMarkerColor(kGreen);
    hDelta3->Scale(1./hDelta3->Integral(), "width");
    hDelta3->Draw("same");


    TLegend *leg = new TLegend(0.68, 0.75, 0.93, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hDeltaDeltaR, "Total", "lp");
    leg->AddEntry(hDelta1, "#Delta#it{R}_{gen} < 0.5", "lp");
    leg->AddEntry(hDelta2, "0.5 < #Delta#it{R}_{gen} < 1.0", "lp");
    leg->AddEntry(hDelta3, "#Delta#it{R}_{gen} > 1.0", "lp");
    leg->Draw();

    c5->SaveAs(Form("%s/plot/DeltaDeltaR_%s.pdf", infnameMCDir, tag.c_str()));

  // draw the single track distribution
  TCanvas *c6 = new TCanvas("c6", "c6", 800, 800);
  c6->SetRightMargin(0.15);
  c6->SetTopMargin(0.05);
  c6->SetLeftMargin(0.15);
  c6->SetLogz();
  singleTrackGenERecoE->SetXTitle("Reco E (GeV)");
  singleTrackGenERecoE->SetYTitle("Gen E (GeV)");
  singleTrackGenERecoE->Draw("colz");
  c6->SaveAs(Form("%s/plot/singleTrackGenERecoE_%s.pdf", infnameMCDir, tag.c_str()));

  // draw teh single track matching distance
  TCanvas *c7 = new TCanvas("c7", "c7", 800, 800);
  c7->SetRightMargin(0.15);
  c7->SetTopMargin(0.05);
  c7->SetLeftMargin(0.15);
  c7->SetLogy();
  singleTrackMatchingDistance->SetXTitle("Matching Distance");
  singleTrackMatchingDistance->SetYTitle("Probability");
  singleTrackMatchingDistance->Draw();
  c7->SaveAs(Form("%s/plot/singleTrackMatchingDistance_%s.pdf", infnameMCDir, tag.c_str()));
  


  int NGen = hTrackPtResp->GetNbinsY();
   int NReco = hTrackPtResp->GetNbinsX();
   std::cout << "NGen = " << NGen << " NReco = " << NReco << std::endl;


   for(int iG = 1; iG <= NGen; iG++)
   {
      double N = 0;
      for(int iR = 1; iR <= NReco; iR++){
         if(iR == NReco){
            std::cout << "iR = " << iR << " iG = " << iG << " N = " << N << std::endl;
         }
         N = N + hTrackPtResp->GetBinContent(iR, iG);
      }
   }

    delete hRecoGenDeltaR;
    delete hRecoGenTrackPt;
    delete hRecoGenTrackE;
    delete hRecoGenEnergyWeighting;
    delete hDeltaDeltaR;
    delete hDelta1;
    delete hDelta2;
    delete hDelta3;
    delete singleTrackGenERecoE;
    delete singleTrackMatchingDistance;
    delete singleTrackMatchingDistanceVsEta;
    delete hTrackPtSmeared;
    delete hTrackPtGen;
    delete hTrackPtResp;
    delete hTrackPtSmeared_SplitMC;
    delete hTrackPtGen_SplitMC;
    delete hTrackPtResp_SplitMC;
    delete hDeltaRSmeared;
    delete hDeltaRGen;
    delete hDeltaRResp;
    delete hDeltaRSmeared_SplitMC;
    delete hDeltaRGen_SplitMC;
    delete hDeltaRResp_SplitMC;
    delete hE1E2Smeared;
    delete hE1E2Gen;
    delete hE1E2Resp;
    delete hE1E2Smeared_SplitMC;
    delete hE1E2Gen_SplitMC;
    delete hE1E2Resp_SplitMC;
    delete rand;
    delete unfoldingRootFile;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete leg;
    delete c6;
    delete c7;
}

int main(int argc, char *argv[])
{
    createUnfoldingHistograms("matchingScheme1/");
    createUnfoldingHistograms("matchingScheme2/");
    createUnfoldingHistograms("matchingScheme3/");
    createUnfoldingHistograms("matchingScheme4/");
}