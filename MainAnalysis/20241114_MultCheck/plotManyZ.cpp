/* plotFromUnfolding.cpp: Macro to plot the final result from the unfolded result.
 * Hannah Bossi, <hannah.bossi@cern.ch>
 * 05/28/2024
 */

#include <iostream>
#include <iomanip> 
#include <limits> 
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TGraph.h"

#include "SetStyle.h"
#include "ProgressBar.h"
#include "CommandLine.h"
#include "Messenger.h"
#include "JetCorrector.h"
#include "alephTrkEfficiency.h"




void MakeCanvasZ(vector<TH1D *> Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 

void SetPad(TPad &P); 
void DivideByBin(TH1D &H, double Bins[]); 
int FindBin(double Value, int NBins, double Bins[]); 
int main(int argc, char *argv[]); 

int main(int argc, char *argv[]){

    // take the parameters from the command line
    CommandLine CL(argc, argv);
    vector<string> FileNames = CL.GetStringVector("Input");
    string output = CL.Get("Output", "plots");
    vector<string> Labels = CL.GetStringVector("Label");
    string Prefix = CL.Get("Prefix", "");
    bool DoRatio = CL.GetBool("DoRatio", true);
    bool DoReflection = CL.GetBool("DoReflection", true);

 
    // SetThesisStyle();
    static vector<int> Colors = GetCVDColors6();


   // theta binning
   const int BinCount = 100;
   double Bins[2*BinCount+1];
   double BinMin = 0.002;
   double BinMax = M_PI / 2;

   // z binning
   double zBins[2*BinCount+1]; 
   double zBinMin = (1- cos(0.002))/2; 
   double zBinMax = 0.5;

   // energy binning
   double EnergyBinMin = 4e-6; 
   double EnergyBinMax = 0.2; 
   double EnergyBins[BinCount+1];
   double logMin = std::log10(EnergyBinMin);
   double logMax = std::log10(EnergyBinMax);
   double logStep = (logMax - logMin) / (BinCount);

   for(int i = 0; i <= BinCount; i++)
   {
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      double logValue = logMin + i * logStep;
      EnergyBins[i] =  std::pow(10, logValue);      //EnergyBins[i] =  exp(log(EnergyBinMin) + (log(EnergyBinMax) - log(EnergyBinMin)) / BinCount * i);
       // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }


    double TotalE = 91.1876; // GeV




    vector<TH1D *> Histograms;

    int N = FileNames.size();

    vector<TFile *> Files(N);
    for(int i = 0; i < N; i++){
      std::cout << "Opening file... " << FileNames[i].c_str() << std::endl;
      Files[i] = new TFile(FileNames[i].c_str());
    }

    for(int i = 0; i < N; i++){
        TH1D* H = (TH1D*)Files[i]->Get("genUnmatched_z");
        TH1D *HN = (TH1D *)Files[i]->Get("HN");

        // divide by the number of events 
        DivideByBin(*H, zBins);
        std::cout << "Number of events " <<  HN->GetBinContent(1) << std::endl;
        std::cout << "Number of entires " << H->GetEntries() << std::endl;
        H->Scale(1.0 /  HN->GetBinContent(1));

    
        H->SetMarkerColor(Colors[i]);
        H->SetLineColor(Colors[i]);
        H->SetMarkerStyle(20);
        H->SetLineWidth(2);

        Histograms.push_back(H);

        // if we do the reflection, also push back another histogram
        TH1D *H2= (TH1D*)H->Clone("H2"); 
        H2->Reset(); 
        if(DoReflection){

            // figure out the bin for the pivot point
            double pivotPoint = 0.5;
            int pivotBin = FindBin(pivotPoint,2*BinCount+1, zBins);

            // now create the flipped histogram
            // Fill the new histogram with the flipped data
            for (int bin = 1; bin <= H->GetNbinsX(); ++bin) {
               // modified to use a pivot point instead of a pivot bin 
               int flippedBin = (2 * pivotBin) + 1 - bin; 
               if (flippedBin >= 1 && flippedBin <= H->GetNbinsX()) {
                    H2->SetBinContent(bin, H->GetBinContent(flippedBin));
                }
            }
            H2->SetMarkerColor(Colors[i+1]);
            H2->SetLineColor(Colors[i+1]);
            H2->SetLineWidth(2);
            H2->SetMarkerStyle(24);
            Histograms.push_back(H2);
            Labels.push_back("Reflection");
        }


    }
 

    double YMax = DoRatio ? 10 : 1;

    // now have a different call in case we are plotting as a function of the theta
    MakeCanvasZ(Histograms, Labels, output + "Z" + Prefix + "EEC2","#it{z} = (1- cos(#theta))/2", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{z}}",2e-3, 3000,
    DoRatio, true);
   
    for(int i = 0; i < N; i++){
        Files[i]->Close();
        delete Files[i];
    }

    return 0;

}

void MakeCanvasZ(vector<TH1D *> Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX){
   

   int NLine = Histograms.size();
   int N = Histograms[0]->GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : 1;
   
   double PadWidth = 1200;
   double PadHeight = DoRatio ? 640 : 640 + 240;
   double PadRHeight = DoRatio ? 240 : 0.001;

   double CanvasWidth = MarginL + PadWidth + MarginR;
   double CanvasHeight = MarginT + PadHeight + PadRHeight + MarginB;

   MarginL = MarginL / CanvasWidth;
   MarginR = MarginR / CanvasWidth;
   MarginT = MarginT / CanvasHeight;
   MarginB = MarginB / CanvasHeight;

   PadWidth   = PadWidth / CanvasWidth;
   PadHeight  = PadHeight / CanvasHeight;
   PadRHeight = PadRHeight / CanvasHeight;

   TCanvas Canvas("Canvas", "", CanvasWidth, CanvasHeight);
   // Canvas.SetLogy();
   // Canvas.SetRightMargin(MarginR);
   // Canvas.SetLeftMargin(MarginL);
   // Canvas.SetTopMargin(MarginT);
   // Canvas.SetBottomMargin(MarginB);
   Canvas.cd(); 

   TPad Pad("Pad", "", MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadHeight + PadRHeight);
   Pad.SetLogy();
   SetPad(Pad);
   
   TPad PadR("PadR", "", MarginL, MarginB, MarginL + PadWidth, MarginB + PadRHeight);
   if(DoRatio)
      SetPad(PadR);

   Pad.cd();

   TH2D HWorld("HWorld", "", N, WorldXMin, WorldXMax, 100, WorldMin, WorldMax);
   HWorld.SetStats(0);
   HWorld.GetXaxis()->SetTickLength(0);
   HWorld.GetXaxis()->SetLabelSize(0);

   HWorld.Draw("axis");
   for(TH1D *H : Histograms){
        H->Draw("exp same");
        std::cout << "On histogram " << H->GetName() << " which has an x axis range of " << H->GetXaxis()->GetXmin() << " to " << H->GetXaxis()->GetXmax() << std::endl;
   }

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : 1 / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : 1/ 2, 10000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.5;
   double WorldRMax = 1.5;
   
   TH2D HWorldR("HWorldR", "", N, WorldXMin, WorldXMax, 100, WorldRMin, WorldRMax);
   TGraph G2;
   
   if(DoRatio)
   {
      HWorldR.SetStats(0);
      HWorldR.GetXaxis()->SetTickLength(0);
      HWorldR.GetXaxis()->SetLabelSize(0);
      HWorldR.GetYaxis()->SetNdivisions(505);

      HWorldR.Draw("axis");
      for(int i = 1; i < NLine; i++)
      {
         TH1D *H = (TH1D *)Histograms[i]->Clone();
         H->Divide(Histograms[0]);
         H->Draw("same");
      }

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }
   
   double BinMin    = (1- cos(0.002))/2;
   double BinMiddle = 0.5;
   double BinMax    = 1 - BinMin;

   Canvas.cd();
   std::cout << "Here 1" << std::endl;
   int nDiv = 505;
   std::cout << "Bin Min " << BinMin << " Bin Middle " << BinMiddle << std::endl;
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, nDiv, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, nDiv, "+-GS");

   std::cout << "Here 2" << std::endl;
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");

   std::cout << "Here 3" << std::endl;
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB,  (1- cos(0.002))/2, 1, 210, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight,  (1- cos(0.002))/2, 1, 210, "+-S");

   Y1.SetLabelFont(42);
   Y2.SetLabelFont(42);
   XL1.SetLabelFont(42);
   XL2.SetLabelFont(42);

   X1.SetLabelSize(0);
   X2.SetLabelSize(0);
   X3.SetLabelSize(0);
   X4.SetLabelSize(0);
   // XL1.SetLabelSize(0);
   XL2.SetLabelSize(0);

   X1.SetTickSize(0.06);
   X2.SetTickSize(0.06);
   X3.SetTickSize(0.06);
   X4.SetTickSize(0.06);
   XL1.SetTickSize(0.03);
   XL2.SetTickSize(0.03);

   if(LogX == true)
   {
      X1.Draw();
      X2.Draw();
      if(DoRatio) X3.Draw();
      if(DoRatio) X4.Draw();
   }
   if(LogX == false)
   {
      XL1.Draw();
      if(DoRatio)
         XL2.Draw();
   }
   if(DoRatio)
      Y1.Draw();
   Y2.Draw();

   TLatex Latex;
   Latex.SetNDC();
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.035);
   Latex.SetTextAlign(23);
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.02, MarginB - 0.01, "10^{-6} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.180, MarginB - 0.01, "10^{-4} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.365, MarginB - 0.01, "10^{-2} ");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.500, MarginB - 0.01, "1/2");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.650, MarginB - 0.01, "1 - 10^{-2}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.823, MarginB - 0.01, "1 - 10^{-4}");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.995, MarginB - 0.01, "1 - 10^{-6}");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#it{z} = 1/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.9, MarginB * 0.4, X.c_str());

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(90);
   Latex.SetTextColor(kBlack);
   if(DoRatio)
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Ratio");
   Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight + PadHeight * 0.5, Y.c_str());

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.DrawLatex(MarginL, MarginB + PadRHeight + PadHeight + 0.012, "ALEPH e^{+}e^{-}, #sqrt{s} = 91.2 GeV, Work-in-progress");

   Latex.SetTextAlign(11);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(19);
   Latex.SetTextSize(0.02);
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress, 2024 August 7th HB");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.04 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
      Legend.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(Histograms[i], Labels[i].c_str(), "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}




void SetPad(TPad &P){
   P.SetLeftMargin(0);
   P.SetTopMargin(0);
   P.SetRightMargin(0);
   P.SetBottomMargin(0);
   P.SetTickx();
   P.SetTicky();
   P.Draw();
}



void DivideByBin(TH1D &H, double Bins[])
{
   int N = H.GetNbinsX();
   for(int i = 1; i <= N; i++)
   {
      double L = Bins[i-1];
      double R = Bins[i];
      H.SetBinContent(i, H.GetBinContent(i) / (R - L));
      H.SetBinError(i, H.GetBinError(i) / (R - L));
   }
}


int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(float(Value) < float(Bins[i])){
         return i - 1;
      }
   return NBins;
}