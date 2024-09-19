#include <iostream>
#include <vector>
#include <map>
using namespace std;

// root includes
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "SetStyle.h"



int main(int argc, char *argv[]);
int FindBin(double Value, int NBins, double Bins[]); 
void MakeCanvasZ(vector<TH1D > Histograms, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void MakeCanvas(vector<TH1D> Histograms, TGraphErrors DataSyst, vector<string> Labels, string Output, string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX); 
void SetPad(TPad &P); 
void DivideByBin(TH1D &H, double Bins[]); 

// can we put this in the CommonCode/ ?
void FillChain(TChain &chain, const vector<string> &files) {
  for (auto file : files) {
    chain.Add(file.c_str());
  }
}

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   SetThesisStyle();
   static vector<int> Colors = GetCVDColors6();


   string InputDataPath          = CL.Get("InputData");
   vector<string> InputMCPaths   = CL.GetStringVector("InputMCs");
   vector<string> InputMCLegs    = CL.GetStringVector("InputMCLegs");
   vector<string> InputSystPaths = CL.GetStringVector("InputSysts");
   vector<string> InputSystNames = CL.GetStringVector("InputSystNames", vector<string>());
   string chooseSyst             = CL.Get("chooseSyst", "All");
   string OutputFileName         = CL.Get("Output", "CompareDataMCEEC.root");

   bool InputSystNamesInitialized= InputSystPaths.size()==InputSystNames.size();
   if (!InputSystNamesInitialized) 
   {
      printf("[INFO] InputSystPaths.size()=%d is not equal to InputSystNames.size()=%d. Will set the InputSystNames to a default enumerated setting. \n", 
               InputSystPaths.size(), InputSystNames.size());
      InputSystNames.clear();
   }

   // 1D histograms
   TFile InputData(InputDataPath.c_str());
   TH1D HEEC2Data( *((TH1D*) InputData.Get("HEEC2")) ); HEEC2Data.SetName("Data");
   TH1D HNData( *((TH1D*) InputData.Get("HN")) );

   auto CopyTH1D = [](TH1D& h_target, TH1D& h_source)
   {
      h_target.SetBins(h_source.GetNbinsX(), h_source.GetXaxis()->GetXmin(), h_source.GetXaxis()->GetXmax());

      for (int i = 1; i <= h_source.GetNbinsX(); ++i) 
      {
        double binContent = h_source.GetBinContent(i);
        double binError = h_source.GetBinError(i);
        
        h_target.SetBinContent(i, binContent);
        h_target.SetBinError(i, binError);
      }
   };

   vector<TH1D> HEEC2MCs;
   vector<TH1D> HNMCs;
   for (int iMCFile = 0; iMCFile < InputMCPaths.size(); ++iMCFile)
   { 
      TFile InputMC(InputMCPaths[iMCFile].c_str());
      TH1D _HEEC2MC( *((TH1D*) InputMC.Get("HEEC2")) );
      TH1D _HNMC( *((TH1D*) InputMC.Get("HN")) );

      HEEC2MCs.push_back( TH1D() );
      HNMCs.push_back( TH1D() );
      
      CopyTH1D(HEEC2MCs[iMCFile], _HEEC2MC);
      CopyTH1D(HNMCs[iMCFile], _HNMC);

      HEEC2MCs[iMCFile].SetName(InputMCLegs[iMCFile].c_str());
   }

   vector<TH1D> HEEC2Systs;
   vector<TH1D> HNSysts;
   for (int iSystFile = 0; iSystFile < InputSystPaths.size(); ++iSystFile)
   {
      TFile InputSyst(InputSystPaths[iSystFile].c_str());
      TH1D _HEEC2Syst( *((TH1D*) InputSyst.Get("HEEC2")) );
      TH1D _HNSyst( *((TH1D*) InputSyst.Get("HN")) );

      HEEC2Systs.push_back( TH1D() );
      HNSysts.push_back( TH1D() );
      
      CopyTH1D(HEEC2Systs[iSystFile], _HEEC2Syst);
      CopyTH1D(HNSysts[iSystFile], _HNSyst);

      if (!InputSystNamesInitialized) InputSystNames.push_back(Form("Syst%d", iSystFile+1));

      HEEC2Systs[iSystFile].SetName(InputSystNames[iSystFile].c_str());
   }

   // divide by number of events
   HEEC2Data.Scale(1/HNData.GetBinContent(1));

   HEEC2Data.SetMarkerColor(Colors[2]);
   HEEC2Data.SetLineColor(Colors[2]);
   HEEC2Data.SetMarkerStyle(20);
   HEEC2Data.SetLineWidth(2);

   for (int iMCFile = 0; iMCFile < InputMCPaths.size(); ++iMCFile)
   {
      HEEC2MCs[iMCFile].Scale(1/HNMCs[iMCFile].GetBinContent(1));
      HEEC2MCs[iMCFile].SetMarkerColor(Colors[3+iMCFile]);
      HEEC2MCs[iMCFile].SetLineColor(Colors[3+iMCFile]);
      HEEC2MCs[iMCFile].SetMarkerStyle(20);
      HEEC2MCs[iMCFile].SetLineWidth(2);

      // HEEC2MCs[iMCFile].Print("all");
   }

   TGraphErrors GEEC2DataSyst( &HEEC2Data ); GEEC2DataSyst.SetName("DataSyst");
   for (int iSystFile = 0; iSystFile < InputSystPaths.size(); ++iSystFile) HEEC2Systs[iSystFile].Scale(1/HNSysts[iSystFile].GetBinContent(1));
   
   for(int i = 1; i <= HEEC2Data.GetNbinsX(); i++)
   {
      int iGraph = i-1;
      // printf("%d/%d %.3f %.3f\n", i, HEEC2Data.GetNbinsX(),
      //        HEEC2Data.GetBinCenter(i),
      //        GEEC2DataSyst.GetPointX(iGraph));

      double systUncertInCurrentBin = 0;
      for (auto ele: HEEC2Systs)
      {
         if (chooseSyst=="All")
         {
            systUncertInCurrentBin += TMath::Power(ele.GetBinContent(i) - HEEC2Data.GetBinContent(i), 2);
         }
         else if (chooseSyst==ele.GetName())
         {
            systUncertInCurrentBin = TMath::Abs(ele.GetBinContent(i) - HEEC2Data.GetBinContent(i));
            break;
         }
      }

      if (chooseSyst=="All") systUncertInCurrentBin = TMath::Sqrt(systUncertInCurrentBin);

      GEEC2DataSyst.SetPointError(iGraph, 0.5, systUncertInCurrentBin);
   }
   // HEEC2Data.Print("all");
   // GEEC2DataSyst.Print("all");
   GEEC2DataSyst.SetLineWidth(0);
   GEEC2DataSyst.SetFillStyle(1001);
   GEEC2DataSyst.SetFillColorAlpha(Colors[2], 0.3);

   vector<TH1D> hists_theta(HEEC2MCs);
   vector<string> legends(InputMCLegs);
   hists_theta.insert( hists_theta.begin(), HEEC2Data );
   legends.insert( legends.begin(), "Data" );

   system("mkdir -p plot/");
   MakeCanvas(hists_theta, GEEC2DataSyst, legends, "plot/EEC2","#theta", "#frac{1}{#it{N}_{event}}#frac{d(Sum E_{i}E_{j}/E^{2})}{d#it{#theta}}",1e-3, 2e0, true, true);

   // write to the output file
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");
   OutputFile.cd(); 
   HEEC2Data.Write(); 
   for(auto ele: HEEC2MCs) ele.Write();

   // cleanup
   InputData.Close(); 

   return 0;
}


int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}

void MakeCanvasZ(vector<TH1D > Histograms, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX){
   

   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

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
   for(TH1D H : Histograms){
      TH1D *HClone = (TH1D *)H.Clone();
      HClone->Draw("exp same");
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
         TH1D *H = (TH1D *)Histograms[i].Clone();
         H->Divide(&Histograms[0]);
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
      Latex.DrawLatex(MarginL * 0.3, MarginB + PadRHeight * 0.5, "Matching Eff.");
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
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
   Legend.Draw();

   TLegend Legend2(0.55, 0.90, 0.8, 0.90 - 0.04 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), "pl");
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


void MakeCanvas(vector<TH1D> Histograms,  TGraphErrors DataSyst, vector<string> Labels, string Output,
   string X, string Y, double WorldMin, double WorldMax, bool DoRatio, bool LogX)
{
   int NLine = Histograms.size();
   int N = Histograms[0].GetNbinsX();

   double MarginL = 180;
   double MarginR = 90;
   double MarginB = 120;
   double MarginT = 90;

   double WorldXMin = LogX ? 0 : 0;
   double WorldXMax = LogX ? N : M_PI;
   
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
   for(TH1D H : Histograms){
      TH1D *H1 = (TH1D *)H.Clone();
      H1->Draw("same");
   }
   // DataSyst.Print("all");
   DataSyst.DrawClone("2 same");

   TGraph G;
   G.SetPoint(0, LogX ? N / 2 : M_PI / 2, 0);
   G.SetPoint(1, LogX ? N / 2 : M_PI / 2, 1000);
   G.SetLineStyle(kDashed);
   G.SetLineColor(kGray);
   G.SetLineWidth(1);
   G.Draw("l");

   if(DoRatio)
      PadR.cd();

   double WorldRMin = 0.6;
   double WorldRMax = 1.4;
   
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
         TH1D *H = (TH1D *)Histograms[i].Clone();
         H->Divide(&Histograms[0]);
         H->Draw("same");
      }
      
      for(int i = 1; i <= Histograms[0].GetNbinsX(); i++)
      {
         int iGraph = i-1;
         DataSyst.SetPoint(iGraph, 
                                 DataSyst.GetPointX(iGraph), 
                                 DataSyst.GetPointY(iGraph)/Histograms[0].GetBinContent(i));
         DataSyst.SetPointError( iGraph, 
                                 DataSyst.GetErrorX(iGraph), 
                                 DataSyst.GetErrorY(iGraph)/Histograms[0].GetBinContent(i));
      }
      DataSyst.DrawClone("2 same");

      G.Draw("l");

      G2.SetPoint(0, 0, 1);
      G2.SetPoint(1, 99999, 1);
      G2.Draw("l");
   }
   
   double BinMin    = 0.002;
   double BinMiddle = M_PI / 2;
   double BinMax    = M_PI - 0.002;

   Canvas.cd();
   TGaxis X1(MarginL, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "GS");
   TGaxis X2(MarginL + PadWidth, MarginB, MarginL + PadWidth / 2, MarginB, BinMin, BinMiddle, 510, "-GS");
   TGaxis X3(MarginL, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis X4(MarginL + PadWidth, MarginB + PadRHeight, MarginL + PadWidth / 2, MarginB + PadRHeight, BinMin, BinMiddle, 510, "+-GS");
   TGaxis Y1(MarginL, MarginB, MarginL, MarginB + PadRHeight, WorldRMin, WorldRMax, 505, "");
   TGaxis Y2(MarginL, MarginB + PadRHeight, MarginL, MarginB + PadRHeight + PadHeight, WorldMin, WorldMax, 510, "G");
   
   TGaxis XL1(MarginL, MarginB, MarginL + PadWidth, MarginB, 0, M_PI, 510, "S");
   TGaxis XL2(MarginL, MarginB + PadRHeight, MarginL + PadWidth, MarginB + PadRHeight, 0, M_PI, 510, "+-S");

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
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.115, MarginB - 0.01, "0.01");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.290, MarginB - 0.01, "0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.465, MarginB - 0.01, "1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.535, MarginB - 0.01, "#pi - 1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.710, MarginB - 0.01, "#pi - 0.1");
   if(LogX) Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.885, MarginB - 0.01, "#pi - 0.01");

   Latex.SetTextAlign(12);
   Latex.SetTextAngle(270);
   Latex.SetTextColor(kGray);
   Latex.DrawLatex(MarginL + (1 - MarginR - MarginL) * 0.5 + 0.0175, 1 - MarginT - 0.015, "#theta_{L} = #pi/2");

   Latex.SetTextAlign(22);
   Latex.SetTextAngle(0);
   Latex.SetTextColor(kBlack);
   Latex.DrawLatex(MarginL + PadWidth * 0.5, MarginB * 0.3, X.c_str());

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
   Latex.DrawLatex(0.01, 0.01, "Work-in-progress Yi Chen + HB/YJL/AB/..., 2024 Feb 12 (26436_FullEventEEC)");

   TLegend Legend(0.15, 0.90, 0.35, 0.90 - 0.05 * min(NLine, 4));
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   for(int i = 0; i < NLine && i < 4; i++)
   {
      if (Labels[i]=="Data")
      {
         Histograms[i].SetFillStyle(DataSyst.GetFillStyle());
         Histograms[i].SetFillColor(DataSyst.GetFillColor());
      }
      Legend.AddEntry(&Histograms[i], Labels[i].c_str(), 
                      (Labels[i]=="Data")? "plf": "pl");
   }
   Legend.Draw();

   TLegend Legend2(0.7, 0.90, 0.9, 0.90 - 0.05 * (NLine - 4));
   Legend2.SetTextFont(42);
   Legend2.SetTextSize(0.035);
   Legend2.SetFillStyle(0);
   Legend2.SetBorderSize(0);
   if(NLine >= 4)
   {
      for(int i = 4; i < NLine; i++)
         Legend2.AddEntry(&Histograms[i], Labels[i].c_str(), 
                      (Labels[i]=="Data")? "plf": "pl");
      Legend2.Draw();
   }

   Canvas.SaveAs((Output + ".pdf").c_str());
}