#include <iostream>
#include <vector>
#include <map>
using namespace std;

#include "TTree.h"
#include "TFile.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "alephTrkEfficiency.h"

#include "TCanvas.h"
#include "TH1D.h"

#define MAX 1000
#define MAXPAIR 10000

int main(int argc, char *argv[]);
double MetricAngle(FourVector A, FourVector B);
double MatchingMetric(FourVector A, FourVector B);
void removeOverlappingTracks(std::vector<FourVector> *reco, std::vector<FourVector> *gen, std::vector<int> genPIDs);
int FindBin(double Value, int NBins, double Bins[]); 

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName  = CL.Get("Input");
   string GenTreeName    = CL.Get("Gen", "tgen");
   string OutputFileName = CL.Get("Output");
   double Fraction       = CL.GetDouble("Fraction", 1.00);
   bool isSherpa         = CL.GetBool("IsSherpa", false); 
   int multCut           = CL.GetInt("MultCut", 100000); 

   TFile* InputFile  = new TFile(InputFileName.c_str());
   TFile* OutputFile = new TFile(OutputFileName.c_str(), "RECREATE");

 

   // z binning
   const int BinCount = 100;
   double zBins[2*BinCount+1]; 
   double zBinMin = (1- cos(0.002))/2; 
   double zBinMax = 0.5;

   for(int i = 0; i <= BinCount; i++){
      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

   }

   // double log binning
   TH1D* genUnmatched_z = new TH1D("genUnmatched_z", "genUnmatched_z", 2 * BinCount, 0, 2 * BinCount);
   ParticleTreeMessenger* MGen = new ParticleTreeMessenger(InputFile, GenTreeName);
   TH1D HN("HN", ";;", 1, 0, 1);

   int EntryCount = MGen->GetEntries() * Fraction;
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1); 
   int nAcceptedEvents = 0; 
   double TotalE = 91.1876; // GeV
   for(int iE = 0; iE < EntryCount; iE++) 
   {

      MGen->GetEntry(iE);
      if(MGen->nChargedHadronsHP < multCut) continue;



      nAcceptedEvents++; 

      vector<FourVector> PGen;
      for(int i = 0; i < MGen->nParticle; i++){
         // charged particle selection 
         if(MGen->charge[i] == 0) continue;
         if(MGen->highPurity[i] == false) continue;
         PGen.push_back(MGen->P[i]);
      }
      //std::cout <<PGen.size() << std::endl;


      for(int i = 0; i < PGen.size(); i++){
      for(int j = i+1; j < PGen.size();j++){
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
           
            // z histograms
            double zGenUnmatched = (1-cos(GetAngle(Gen1, Gen2)))/2; 
            int BinZGen = FindBin(zGenUnmatched, 2*BinCount, zBins); 
            genUnmatched_z->Fill(BinZGen, Gen1[0]*Gen2[0]/(TotalE*TotalE));
           
         }
      }

   }

   std::cout << "The number of accepted events is " << nAcceptedEvents << std::endl;
   HN.SetBinContent(1, nAcceptedEvents);
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   OutputFile->cd();
   genUnmatched_z->Write();
   HN.Write(); 
   OutputFile->Close();
   InputFile->Close();

   return 0;
}



int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}
