#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "Matching.h"
#include "ProgressBar.h"
#include "TauHelperFunctions3.h"
#include "alephTrkEfficiency.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLatex.h"

#define MAX 1000
#define MAXPAIR 10000

int main(int argc, char *argv[]);
double MetricAngle(FourVector A, FourVector B);

// [Warning] A global variable used in *this* cpp script! //
int MatchingSchemeChoice = 2;
double MatchingMetricCore( double deltaTheta, double deltaPhi, double deltaE,
                           double meanE,     // the average btw the GenE and RecoE,
                                             // this would be used to model the energy-dependency in the resolution assignment
                           int schemeChoice, // 1: 10% energy as resolution, flat theta and phi resolutions
                                             // 2: resolutions are cast according to ALEPH tracker performance (energy-dependent description)
                                             // 3: scale btw the importance of angular-match versus energy-match ( 5x)
                                             // 4: scale btw the importance of angular-match versus energy-match (15x)
                           double& chiTheta, double& chiPhi, double& chiE);
double MatchingMetric(FourVector A, FourVector B);
int FindBin(double Value, int NBins, double Bins[]);
void MatchingPerformance(string& matchedRstRoot, string& rstDirName,
                         ParticleTreeMessenger& MGen,
                         ParticleTreeMessenger& MReco);

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   string InputFileName  = CL.Get("Input");
   string GenTreeName    = CL.Get("Gen", "tgen");
   string RecoTreeName   = CL.Get("Reco", "t");
   string OutputFileName = CL.Get("Output");
   double Fraction       = CL.GetDouble("Fraction", 1.00);
   // [Warning] A global variable used in *this* cpp script! //
   MatchingSchemeChoice  = CL.GetInt("MatchingSchemeChoice", 2);

   if (MatchingSchemeChoice!=1 &&
       MatchingSchemeChoice!=2 &&
       MatchingSchemeChoice!=3 &&
       MatchingSchemeChoice!=4 )
   {
      printf("MatchingSchemeChoice %d not supported. Exiting\n", MatchingSchemeChoice);
      exit(1);
   }
   string rstDirName     = (MatchingSchemeChoice==1)? "matchingScheme1/":
                           (MatchingSchemeChoice==2)? "matchingScheme2/":
                           (MatchingSchemeChoice==3)? "matchingScheme3/": "matchingScheme4/";
   system(("mkdir -p "+rstDirName).c_str());
   OutputFileName = rstDirName + OutputFileName;

   TFile InputFile(InputFileName.c_str());
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   //------------------------------------
   // define the binning
   //------------------------------------

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
   double EnergyBins[BinCount+1];
   double EnergyBinMin = 4e-6;
   double EnergyBinMax = 0.2;
   double logMin = std::log10(EnergyBinMin);
   double logMax = std::log10(EnergyBinMax);
   double logStep = (logMax - logMin) / (BinCount);

   for(int i = 0; i <= BinCount; i++){
      // theta double log binning
      Bins[i] = exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);
      Bins[2*BinCount-i] = BinMax * 2 - exp(log(BinMin) + (log(BinMax) - log(BinMin)) / BinCount * i);

      // z double log binning
      zBins[i] = exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);
      zBins[2*BinCount-i] = zBinMax * 2 - exp(log(zBinMin) + (log(zBinMax) - log(zBinMin)) / BinCount * i);

      double logValue = logMin + i * logStep;
      EnergyBins[i] =  std::pow(10, logValue);//exp(log(EnergyBinMin) + (log(EnergyBinMax) - log(EnergyBinMin)) / BinCount * i);
      //std::cout << "Adding energy bin " << EnergyBins[i] << std::endl;
   }

   //------------------------------------
   // define the trees
   //------------------------------------

   // tree for single track matching
   TTree OutputTree("MatchedTree", "");
   // tree for matched pairs
   TTree OutputPairTree("PairTree", "");
   // tree that just takes pairs, not matching taken into account
   TTree OutputUnmatchedTree("UnmatchedPairTree", "");

   //------------------------------------
   // define the histograms
   //------------------------------------


   // as a function of theta
   TH1D recoUnmatched("recoUnmatched", "recoUnmatched", 2 * BinCount, 0, 2 * BinCount);
   TH1D genUnmatched("genUnmatched", "genUnmatched", 2 * BinCount, 0, 2 * BinCount);
   TH1D recoMatched("recoMatched", "recoMatched", 2 * BinCount, 0, 2 * BinCount);
   TH1D genMatched("genMatched", "genMatched", 2 * BinCount, 0, 2 * BinCount);

   // as a function of z
   TH1D recoUnmatched_z("recoUnmatched_z", "recoUnmatched_z", 2 * BinCount, 0, 2 * BinCount);
   TH1D genUnmatched_z("genUnmatched_z", "genUnmatched_z", 2 * BinCount, 0, 2 * BinCount);
   TH1D recoMatched_z("recoMatched_z", "recoMatched_z", 2 * BinCount, 0, 2 * BinCount);
   TH1D genMatched_z("genMatched_z", "genMatched_z", 2 * BinCount, 0, 2 * BinCount);

   // as a function of e1e2
   TH1D e1e2RecoUnmatched("e1e2RecoUnmatched", "e1e2RecoUnmatched", BinCount, EnergyBins);
   TH1D e1e2RecoMatched("e1e2RecoMatched", "e1e2RecoMatched", BinCount, EnergyBins);
   TH1D e1e2GenUnmatched("e1e2GenUnmatched", "e1e2GenUnmatched", BinCount, EnergyBins);
   TH1D e1e2GenMatched("e1e2GenMatched", "e1e2GenMatched", BinCount, EnergyBins);


   alephTrkEfficiency efficiencyCorrector;
   // variables for the matched tree
   int NParticle, eventID;
   double GenE[MAX], GenX[MAX], GenY[MAX], GenZ[MAX];
   double RecoE[MAX], RecoX[MAX], RecoY[MAX], RecoZ[MAX];
   double Distance[MAX], DeltaPhi[MAX], Metric[MAX], DeltaE[MAX], DeltaTheta[MAX];
   double RecoEta[MAX], RecoPhi[MAX], RecoTheta[MAX], RecoRapidity[MAX];
   int RecoPWFlag[MAX];
   double GenEta[MAX], GenPhi[MAX], GenTheta[MAX], GenRapidity[MAX];
   double RecoEfficiency[MAX];
   OutputTree.Branch("EventID", &eventID, "EventID/I");
   OutputTree.Branch("NParticle", &NParticle, "NParticle/I");
   OutputTree.Branch("GenE", &GenE, "GenE[NParticle]/D");
   OutputTree.Branch("GenX", &GenX, "GenX[NParticle]/D");
   OutputTree.Branch("GenY", &GenY, "GenY[NParticle]/D");
   OutputTree.Branch("GenZ", &GenZ, "GenZ[NParticle]/D");
   OutputTree.Branch("RecoE", &RecoE, "RecoE[NParticle]/D");
   OutputTree.Branch("RecoX", &RecoX, "RecoX[NParticle]/D");
   OutputTree.Branch("RecoY", &RecoY, "RecoY[NParticle]/D");
   OutputTree.Branch("RecoZ", &RecoZ, "RecoZ[NParticle]/D");
   OutputTree.Branch("Distance", &Distance, "Distance[NParticle]/D");
   OutputTree.Branch("DeltaPhi", &DeltaPhi, "DeltaPhi[NParticle]/D");
   OutputTree.Branch("DeltaTheta", &DeltaTheta, "DeltaTheta[NParticle]/D");
   OutputTree.Branch("DeltaE", &DeltaE, "DeltaE[NParticle]/D");
   OutputTree.Branch("Metric", &Metric, "Metric[NParticle]/D");
   OutputTree.Branch("RecoEta", &RecoEta, "RecoEta[NParticle]/D");
   OutputTree.Branch("RecoPhi", &RecoPhi, "RecoPhi[NParticle]/D");
   OutputTree.Branch("RecoTheta", &RecoTheta, "RecoTheta[NParticle]/D");
   OutputTree.Branch("GenEta", &GenEta, "GenEta[NParticle]/D");
   OutputTree.Branch("GenPhi", &GenPhi, "GenPhi[NParticle]/D");
   OutputTree.Branch("GenTheta", &GenTheta, "GenTheta[NParticle]/D");
   OutputTree.Branch("RecoRapidity", &RecoRapidity, "RecoRapidity[NParticle]/D");
   OutputTree.Branch("GenRapidity", &GenRapidity, "GenRapidity[NParticle]/D");
   OutputTree.Branch("RecoPwFlag", &RecoPWFlag, "RecoPWFlag[NParticle]/I");
   OutputTree.Branch("RecoEfficiency", &RecoEfficiency, "RecoEfficiency[NParticle]/D");

   // valiables for the pair tree
   int NPair;
   double GenE1[MAXPAIR], GenX1[MAXPAIR], GenY1[MAXPAIR], GenZ1[MAXPAIR];
   double GenE2[MAXPAIR], GenX2[MAXPAIR], GenY2[MAXPAIR], GenZ2[MAXPAIR];
   double RecoE1[MAXPAIR], RecoX1[MAXPAIR], RecoY1[MAXPAIR], RecoZ1[MAXPAIR];
   double RecoE2[MAXPAIR], RecoX2[MAXPAIR], RecoY2[MAXPAIR], RecoZ2[MAXPAIR];
   double DistanceGen[MAXPAIR], DistanceReco[MAXPAIR], Distance1[MAXPAIR], Distance2[MAXPAIR];
   double E1E2Gen[MAXPAIR], E1E2Reco[MAXPAIR];
   double RecoEfficiency1[MAXPAIR], RecoEfficiency2[MAXPAIR];
   OutputPairTree.Branch("NPair", &NPair);
   OutputPairTree.Branch("GenE1", &GenE1, "GenE1[NPair]/D");
   OutputPairTree.Branch("GenX1", &GenX1, "GenX1[NPair]/D");
   OutputPairTree.Branch("GenY1", &GenY1, "GenY1[NPair]/D");
   OutputPairTree.Branch("GenZ1", &GenZ1, "GenZ1[NPair]/D");
   OutputPairTree.Branch("GenE2", &GenE2, "GenE2[NPair]/D");
   OutputPairTree.Branch("GenX2", &GenX2, "GenX2[NPair]/D");
   OutputPairTree.Branch("GenY2", &GenY2, "GenY2[NPair]/D");
   OutputPairTree.Branch("GenZ2", &GenZ2, "GenZ2[NPair]/D");
   OutputPairTree.Branch("RecoE1", &RecoE1, "RecoE1[NPair]/D");
   OutputPairTree.Branch("RecoX1", &RecoX1, "RecoX1[NPair]/D");
   OutputPairTree.Branch("RecoY1", &RecoY1, "RecoY1[NPair]/D");
   OutputPairTree.Branch("RecoZ1", &RecoZ1, "RecoZ1[NPair]/D");
   OutputPairTree.Branch("RecoE2", &RecoE2, "RecoE2[NPair]/D");
   OutputPairTree.Branch("RecoX2", &RecoX2, "RecoX2[NPair]/D");
   OutputPairTree.Branch("RecoY2", &RecoY2, "RecoY2[NPair]/D");
   OutputPairTree.Branch("RecoZ2", &RecoZ2, "RecoZ2[NPair]/D");
   OutputPairTree.Branch("DistanceGen", &DistanceGen, "DistanceGen[NPair]/D");
   OutputPairTree.Branch("DistanceReco", &DistanceReco, "DistanceReco[NPair]/D");
   OutputPairTree.Branch("Distance1", &Distance1, "Distance1[NPair]/D");
   OutputPairTree.Branch("Distance2", &Distance2, "Distance2[NPair]/D");
   OutputPairTree.Branch("E1E2Gen", &E1E2Gen, "E1E2Gen[NPair]/D");
   OutputPairTree.Branch("E1E2Reco", &E1E2Reco, "E1E2Reco[NPair]/D");
   OutputPairTree.Branch("RecoEfficiency1", &RecoEfficiency1, "RecoEfficiency1[NPair]/D");
   OutputPairTree.Branch("RecoEfficiency2", &RecoEfficiency2, "RecoEfficiency2[NPair]/D");

   // variables for the unmatched pair tree
   int NUnmatchedPairGen;
   int NUnmatchedPairReco;
   double GenE1Unmatched[MAXPAIR], GenX1Unmatched[MAXPAIR], GenY1Unmatched[MAXPAIR], GenZ1Unmatched[MAXPAIR];
   double GenE2Unmatched[MAXPAIR], GenX2Unmatched[MAXPAIR], GenY2Unmatched[MAXPAIR], GenZ2Unmatched[MAXPAIR];
   double RecoE1Unmatched[MAXPAIR], RecoX1Unmatched[MAXPAIR], RecoY1Unmatched[MAXPAIR], RecoZ1Unmatched[MAXPAIR];
   double RecoE2Unmatched[MAXPAIR], RecoX2Unmatched[MAXPAIR], RecoY2Unmatched[MAXPAIR], RecoZ2Unmatched[MAXPAIR];
   double DistanceUnmatchedGen[MAXPAIR], DistanceUnmatchedReco[MAXPAIR], DeltaPhiUnmatchedGen[MAXPAIR], DeltaPhiUnmatchedReco[MAXPAIR], DeltaEUnmatchedReco[MAXPAIR], DeltaEUnmatchedGen[MAXPAIR], DeltaThetaUnmatchedGen[MAXPAIR], DeltaThetaUnmatchedReco[MAXPAIR];
   double E1E2GenUnmatched[MAXPAIR], E1E2RecoUnmatched[MAXPAIR];
   OutputUnmatchedTree.Branch("NUnmatchedPairGen", &NUnmatchedPairGen);
   OutputUnmatchedTree.Branch("NUnmatchedPairReco", &NUnmatchedPairReco);
   OutputUnmatchedTree.Branch("GenE1Unmatched", &GenE1Unmatched, "GenE1Unmatched[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("GenE2Unmatched", &GenE2Unmatched, "GenE2Unmatched[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("RecoE1Unmatched", &RecoE1Unmatched, "RecoE1Unmatched[NUnmatchedPairReco]/D");
   OutputUnmatchedTree.Branch("RecoE2Unmatched", &RecoE2Unmatched, "RecoE2Unmatched[NUnmatchedPairReco]/D");
   OutputUnmatchedTree.Branch("DistanceUnmatchedGen", &DistanceUnmatchedGen, "DistanceUnmatchedGen[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("DistanceUnmatchedReco", &DistanceUnmatchedReco, "DistanceUnmatchedReco[NUnmatchedPairReco]/D");
   OutputUnmatchedTree.Branch("DeltaPhiUnmatchedGen", &DeltaPhiUnmatchedGen, "DeltaPhiUnmatchedGen[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("DeltaPhiUnmatchedReco", &DeltaPhiUnmatchedReco, "DeltaPhiUnmatchedReco[NUnmatchedPairReco]/D");
   OutputUnmatchedTree.Branch("DeltaEUnmatchedGen", &DeltaEUnmatchedGen, "DeltaEUnmatchedGen[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("DeltaEUnmatchedReco", &DeltaEUnmatchedReco, "DeltaEUnmatchedReco[NUnmatchedPairReco]/D");
   OutputUnmatchedTree.Branch("DeltaThetaUnmatchedGen", &DeltaThetaUnmatchedGen, "DeltaThetaUnmatchedGen[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("DeltaThetaUnmatchedReco", &DeltaThetaUnmatchedReco, "DeltaThetaUnmatchedReco[NUnmatchedPairReco]/D");
   OutputUnmatchedTree.Branch("E1E2GenUnmatched", &E1E2GenUnmatched, "E1E2GenUnmatched[NUnmatchedPairGen]/D");
   OutputUnmatchedTree.Branch("E1E2RecoUnmatched", &E1E2RecoUnmatched, "E1E2RecoUnmatched[NUnmatchedPairReco]/D");

   ParticleTreeMessenger MGen(InputFile, GenTreeName);
   ParticleTreeMessenger MReco(InputFile, RecoTreeName);

   int EntryCount = MReco.GetEntries() * Fraction;
   int nAcceptedEvents = 0;
   ProgressBar Bar(cout, EntryCount);
   Bar.SetStyle(-1);
   for(int iE = 0; iE < EntryCount; iE++)
   {

      MGen.GetEntry(iE);
      MReco.GetEntry(iE);

      // event selection (not applied for the moment as this is assumed to be applied to the tree)
      /*
      if(MReco.STheta < (7*M_PI)/36 || MReco.STheta >(29*M_PI)/36)continue;
      if(MGen.STheta < (7*M_PI)/36 || MGen.STheta > (29*M_PI)/36)continue;
      if(MGen.passesAll == false) continue;
      if(MReco.passesAll == false) continue;
      if(MReco.PassBaselineCut() == false) continue;
      if(MReco.passesLEP1TwoPC == false) continue;
      */


      // place a cut on the number of charged hadrons if you want to make the mult plot
      // if( MGen.nChargedHadronsHP < 40)continue;


      vector<FourVector> PGen, PReco;
      for(int i = 0; i < MGen.nParticle; i++){
         if(MGen.charge[i] == 0) continue;
         if(MGen.highPurity[i] == false) continue;
         PGen.push_back(MGen.P[i]);
      }

      for(int i = 0; i < MReco.nParticle; i++){
         if(MReco.charge[i] == 0) continue;
         if(MReco.highPurity[i] == false) continue;
         // place cut on the reco energy, not included at gen level
         if(MReco.P[i][0] < 0.2) continue;
         PReco.push_back(MReco.P[i]);
      }


      // now fill the unmatched tree
      int index_counter = 0;
      double TotalE = 91.1876;
      int NunmatchedRecoPairs = 0;
      for(int i = 0; i < PReco.size(); i++){
         for(int j = i+1; j < PReco.size();j++){
            if(i == j) continue; // don't match particles with themselves
            // fill the tree
            FourVector Reco1 = PReco.at(i);
            FourVector Reco2 = PReco.at(j);
            RecoE1Unmatched[index_counter] = Reco1[0];
            RecoE2Unmatched[index_counter] = Reco2[0];
            DistanceUnmatchedReco[index_counter] = GetAngle(Reco1,Reco2);
            DeltaPhiUnmatchedReco[index_counter] = GetDPhi(Reco1, Reco2);
            DeltaEUnmatchedReco[index_counter] = Reco1[0]-Reco2[0];
            DeltaThetaUnmatchedReco[index_counter] = Reco1.GetTheta() - Reco2.GetTheta();
            E1E2RecoUnmatched[index_counter] =  Reco1[0]*Reco2[0];

            // fill the theta distributions
            int BinThetaMeasuredMC = FindBin(GetAngle(Reco1,Reco2), 2 * BinCount, Bins);
            recoUnmatched.Fill(BinThetaMeasuredMC,Reco1[0]*Reco2[0]/(TotalE*TotalE));

            // fille in the z distributions
            double zUnmatched = (1-cos(GetAngle(Reco1, Reco2)))/2;
            int BinZMeasured = FindBin(zUnmatched, 2*BinCount, zBins);
            recoUnmatched_z.Fill(BinZMeasured, Reco1[0]*Reco2[0]/(TotalE*TotalE));

            // fill the energy distributions
            int BinEnergyMeasured = FindBin(Reco1[0]*Reco2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            e1e2RecoUnmatched.Fill(Reco1[0]*Reco2[0]/(TotalE*TotalE));

            NunmatchedRecoPairs++;
            index_counter++;
         }
      }

      NUnmatchedPairReco = NunmatchedRecoPairs;

      int index_gen = 0;
      int NunmatchedGenPairs = 0;
      for(int i = 0; i < PGen.size(); i++){
      for(int j = i+1; j < PGen.size();j++){
            FourVector Gen1 = PGen.at(i);
            FourVector Gen2 = PGen.at(j);
            GenE1Unmatched[index_gen] = Gen1[0];
            GenE2Unmatched[index_gen] = Gen2[0];
            DistanceUnmatchedGen[index_gen] = GetAngle(Gen1,Gen2);
            DeltaPhiUnmatchedGen[index_gen] = GetDPhi(Gen1, Gen2);
            DeltaEUnmatchedGen[index_gen] = Gen1[0]-Gen2[0];
            DeltaThetaUnmatchedGen[index_gen] = Gen1.GetTheta() - Gen2.GetTheta();
            E1E2GenUnmatched[index_gen] =  Gen1[0]*Gen2[0];

            // theta histograms
            int BinThetaGenMC = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
            genUnmatched.Fill(BinThetaGenMC,Gen1[0]*Gen2[0]/(TotalE*TotalE));

            // energy histograms
            int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), BinCount, EnergyBins);
            e1e2GenUnmatched.Fill(Gen1[0]*Gen2[0]/(TotalE*TotalE));

            // z histograms
            double zGenUnmatched = (1-cos(GetAngle(Gen1, Gen2)))/2;
            int BinZGen = FindBin(zGenUnmatched, 2*BinCount, zBins);
            genUnmatched_z.Fill(BinZGen, Gen1[0]*Gen2[0]/(TotalE*TotalE));


            index_gen++;
            NunmatchedGenPairs++;
         }
      }

      NUnmatchedPairGen = NunmatchedGenPairs;


      nAcceptedEvents++;
      OutputUnmatchedTree.Fill();

      // perform the matching
      map<int, int> Matching = MatchJetsHungarian(MatchingMetric, PGen, PReco);
      int Count = 0;
      NParticle = Matching.size();
      eventID = MGen.EventNo;
      for(auto iter : Matching)
      {
         FourVector Gen = iter.first >= 0 ? PGen[iter.first] : FourVector(-1, 0, 0, 0);
         FourVector Reco = iter.second >= 0 ? PReco[iter.second] : FourVector(-1, 0, 0, 0);
         // if(GetAngle(Gen, Reco) > 1.0) {
         //    Gen = FourVector(-1, 0, 0, 0);
         //    Reco = FourVector(-1, 0, 0, 0);
         // }
         GenE[Count] = Gen[0];
         GenX[Count] = Gen[1];
         GenY[Count] = Gen[2];
         GenZ[Count] = Gen[3];
         GenEta[Count] = Gen.GetEta();
         GenRapidity[Count] = Gen.GetRapidity();
         GenPhi[Count] = Gen.GetPhi();
         GenTheta[Count] = Gen.GetTheta();
         RecoE[Count] = Reco[0];
         RecoX[Count] = Reco[1];
         RecoY[Count] = Reco[2];
         RecoZ[Count] = Reco[3];
         RecoRapidity[Count] = Reco.GetRapidity();
         RecoEta[Count] = Reco.GetEta();
         RecoPhi[Count] = Reco.GetPhi();
         RecoTheta[Count] = Reco.GetTheta();
         Distance[Count] = GetAngle(Gen, Reco);
         DeltaPhi[Count] = GetDPhi(Gen,Reco);
         DeltaE[Count] = Gen[0] - Reco[0];
         DeltaTheta[Count] = Gen.GetTheta() - Reco.GetTheta();
         double Efficiency;
         if (Reco.GetPT() <  0.2) Efficiency = 1;
         else Efficiency = efficiencyCorrector.efficiency(Reco.GetTheta(), Reco.GetPhi(), Reco.GetPT(), MReco.nChargedHadronsHP);
         RecoEfficiency[Count] = 1/Efficiency;
         Count = Count + 1;
      }

      OutputTree.Fill(); // fill the tree

      // now fill the tree for the matched pairs
      NPair = 0;
      int NMatchedPairs = 0;
      for(int i = 0; i < NParticle; i++)
      {
         for(int j = i + 1; j < NParticle; j++)
         {
            GenE1[NPair] = GenE[i];
            GenX1[NPair] = GenX[i];
            GenY1[NPair] = GenY[i];
            GenZ1[NPair] = GenZ[i];
            GenE2[NPair] = GenE[j];
            GenX2[NPair] = GenX[j];
            GenY2[NPair] = GenY[j];
            GenZ2[NPair] = GenZ[j];
            RecoE1[NPair] = RecoE[i];
            RecoX1[NPair] = RecoX[i];
            RecoY1[NPair] = RecoY[i];
            RecoZ1[NPair] = RecoZ[i];
            RecoE2[NPair] = RecoE[j];
            RecoX2[NPair] = RecoX[j];
            RecoY2[NPair] = RecoY[j];
            RecoZ2[NPair] = RecoZ[j];
            RecoEfficiency1[NPair] = RecoEfficiency[i];
            RecoEfficiency2[NPair] = RecoEfficiency[j];

            Distance1[NPair] = Distance[i];
            Distance2[NPair] = Distance[j];

            FourVector Gen1(GenE[i], GenX[i], GenY[i], GenZ[i]);
            FourVector Gen2(GenE[j], GenX[j], GenY[j], GenZ[j]);
            FourVector Reco1(RecoE[i], RecoX[i], RecoY[i], RecoZ[i]);
            FourVector Reco2(RecoE[j], RecoX[j], RecoY[j], RecoZ[j]);

            DistanceGen[NPair] = GetAngle(Gen1, Gen2);
            DistanceReco[NPair] = GetAngle(Reco1, Reco2);

            E1E2Gen[NPair] = (GenE[i]*GenE[j])/(TotalE*TotalE);
            E1E2Reco[NPair] = (RecoE[i]*RecoE[j])/(TotalE*TotalE);

            if(RecoE[i] > 0 && RecoE[j] > 0 && GenE[i] > 0 && GenE[j] > 0){
               // theta histograms
               int BinThetaMeasuredMC = FindBin(GetAngle(Reco1,Reco2), 2 * BinCount, Bins);
               recoMatched.Fill(BinThetaMeasuredMC,Reco1[0]*Reco2[0]/(TotalE*TotalE));
               int BinThetaGenMC = FindBin(GetAngle(Gen1,Gen2), 2 * BinCount, Bins);
               genMatched.Fill(BinThetaGenMC,Gen1[0]*Gen2[0]/(TotalE*TotalE));

               // energy histograms
               int BinEnergyMeasured = FindBin(Reco1[0]*Reco2[0]/(TotalE*TotalE), BinCount, EnergyBins);
               e1e2RecoMatched.Fill(Reco1[0]*Reco2[0]/(TotalE*TotalE));
               int BinEnergyGen = FindBin(Gen1[0]*Gen2[0]/(TotalE*TotalE), BinCount, EnergyBins);
               e1e2GenMatched.Fill(Gen1[0]*Gen2[0]/(TotalE*TotalE));

               // z histograms
               double zRecoMatched = (1-cos(GetAngle(Reco1, Reco2)))/2;
               int BinZMeasured = FindBin(zRecoMatched, 2*BinCount, zBins);
               recoMatched_z.Fill(BinZMeasured, Reco1[0]*Reco2[0]/(TotalE*TotalE));

               double zGenMatched = (1-cos(GetAngle(Gen1, Gen2)))/2;
               int BinZMC = FindBin(zGenMatched, 2*BinCount, zBins);
               genMatched_z.Fill(BinZMC, Gen1[0]*Gen2[0]/(TotalE*TotalE));

               // increment the number of matched pairs
               NMatchedPairs++;
            }

            NPair = NPair + 1;
         }
      }

      OutputPairTree.Fill();
   }
   std::cout << "Number of accepted events is " << nAcceptedEvents << std::endl;
   Bar.Update(EntryCount);
   Bar.Print();
   Bar.PrintLine();

   OutputFile.cd();

   // write the trees to the output file
   OutputTree.Write();
   OutputPairTree.Write();
   OutputUnmatchedTree.Write();

   // -------------------------------------------------------------------
   // histograms for the matching efficiency and fake fraction
   // -------------------------------------------------------------------

   //function of theta
   recoUnmatched.Write();
   recoMatched.Write();
   genUnmatched.Write();
   genMatched.Write();

   // function of z
   recoUnmatched_z.Write();
   recoMatched_z.Write();
   genUnmatched_z.Write();
   genMatched_z.Write();

   // function of energy
   e1e2GenMatched.Write();
   e1e2RecoMatched.Write();
   e1e2GenUnmatched.Write();
   e1e2RecoUnmatched.Write();
   // -------------------------------------------------------------------

   // write the output files
   OutputFile.Close();
   MatchingPerformance(OutputFileName, rstDirName, MGen, MReco);
   InputFile.Close();


   return 0;
}

// metric angle used as the first guess of the matching metric
double MetricAngle(FourVector A, FourVector B)
{
   double Angle = GetAngle(A, B);
   return Angle;
}

// core function to calculate the chisquare(theta,phi,E) based on different assignments of resolutions
double MatchingMetricCore( double deltaTheta, double deltaPhi, double deltaE,
                           double meanE,     // the average btw the GenE and RecoE,
                                             // this would be used to model the energy-dependency in the resolution assignment
                           int schemeChoice, // 1: 10% energy as resolution, flat theta and phi resolutions
                                             // 2: resolutions are cast according to ALEPH tracker performance (energy-dependent description)
                                             // 3: scale btw the importance of angular-match versus energy-match ( 5x)
                                             // 4: scale btw the importance of angular-match versus energy-match (15x)
                           double& chiTheta, double& chiPhi, double& chiE)
{
   if (schemeChoice==1)
   {
      double phiRes  = 0.002; // somewhere in the ballpark of 0.002-0.005
      double Eres    = 0.1*meanE; // use 10% energy resolution //0.85/sqrt(meanE); // take energy resolution based on the mean energy
      double AngleRes= (0.01*0.01); // take angular resolution based on the detector resolution
      chiTheta       = deltaTheta/AngleRes;
      chiPhi         = deltaPhi  /phiRes;
      chiE           = deltaE    /Eres;
   }
   else if (schemeChoice==2 || schemeChoice==3 || schemeChoice==4)
   {
      double scaleFactorTheta = 2.8; // sigma(rz)   = 28 µm
      double scaleFactorPhi   = 2.3; // sigma(rphi) = 23 µm
      double scaleFactorE     = (schemeChoice==2)? 1:
                                (schemeChoice==3)? 5: 15;
      double sigmaDelta = 25e-6 + 95e-6 / meanE;
      double sigmaTheta = sigmaDelta / 6e-2     // inner vertex detector radius
                        * scaleFactorTheta;     // considering the projection (average) of minimal distance to the r-z
      double sigmaPhi   = sigmaDelta / 6e-2     // inner vertex detector radius
                        * scaleFactorPhi;       // considering the projection (average) of minimal distance to the r-z
      double sigmaE = TMath::Sqrt( (6e-4*meanE)*(6e-4*meanE) + 0.005 * 0.005 ) * meanE
                        * scaleFactorE;
      chiTheta       = deltaTheta/sigmaTheta;
      chiPhi         = deltaPhi  /sigmaPhi;
      chiE           = deltaE    /sigmaE;
   }
   else
   {
      printf("[Error] schemeChoice %d is not supported in MatchingMetricCore. Please choose a number btw 1-4. Exiting...\n", schemeChoice);
      exit(1);
   }

   double chi2metric = chiTheta*chiTheta + chiPhi*chiPhi + chiE*chiE;
   return (chi2metric);
}

double MatchingMetric(FourVector A, FourVector B){
   double Angle = GetAngle(A,B);
   double dPhi = GetDPhi(A,B);
   double Ediff = (A[0] - B[0]);
   double meanE = (A[0] + B[0])/2;

   double _chiTheta, _chiPhi, _chiE;
   double chi2metric = MatchingMetricCore( Angle, dPhi, Ediff,
                                           meanE, MatchingSchemeChoice,
                                           _chiTheta, _chiPhi, _chiE);
   return chi2metric;
}

// function for filling the histograms
int FindBin(double Value, int NBins, double Bins[])
{
   for(int i = 0; i < NBins; i++)
      if(Value < Bins[i])
         return i - 1;
   return NBins;
}

// can we put this in the CommonCode/ ?
void FillChain(TChain &chain, const vector<string> &files) {
  for (auto file : files) {
    chain.Add(file.c_str());
  }
}

// Matching criteria performance check
void MatchingPerformance(string& matchedRstRootName, string& rstDirName,
                         ParticleTreeMessenger& MGen,
                         ParticleTreeMessenger& MReco)
{
   //------------------------------------
   // Single-track level
   //------------------------------------
   TChain particleChain("MatchedTree");
   FillChain(particleChain, {matchedRstRootName});
   TTreeReader particleReader(&particleChain);

   TTreeReaderValue<int> nParticles(   particleReader, "NParticle");
   TTreeReaderArray<double> GenE(      particleReader, "GenE");
   TTreeReaderArray<double> GenX(      particleReader, "GenX");
   TTreeReaderArray<double> GenY(      particleReader, "GenY");
   TTreeReaderArray<double> GenZ(      particleReader, "GenZ");
   TTreeReaderArray<double> GenTheta(  particleReader, "GenTheta");
   TTreeReaderArray<double> GenPhi(    particleReader, "GenPhi");
   TTreeReaderArray<double> RecoE(     particleReader, "RecoE");
   TTreeReaderArray<double> RecoX(     particleReader, "RecoX");
   TTreeReaderArray<double> RecoY(     particleReader, "RecoY");
   TTreeReaderArray<double> RecoZ(     particleReader, "RecoZ");
   TTreeReaderArray<double> RecoTheta( particleReader, "RecoTheta");
   TTreeReaderArray<double> RecoPhi(   particleReader, "RecoPhi");
   TTreeReaderArray<double> Distance(  particleReader, "Distance");

   //------------------------------------
   // define the histograms
   //------------------------------------
   string PerformanceFileName(matchedRstRootName);
   PerformanceFileName.replace(PerformanceFileName.find(".root"), 5, "_performance.root");
   TFile PerformanceFile(PerformanceFileName.c_str(), "RECREATE");

   PerformanceFile.cd();

   TH1D trkChi2Tot    ("trkChi2Tot",    "trkChi2Tot",  50, -5, 5);
   TH1D trkChiTheta   ("trkChiTheta",  "trkChiTheta",50, -5, 5);
   TH1D trkChiPhi     ("trkChiPhi",    "trkChiPhi",  50, -5, 5);
   TH1D trkChiE       ("trkChiE",      "trkChiE",    50, -5, 5);
   TH1D trkDeltaTheta ("trkDeltaTheta",  "trkDeltaTheta",50, -4, 4);
   TH1D trkDeltaPhi   ("trkDeltaPhi",    "trkDeltaPhi",  50, -4, 4);
   TH1D trkDeltaE     ("trkDeltaE",      "trkDeltaE",    50, -10, 10);
   TH1D trkDeltaP     ("trkDeltaP",      "trkDeltaP",    50, -10, 10);
   TH1D trkDeltaTheta_ge5 ("trkDeltaTheta_ge5",  "trkDeltaTheta (E>5GeV)",50, -4, 4);
   TH1D trkDeltaPhi_ge5   ("trkDeltaPhi_ge5",    "trkDeltaPhi (E>5GeV)",  50, -4, 4);
   TH1D trkDeltaE_ge5     ("trkDeltaE_ge5",      "trkDeltaE (E>5GeV)",    50, -10, 10);
   TH1D trkDeltaP_ge5     ("trkDeltaP_ge5",      "trkDeltaP (E>5GeV)",    50, -10, 10);
   TH1D trkDeltaTheta_2to5 ("trkDeltaTheta_2to5",  "trkDeltaTheta (2<E<5GeV)",50, -4, 4);
   TH1D trkDeltaPhi_2to5   ("trkDeltaPhi_2to5",    "trkDeltaPhi (2<E<5GeV)",  50, -4, 4);
   TH1D trkDeltaE_2to5     ("trkDeltaE_2to5",      "trkDeltaE (2<E<5GeV)",    50, -10, 10);
   TH1D trkDeltaP_2to5     ("trkDeltaP_2to5",      "trkDeltaP (2<E<5GeV)",    50, -10, 10);
   TH1D trkDeltaTheta_1to2 ("trkDeltaTheta_1to2",  "trkDeltaTheta (1<E<2GeV)",50, -4, 4);
   TH1D trkDeltaPhi_1to2   ("trkDeltaPhi_1to2",    "trkDeltaPhi (1<E<2GeV)",  50, -4, 4);
   TH1D trkDeltaE_1to2     ("trkDeltaE_1to2",      "trkDeltaE (1<E<2GeV)",    50, -10, 10);
   TH1D trkDeltaP_1to2     ("trkDeltaP_1to2",      "trkDeltaP (1<E<2GeV)",    50, -10, 10);
   TH1D trkDeltaTheta_le1 ("trkDeltaTheta_le1",  "trkDeltaTheta (E<1GeV)",50, -4, 4);
   TH1D trkDeltaPhi_le1   ("trkDeltaPhi_le1",    "trkDeltaPhi (E<1GeV)",  50, -4, 4);
   TH1D trkDeltaE_le1     ("trkDeltaE_le1",      "trkDeltaE (E<1GeV)",    50, -10, 10);
   TH1D trkDeltaP_le1     ("trkDeltaP_le1",      "trkDeltaP (E<1GeV)",    50, -10, 10);

   int nGenTracks       = 0;
   int nRecoTracks      = 0;
   int nMatchedTracks   = 0;
   int nUnmatchedTracks = 0;
   int nMatchedTracks_1p0   = 0;
   int nUnmatchedTracks_1p0 = 0;
   int nMatchedTracks_0p4   = 0;
   int nUnmatchedTracks_0p4 = 0;
   int nMatchedTracks_0p2   = 0;
   int nUnmatchedTracks_0p2 = 0;
   int nMatchedTracks_0p1   = 0;
   int nUnmatchedTracks_0p1 = 0;
   Long64_t totEvents = particleReader.GetEntries(true);
   for (Long64_t iEvent = 0; iEvent < totEvents; iEvent++)
   {
      particleReader.Next();

      int nMatchedTracksInOneEvent = 0;
      int nMatchedTracksInOneEvent_1p0 = 0;
      int nMatchedTracksInOneEvent_0p4 = 0;
      int nMatchedTracksInOneEvent_0p2 = 0;
      int nMatchedTracksInOneEvent_0p1 = 0;
      for (int iParticle = 0; iParticle < *nParticles; iParticle++)
      {
         if (RecoE[iParticle]<0 || GenE[iParticle]<0) continue; // remove non-matched pairs

         double deltaTheta = RecoTheta[iParticle]-GenTheta[iParticle];
         double deltaPhi   = RecoPhi[iParticle]-GenPhi[iParticle];
         double deltaE     = RecoE[iParticle]-GenE[iParticle];
         double deltaP     = TMath::Sqrt(RecoX[iParticle]*RecoX[iParticle]+RecoY[iParticle]*RecoY[iParticle]+RecoZ[iParticle]*RecoZ[iParticle])-
                             TMath::Sqrt(GenX[iParticle]*GenX[iParticle]+GenY[iParticle]*GenY[iParticle]+GenZ[iParticle]*GenZ[iParticle]);

         trkChi2Tot.Fill(Distance[iParticle]);
         trkDeltaTheta.Fill(deltaTheta);
         trkDeltaPhi.Fill(deltaPhi);
         trkDeltaE.Fill(deltaE);
         trkDeltaP.Fill(deltaP);

         if (RecoE[iParticle] >= 5)
         {
            trkDeltaTheta_ge5.Fill(deltaTheta);
            trkDeltaPhi_ge5.Fill(deltaPhi);
            trkDeltaE_ge5.Fill(deltaE);
            trkDeltaP_ge5.Fill(deltaP);
         }
         else if (2 <= RecoE[iParticle] && RecoE[iParticle] < 5)
         {
            trkDeltaTheta_2to5.Fill(deltaTheta);
            trkDeltaPhi_2to5.Fill(deltaPhi);
            trkDeltaE_2to5.Fill(deltaE);
            trkDeltaP_2to5.Fill(deltaP);
         }
         else if (1 <= RecoE[iParticle] && RecoE[iParticle] < 2)
         {
            trkDeltaTheta_1to2.Fill(deltaTheta);
            trkDeltaPhi_1to2.Fill(deltaPhi);
            trkDeltaE_1to2.Fill(deltaE);
            trkDeltaP_1to2.Fill(deltaP);
         }
         else // recoE < 1
         {
            trkDeltaTheta_le1.Fill(deltaTheta);
            trkDeltaPhi_le1.Fill(deltaPhi);
            trkDeltaE_le1.Fill(deltaE);
            trkDeltaP_le1.Fill(deltaP);
         }

         double meanE = (RecoE[iParticle] + GenE[iParticle])/2;

         double chiTheta, chiPhi, chiE;
         double chi2metric = MatchingMetricCore( deltaTheta, deltaPhi, deltaE,
                                                 meanE, MatchingSchemeChoice,
                                                 chiTheta, chiPhi, chiE);

         trkChiTheta.Fill(chiTheta);
         trkChiPhi.Fill(chiPhi);
         trkChiE.Fill(chiE);

         nMatchedTracksInOneEvent++;
         if (Distance[iParticle]<=1.0) nMatchedTracksInOneEvent_1p0++;
         if (Distance[iParticle]<=0.4) nMatchedTracksInOneEvent_0p4++;
         if (Distance[iParticle]<=0.2) nMatchedTracksInOneEvent_0p2++;
         if (Distance[iParticle]<=0.1) nMatchedTracksInOneEvent_0p1++;
      }

      MGen.GetEntry(iEvent);
      MReco.GetEntry(iEvent);

      // printf("%d %d %d %d %d %d\n", MGen.nParticle, MGen.nChargedHadronsHP, MReco.nParticle, MReco.nChargedHadronsHP, (*nParticles), nMatchedTracksInOneEvent);
      nGenTracks       += MGen.nChargedHadronsHP;
      nRecoTracks      += MReco.nChargedHadronsHP;
      nMatchedTracks   += nMatchedTracksInOneEvent;
      nUnmatchedTracks += (MReco.nChargedHadronsHP-nMatchedTracksInOneEvent);
      nMatchedTracks_1p0   += nMatchedTracksInOneEvent_1p0;
      nUnmatchedTracks_1p0 += (MReco.nChargedHadronsHP-nMatchedTracksInOneEvent_1p0);
      nMatchedTracks_0p4   += nMatchedTracksInOneEvent_0p4;
      nUnmatchedTracks_0p4 += (MReco.nChargedHadronsHP-nMatchedTracksInOneEvent_0p4);
      nMatchedTracks_0p2   += nMatchedTracksInOneEvent_0p2;
      nUnmatchedTracks_0p2 += (MReco.nChargedHadronsHP-nMatchedTracksInOneEvent_0p2);
      nMatchedTracks_0p1   += nMatchedTracksInOneEvent_0p1;
      nUnmatchedTracks_0p1 += (MReco.nChargedHadronsHP-nMatchedTracksInOneEvent_0p1);

   } // end loop over the number of events

   printf("Summary >>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
   printf("Efficiency = %d/%d = %.3f\n", nMatchedTracks, nGenTracks, nMatchedTracks/((double) nGenTracks));
   printf("Fake rate  = %d/%d = %.3f\n", nUnmatchedTracks, nRecoTracks, nUnmatchedTracks/((double) nRecoTracks));
   printf("Efficiency = %d/%d = %.3f | cutoff = 1.0\n", nMatchedTracks_1p0, nGenTracks, nMatchedTracks_1p0/((double) nGenTracks));
   printf("Fake rate  = %d/%d = %.3f | cutoff = 1.0\n", nUnmatchedTracks_1p0, nRecoTracks, nUnmatchedTracks_1p0/((double) nRecoTracks));
   printf("Efficiency = %d/%d = %.3f | cutoff = 0.4\n", nMatchedTracks_0p4, nGenTracks, nMatchedTracks_0p4/((double) nGenTracks));
   printf("Fake rate  = %d/%d = %.3f | cutoff = 0.4\n", nUnmatchedTracks_0p4, nRecoTracks, nUnmatchedTracks_0p4/((double) nRecoTracks));
   printf("Efficiency = %d/%d = %.3f | cutoff = 0.2\n", nMatchedTracks_0p2, nGenTracks, nMatchedTracks_0p2/((double) nGenTracks));
   printf("Fake rate  = %d/%d = %.3f | cutoff = 0.2\n", nUnmatchedTracks_0p2, nRecoTracks, nUnmatchedTracks_0p2/((double) nRecoTracks));
   printf("Efficiency = %d/%d = %.3f | cutoff = 0.1\n", nMatchedTracks_0p1, nGenTracks, nMatchedTracks_0p1/((double) nGenTracks));
   printf("Fake rate  = %d/%d = %.3f | cutoff = 0.1\n", nUnmatchedTracks_0p1, nRecoTracks, nUnmatchedTracks_0p1/((double) nRecoTracks));
   printf("End Summary <<<<<<<<<<<<<<<<<<<<<<<<<\n");

   trkChi2Tot.Write();
   trkChiTheta.Write();
   trkChiPhi.Write();
   trkChiE.Write();
   trkDeltaTheta.Write();
   trkDeltaPhi.Write();
   trkDeltaE.Write();
   trkDeltaP.Write();
   trkDeltaTheta_ge5.Write();
   trkDeltaPhi_ge5.Write();
   trkDeltaE_ge5.Write();
   trkDeltaP_ge5.Write();
   trkDeltaTheta_2to5.Write();
   trkDeltaPhi_2to5.Write();
   trkDeltaE_2to5.Write();
   trkDeltaP_2to5.Write();
   trkDeltaTheta_1to2.Write();
   trkDeltaPhi_1to2.Write();
   trkDeltaE_1to2.Write();
   trkDeltaP_1to2.Write();
   trkDeltaTheta_le1.Write();
   trkDeltaPhi_le1.Write();
   trkDeltaE_le1.Write();
   trkDeltaP_le1.Write();

   auto plotResolution = [](TH1D& h, string rstDirName,
                            bool doFit=true, bool logy=false)
   {
      TCanvas *c = new TCanvas("c", "", 600, 600);
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.05);

      TF1 f("f","[p0]/sqrt(2*pi)/[p2]*exp(-0.5*((x-[p1])/[p2])*((x-[p1])/[p2]))+[p3]/sqrt(2*pi)/[p4]*exp(-0.5*((x-[p1])/[p4])*((x-[p1])/[p4]))",-5,5);
      double hIntegral = h.Integral(h.GetXaxis()->FindBin(-5+0.001),
                                    h.GetXaxis()->FindBin( 5-0.001)) * h.GetXaxis()->GetBinWidth(1);
      f.SetParameters(hIntegral*0.7,
                      0,
                      1,
                      hIntegral*0.3,
                      2.5);
      f.SetParLimits(0, hIntegral*0.5, hIntegral*1.);
      f.SetParLimits(2, 0.01, 1.5);
      f.SetParLimits(3, 0, hIntegral*0.5);
      f.SetParLimits(4, 1.5, 5);

      h.GetYaxis()->SetRangeUser((logy)? 7e-1: 0, h.GetMaximum()*1.6);
      h.Draw();

      if (logy) gPad->SetLogy();

      if (doFit)
      {
         h.Fit("f","RSLM");

         TLatex latex;
         latex.SetNDC();  // Use normalized coordinates [0,1] for positioning
         latex.SetTextSize(0.04);

         // Write the fit parameters on the canvas
         latex.DrawLatex(0.15, 0.85, Form("Ampl. 1st Gaus = %.2f #pm %.2f", f.GetParameter(0), f.GetParError(0)));
         latex.DrawLatex(0.15, 0.80, Form("Mean Double Gaus = %.2f #pm %.2f", f.GetParameter(1), f.GetParError(1)));
         latex.DrawLatex(0.15, 0.75, Form("Width 1st Gaus = %.2f #pm %.2f", f.GetParameter(2), f.GetParError(2)));
         latex.DrawLatex(0.15, 0.70, Form("Ampl. 2nd Gaus = %.2f #pm %.2f", f.GetParameter(3), f.GetParError(3)));
         latex.DrawLatex(0.15, 0.65, Form("Width 2nd Gaus = %.2f #pm %.2f", f.GetParameter(4), f.GetParError(4)));
         latex.DrawLatex(0.15, 0.60, Form("(Total normalization = %.2f)", hIntegral));

         TF1 f_sub("f_sub","[p0]/sqrt(2*pi)/[p2]*exp(-0.5*((x-[p1])/[p2])*((x-[p1])/[p2]))+[p3]/sqrt(2*pi)/[p4]*exp(-0.5*((x-[p1])/[p4])*((x-[p1])/[p4]))",-5,5);
         f_sub.SetParameter(0, f.GetParameter(0));
         f_sub.SetParameter(1, f.GetParameter(1));
         f_sub.SetParameter(2, f.GetParameter(2));
         f_sub.SetParameter(3, 0);
         f_sub.SetParameter(4, f.GetParameter(4));
         f_sub.SetLineStyle(2);
         f_sub.SetLineColor(kRed);
         f_sub.DrawClone("same");

         f_sub.SetParameter(0, 0);
         f_sub.SetParameter(1, f.GetParameter(1));
         f_sub.SetParameter(2, f.GetParameter(2));
         f_sub.SetParameter(3, f.GetParameter(3));
         f_sub.SetParameter(4, f.GetParameter(4));
         f_sub.SetLineStyle(3);
         f_sub.SetLineColor(kRed);
         f_sub.DrawClone("same");
      }

      system(Form("mkdir -p %s/plot/", rstDirName.c_str()));
      c->SaveAs(Form("%s/plot/%s.pdf", rstDirName.c_str(), h.GetName()));
      // system(Form("dropbox_uploader.sh upload %s/plot/%s.pdf /tmp/", rstDirName.c_str(), h.GetName()));
      delete c;
   };

   plotResolution(trkChiTheta, rstDirName);
   plotResolution(trkChiPhi, rstDirName);
   plotResolution(trkChiE, rstDirName);

   plotResolution(trkDeltaTheta, rstDirName, false, true);
   plotResolution(trkDeltaPhi, rstDirName, false, true);
   plotResolution(trkDeltaE, rstDirName, false, true);
   plotResolution(trkDeltaP, rstDirName, false, true);

   plotResolution(trkDeltaTheta_ge5, rstDirName, false, true);
   plotResolution(trkDeltaPhi_ge5, rstDirName, false, true);
   plotResolution(trkDeltaE_ge5, rstDirName, false, true);
   plotResolution(trkDeltaP_ge5, rstDirName, false, true);

   plotResolution(trkDeltaTheta_2to5, rstDirName, false, true);
   plotResolution(trkDeltaPhi_2to5, rstDirName, false, true);
   plotResolution(trkDeltaE_2to5, rstDirName, false, true);
   plotResolution(trkDeltaP_2to5, rstDirName, false, true);

   plotResolution(trkDeltaTheta_1to2, rstDirName, false, true);
   plotResolution(trkDeltaPhi_1to2, rstDirName, false, true);
   plotResolution(trkDeltaE_1to2, rstDirName, false, true);
   plotResolution(trkDeltaP_1to2, rstDirName, false, true);

   plotResolution(trkDeltaTheta_le1, rstDirName, false, true);
   plotResolution(trkDeltaPhi_le1, rstDirName, false, true);
   plotResolution(trkDeltaE_le1, rstDirName, false, true);
   plotResolution(trkDeltaP_le1, rstDirName, false, true);

   PerformanceFile.Close();
}
