#include <algorithm>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "Messenger.h"

JetTreeMessenger::JetTreeMessenger()
{
   Tree = nullptr;
}

JetTreeMessenger::JetTreeMessenger(TFile &file, std::string name)
{
   TTree *tree = (TTree *)file.Get(name.c_str());
   Initialize(tree);
}

JetTreeMessenger::JetTreeMessenger(TFile *file, std::string name)
{
   if(file == nullptr)
   {
      Tree = nullptr;
      return;
   }
   TTree *tree = (TTree *)file->Get(name.c_str());
   Initialize(tree);
}

JetTreeMessenger::JetTreeMessenger(TTree *tree)
{
   Initialize(tree);
}

bool JetTreeMessenger::Initialize(TTree *tree)
{
   Tree = tree;
   return Initialize();
}

bool JetTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   Tree->SetBranchAddress("nref",                     &nref);
   Tree->SetBranchAddress("jtpt",                     &jtpt);
   Tree->SetBranchAddress("jteta",                    &jteta);
   Tree->SetBranchAddress("jtphi",                    &jtphi);
   Tree->SetBranchAddress("jtm",                      &jtm);
   Tree->SetBranchAddress("jtN",                      &jtN);
   Tree->SetBranchAddress("jtNPW",                    &jtNPW);
   Tree->SetBranchAddress("jtptFracPW",               &jtptFracPW);
   Tree->SetBranchAddress("zgJtPt_Beta0p00ZCut0p10",  &zgJtPt_Beta0p00ZCut0p10);
   Tree->SetBranchAddress("zgJtPhi_Beta0p00ZCut0p10", &zgJtPhi_Beta0p00ZCut0p10);
   Tree->SetBranchAddress("zgJtEta_Beta0p00ZCut0p10", &zgJtEta_Beta0p00ZCut0p10);
   Tree->SetBranchAddress("zg_Beta0p00ZCut0p10",      &zg_Beta0p00ZCut0p10);
   Tree->SetBranchAddress("rg_Beta0p00ZCut0p10",      &rg_Beta0p00ZCut0p10);

   return true;
}

bool JetTreeMessenger::GetEntry(int iEntry)
{
   Jet.clear();

   if(Tree == nullptr)
      return false;
   if(iEntry < 0)
      return false;
   if(iEntry >= GetEntries())
      return false;

   Tree->GetEntry(iEntry);

   Jet.resize(nref);
   for(int i = 0; i < nref; i++)
      Jet[i].SetPtEtaPhiMass(jtpt[i], jteta[i], jtphi[i], jtm[i]);

   return true;
}

int JetTreeMessenger::GetEntries()
{
   if(Tree == nullptr)
      return 0;
   return Tree->GetEntries();
}

ParticleTreeMessenger::ParticleTreeMessenger()
{
   Tree = nullptr;
}

ParticleTreeMessenger::ParticleTreeMessenger(TFile &file, std::string name)
{
   TTree *tree = (TTree *)file.Get(name.c_str());
   Initialize(tree);
}

ParticleTreeMessenger::ParticleTreeMessenger(TFile *file, std::string name)
{
   if(file == nullptr)
   {
      Tree = nullptr;
      return;
   }
   TTree *tree = (TTree *)file->Get(name.c_str());
   Initialize(tree);
}

ParticleTreeMessenger::ParticleTreeMessenger(TTree *tree)
{
   Initialize(tree);
}

bool ParticleTreeMessenger::Initialize(TTree *tree)
{
   Tree = tree;
   return Initialize();
}

bool ParticleTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   if(Tree->GetBranch("EventNo") != nullptr)
      Tree->SetBranchAddress("EventNo", &EventNo);
   if(Tree->GetBranch("RunNo") != nullptr)
      Tree->SetBranchAddress("RunNo", &RunNo);
   if(Tree->GetBranch("year") != nullptr)
      Tree->SetBranchAddress("year", &year);
   if(Tree->GetBranch("subDir") != nullptr)
      Tree->SetBranchAddress("subDir", &subDir);
   if(Tree->GetBranch("process") != nullptr)
      Tree->SetBranchAddress("process", &process);
   if(Tree->GetBranch("isMC") != nullptr)
      Tree->SetBranchAddress("isMC", &isMC);
   if(Tree->GetBranch("uniqueID") != nullptr)
      Tree->SetBranchAddress("uniqueID", &uniqueID);
   if(Tree->GetBranch("Energy") != nullptr)
      Tree->SetBranchAddress("Energy", &Energy);
   if(Tree->GetBranch("bFlag") != nullptr)
      Tree->SetBranchAddress("bFlag", &bFlag);
   if(Tree->GetBranch("particleWeight") != nullptr)
      Tree->SetBranchAddress("particleWeight", &particleWeight);
   if(Tree->GetBranch("bx") != nullptr)
      Tree->SetBranchAddress("bx", &bx);
   if(Tree->GetBranch("by") != nullptr)
      Tree->SetBranchAddress("by", &by);
   if(Tree->GetBranch("ebx") != nullptr)
      Tree->SetBranchAddress("ebx", &ebx);
   if(Tree->GetBranch("eby") != nullptr)
      Tree->SetBranchAddress("eby", &eby);
   if(Tree->GetBranch("nParticle") != nullptr)
      Tree->SetBranchAddress("nParticle", &nParticle);
   if(Tree->GetBranch("px") != nullptr)
      Tree->SetBranchAddress("px", &px);
   if(Tree->GetBranch("py") != nullptr)
      Tree->SetBranchAddress("py", &py);
   if(Tree->GetBranch("pz") != nullptr)
      Tree->SetBranchAddress("pz", &pz);
   if(Tree->GetBranch("pt") != nullptr)
      Tree->SetBranchAddress("pt", &pt);
   if(Tree->GetBranch("pmag") != nullptr)
      Tree->SetBranchAddress("pmag", &pmag);
   if(Tree->GetBranch("rap") != nullptr)
      Tree->SetBranchAddress("rap", &rap);
   if(Tree->GetBranch("eta") != nullptr)
      Tree->SetBranchAddress("eta", &eta);
   if(Tree->GetBranch("theta") != nullptr)
      Tree->SetBranchAddress("theta", &theta);
   if(Tree->GetBranch("phi") != nullptr)
      Tree->SetBranchAddress("phi", &phi);
   if(Tree->GetBranch("mass") != nullptr)
      Tree->SetBranchAddress("mass", &mass);
   if(Tree->GetBranch("charge") != nullptr)
      Tree->SetBranchAddress("charge", &charge);
   if(Tree->GetBranch("isCharged") != nullptr)
      Tree->SetBranchAddress("isCharged", &isCharged);
   if(Tree->GetBranch("pwflag") != nullptr)
      Tree->SetBranchAddress("pwflag", &pwflag);
   if(Tree->GetBranch("pid") != nullptr)
      Tree->SetBranchAddress("pid", &pid);
   if(Tree->GetBranch("d0") != nullptr)
      Tree->SetBranchAddress("d0", &d0);
   if(Tree->GetBranch("z0") != nullptr)
      Tree->SetBranchAddress("z0", &z0);
   if(Tree->GetBranch("highPurity") != nullptr)
      Tree->SetBranchAddress("highPurity", &highPurity);
   if(Tree->GetBranch("ntpc") != nullptr)
      Tree->SetBranchAddress("ntpc", &ntpc);
   if(Tree->GetBranch("nitc") != nullptr)
      Tree->SetBranchAddress("nitc", &nitc);
   if(Tree->GetBranch("nvdet") != nullptr)
      Tree->SetBranchAddress("nvdet", &nvdet);
   if(Tree->GetBranch("vx") != nullptr)
      Tree->SetBranchAddress("vx", &vx);
   if(Tree->GetBranch("vy") != nullptr)
      Tree->SetBranchAddress("vy", &vy);
   if(Tree->GetBranch("vz") != nullptr)
      Tree->SetBranchAddress("vz", &vz);
   if(Tree->GetBranch("weight") != nullptr)
      Tree->SetBranchAddress("weight", &weight);
   if(Tree->GetBranch("pt_wrtThr") != nullptr)
      Tree->SetBranchAddress("pt_wrtThr", &pt_wrtThr);
   if(Tree->GetBranch("eta_wrtThr") != nullptr)
      Tree->SetBranchAddress("eta_wrtThr", &eta_wrtThr);
   if(Tree->GetBranch("rap_wrtThr") != nullptr)
      Tree->SetBranchAddress("rap_wrtThr", &rap_wrtThr);
   if(Tree->GetBranch("theta_wrtThr") != nullptr)
      Tree->SetBranchAddress("theta_wrtThr", &theta_wrtThr);
   if(Tree->GetBranch("phi_wrtThr") != nullptr)
      Tree->SetBranchAddress("phi_wrtThr", &phi_wrtThr);
   if(Tree->GetBranch("pt_wrtThrPerp") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrPerp", &pt_wrtThrPerp);
   if(Tree->GetBranch("eta_wrtThrPerp") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrPerp", &eta_wrtThrPerp);
   if(Tree->GetBranch("rap_wrtThrPerp") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrPerp", &rap_wrtThrPerp);
   if(Tree->GetBranch("theta_wrtThrPerp") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrPerp", &theta_wrtThrPerp);
   if(Tree->GetBranch("phi_wrtThrPerp") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrPerp", &phi_wrtThrPerp);
   if(Tree->GetBranch("pt_wrtChThr") != nullptr)
      Tree->SetBranchAddress("pt_wrtChThr", &pt_wrtChThr);
   if(Tree->GetBranch("eta_wrtChThr") != nullptr)
      Tree->SetBranchAddress("eta_wrtChThr", &eta_wrtChThr);
   if(Tree->GetBranch("rap_wrtChThr") != nullptr)
      Tree->SetBranchAddress("rap_wrtChThr", &rap_wrtChThr);
   if(Tree->GetBranch("theta_wrtChThr") != nullptr)
      Tree->SetBranchAddress("theta_wrtChThr", &theta_wrtChThr);
   if(Tree->GetBranch("phi_wrtChThr") != nullptr)
      Tree->SetBranchAddress("phi_wrtChThr", &phi_wrtChThr);
   if(Tree->GetBranch("pt_wrtChThrPerp") != nullptr)
      Tree->SetBranchAddress("pt_wrtChThrPerp", &pt_wrtChThrPerp);
   if(Tree->GetBranch("eta_wrtChThrPerp") != nullptr)
      Tree->SetBranchAddress("eta_wrtChThrPerp", &eta_wrtChThrPerp);
   if(Tree->GetBranch("rap_wrtChThrPerp") != nullptr)
      Tree->SetBranchAddress("rap_wrtChThrPerp", &rap_wrtChThrPerp);
   if(Tree->GetBranch("theta_wrtChThrPerp") != nullptr)
      Tree->SetBranchAddress("theta_wrtChThrPerp", &theta_wrtChThrPerp);
   if(Tree->GetBranch("phi_wrtChThrPerp") != nullptr)
      Tree->SetBranchAddress("phi_wrtChThrPerp", &phi_wrtChThrPerp);
   if(Tree->GetBranch("pt_wrtNeuThr") != nullptr)
      Tree->SetBranchAddress("pt_wrtNeuThr", &pt_wrtNeuThr);
   if(Tree->GetBranch("eta_wrtNeuThr") != nullptr)
      Tree->SetBranchAddress("eta_wrtNeuThr", &eta_wrtNeuThr);
   if(Tree->GetBranch("rap_wrtNeuThr") != nullptr)
      Tree->SetBranchAddress("rap_wrtNeuThr", &rap_wrtNeuThr);
   if(Tree->GetBranch("theta_wrtNeuThr") != nullptr)
      Tree->SetBranchAddress("theta_wrtNeuThr", &theta_wrtNeuThr);
   if(Tree->GetBranch("phi_wrtNeuThr") != nullptr)
      Tree->SetBranchAddress("phi_wrtNeuThr", &phi_wrtNeuThr);
   if(Tree->GetBranch("pt_wrtNeuThrPerp") != nullptr)
      Tree->SetBranchAddress("pt_wrtNeuThrPerp", &pt_wrtNeuThrPerp);
   if(Tree->GetBranch("eta_wrtNeuThrPerp") != nullptr)
      Tree->SetBranchAddress("eta_wrtNeuThrPerp", &eta_wrtNeuThrPerp);
   if(Tree->GetBranch("theta_wrtNeuThrPerp") != nullptr)
      Tree->SetBranchAddress("theta_wrtNeuThrPerp", &theta_wrtNeuThrPerp);
   if(Tree->GetBranch("phi_wrtNeuThrPerp") != nullptr)
      Tree->SetBranchAddress("phi_wrtNeuThrPerp", &phi_wrtNeuThrPerp);
   if(Tree->GetBranch("rap_wrtNeuThrPerp") != nullptr)
      Tree->SetBranchAddress("rap_wrtNeuThrPerp", &rap_wrtNeuThrPerp);
   if(Tree->GetBranch("pt_wrtThrCorr") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrCorr", &pt_wrtThrCorr);
   if(Tree->GetBranch("eta_wrtThrCorr") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrCorr", &eta_wrtThrCorr);
   if(Tree->GetBranch("rap_wrtThrCorr") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrCorr", &rap_wrtThrCorr);
   if(Tree->GetBranch("theta_wrtThrCorr") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrCorr", &theta_wrtThrCorr);
   if(Tree->GetBranch("phi_wrtThrCorr") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrCorr", &phi_wrtThrCorr);
   if(Tree->GetBranch("pt_wrtThrCorrPerp") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrCorrPerp", &pt_wrtThrCorrPerp);
   if(Tree->GetBranch("eta_wrtThrCorrPerp") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrCorrPerp", &eta_wrtThrCorrPerp);
   if(Tree->GetBranch("rap_wrtThrCorrPerp") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrCorrPerp", &rap_wrtThrCorrPerp);
   if(Tree->GetBranch("theta_wrtThrCorrPerp") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrCorrPerp", &theta_wrtThrCorrPerp);
   if(Tree->GetBranch("phi_wrtThrCorrPerp") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrCorrPerp", &phi_wrtThrCorrPerp);
   if(Tree->GetBranch("pt_wrtThrCorrInverse") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrCorrInverse", &pt_wrtThrCorrInverse);
   if(Tree->GetBranch("eta_wrtThrCorrInverse") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrCorrInverse", &eta_wrtThrCorrInverse);
   if(Tree->GetBranch("rap_wrtThrCorrInverse") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrCorrInverse", &rap_wrtThrCorrInverse);
   if(Tree->GetBranch("theta_wrtThrCorrInverse") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrCorrInverse", &theta_wrtThrCorrInverse);
   if(Tree->GetBranch("phi_wrtThrCorrInverse") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrCorrInverse", &phi_wrtThrCorrInverse);
   if(Tree->GetBranch("pt_wrtThrCorrInversePerp") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrCorrInversePerp", &pt_wrtThrCorrInversePerp);
   if(Tree->GetBranch("eta_wrtThrCorrInversePerp") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrCorrInversePerp", &eta_wrtThrCorrInversePerp);
   if(Tree->GetBranch("rap_wrtThrCorrInversePerp") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrCorrInversePerp", &rap_wrtThrCorrInversePerp);
   if(Tree->GetBranch("theta_wrtThrCorrInversePerp") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrCorrInversePerp", &theta_wrtThrCorrInversePerp);
   if(Tree->GetBranch("phi_wrtThrCorrInversePerp") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrCorrInversePerp", &phi_wrtThrCorrInversePerp);
   if(Tree->GetBranch("pt_wrtThrMissP") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrMissP", &pt_wrtThrMissP);
   if(Tree->GetBranch("eta_wrtThrMissP") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrMissP", &eta_wrtThrMissP);
   if(Tree->GetBranch("rap_wrtThrMissP") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrMissP", &rap_wrtThrMissP);
   if(Tree->GetBranch("theta_wrtThrMissP") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrMissP", &theta_wrtThrMissP);
   if(Tree->GetBranch("phi_wrtThrMissP") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrMissP", &phi_wrtThrMissP);
   if(Tree->GetBranch("pt_wrtThrMissPPerp") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrMissPPerp", &pt_wrtThrMissPPerp);
   if(Tree->GetBranch("eta_wrtThrMissPPerp") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrMissPPerp", &eta_wrtThrMissPPerp);
   if(Tree->GetBranch("rap_wrtThrMissPPerp") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrMissPPerp", &rap_wrtThrMissPPerp);
   if(Tree->GetBranch("theta_wrtThrMissPPerp") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrMissPPerp", &theta_wrtThrMissPPerp);
   if(Tree->GetBranch("phi_wrtThrMissPPerp") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrMissPPerp", &phi_wrtThrMissPPerp);
   if(Tree->GetBranch("passesArtificAccept") != nullptr)
      Tree->SetBranchAddress("passesArtificAccept", &passesArtificAccept);
   if(Tree->GetBranch("artificAcceptEffCorrection") != nullptr)
      Tree->SetBranchAddress("artificAcceptEffCorrection", &artificAcceptEffCorrection);
   if(Tree->GetBranch("passesNTupleAfterCut") != nullptr)
      Tree->SetBranchAddress("passesNTupleAfterCut", &passesNTupleAfterCut);
   if(Tree->GetBranch("passesTotalChgEnergyMin") != nullptr)
      Tree->SetBranchAddress("passesTotalChgEnergyMin", &passesTotalChgEnergyMin);
   if(Tree->GetBranch("passesNTrkMin") != nullptr)
      Tree->SetBranchAddress("passesNTrkMin", &passesNTrkMin);
   if(Tree->GetBranch("passesSTheta") != nullptr)
      Tree->SetBranchAddress("passesSTheta", &passesSTheta);
   else if(Tree->GetBranch("passSphericity") != nullptr)   // Hannah pythia8
      Tree->SetBranchAddress("passSphericity", &passesSTheta);
   if(Tree->GetBranch("passesMissP") != nullptr)
      Tree->SetBranchAddress("passesMissP", &passesMissP);
   if(Tree->GetBranch("passesISR") != nullptr)
      Tree->SetBranchAddress("passesISR", &passesISR);
   if(Tree->GetBranch("passesWW") != nullptr)
      Tree->SetBranchAddress("passesWW", &passesWW);
   if(Tree->GetBranch("passesNeuNch") != nullptr)
      Tree->SetBranchAddress("passesNeuNch", &passesNeuNch);
   if(Tree->GetBranch("passesAll") != nullptr)
      Tree->SetBranchAddress("passesAll", &passesAll);
   if(Tree->GetBranch("missP") != nullptr)
      Tree->SetBranchAddress("missP", &missP);
   if(Tree->GetBranch("missPt") != nullptr)
      Tree->SetBranchAddress("missPt", &missPt);
   if(Tree->GetBranch("missTheta") != nullptr)
      Tree->SetBranchAddress("missTheta", &missTheta);
   if(Tree->GetBranch("missPhi") != nullptr)
      Tree->SetBranchAddress("missPhi", &missPhi);
   if(Tree->GetBranch("missChargedP") != nullptr)
      Tree->SetBranchAddress("missChargedP", &missChargedP);
   if(Tree->GetBranch("missChargedPt") != nullptr)
      Tree->SetBranchAddress("missChargedPt", &missChargedPt);
   if(Tree->GetBranch("missChargedTheta") != nullptr)
      Tree->SetBranchAddress("missChargedTheta", &missChargedTheta);
   if(Tree->GetBranch("missChargedPhi") != nullptr)
      Tree->SetBranchAddress("missChargedPhi", &missChargedPhi);
   if(Tree->GetBranch("nChargedHadrons") != nullptr)
      Tree->SetBranchAddress("nChargedHadrons", &nChargedHadrons);
   if(Tree->GetBranch("nChargedHadronsHP") != nullptr)
      Tree->SetBranchAddress("nChargedHadronsHP", &nChargedHadronsHP);
   if(Tree->GetBranch("nChargedHadronsHP_Corrected") != nullptr)
      Tree->SetBranchAddress("nChargedHadronsHP_Corrected", &nChargedHadronsHP_Corrected);
   if(Tree->GetBranch("nChargedHadrons_GT0p4") != nullptr)
      Tree->SetBranchAddress("nChargedHadrons_GT0p4", &nChargedHadrons_GT0p4);
   if(Tree->GetBranch("nChargedHadrons_GT0p4Thrust") != nullptr)
      Tree->SetBranchAddress("nChargedHadrons_GT0p4Thrust", &nChargedHadrons_GT0p4Thrust);
   if(Tree->GetBranch("Thrust") != nullptr)
      Tree->SetBranchAddress("Thrust", &Thrust);
   if(Tree->GetBranch("TTheta") != nullptr)
      Tree->SetBranchAddress("TTheta", &TTheta);
   if(Tree->GetBranch("TPhi") != nullptr)
      Tree->SetBranchAddress("TPhi", &TPhi);
   if(Tree->GetBranch("Thrust_charged") != nullptr)
      Tree->SetBranchAddress("Thrust_charged", &Thrust_charged);
   if(Tree->GetBranch("TTheta_charged") != nullptr)
      Tree->SetBranchAddress("TTheta_charged", &TTheta_charged);
   if(Tree->GetBranch("TPhi_charged") != nullptr)
      Tree->SetBranchAddress("TPhi_charged", &TPhi_charged);
   if(Tree->GetBranch("Thrust_neutral") != nullptr)
      Tree->SetBranchAddress("Thrust_neutral", &Thrust_neutral);
   if(Tree->GetBranch("TTheta_neutral") != nullptr)
      Tree->SetBranchAddress("TTheta_neutral", &TTheta_neutral);
   if(Tree->GetBranch("TPhi_neutral") != nullptr)
      Tree->SetBranchAddress("TPhi_neutral", &TPhi_neutral);
   if(Tree->GetBranch("ThrustCorr") != nullptr)
      Tree->SetBranchAddress("ThrustCorr", &ThrustCorr);
   if(Tree->GetBranch("TThetaCorr") != nullptr)
      Tree->SetBranchAddress("TThetaCorr", &TThetaCorr);
   if(Tree->GetBranch("TPhiCorr") != nullptr)
      Tree->SetBranchAddress("TPhiCorr", &TPhiCorr);
   if(Tree->GetBranch("ThrustCorrInverse") != nullptr)
      Tree->SetBranchAddress("ThrustCorrInverse", &ThrustCorrInverse);
   if(Tree->GetBranch("TThetaCorrInverse") != nullptr)
      Tree->SetBranchAddress("TThetaCorrInverse", &TThetaCorrInverse);
   if(Tree->GetBranch("TPhiCorrInverse") != nullptr)
      Tree->SetBranchAddress("TPhiCorrInverse", &TPhiCorrInverse);
   if(Tree->GetBranch("ThrustWithMissP") != nullptr)
      Tree->SetBranchAddress("ThrustWithMissP", &ThrustWithMissP);
   if(Tree->GetBranch("TThetaWithMissP") != nullptr)
      Tree->SetBranchAddress("TThetaWithMissP", &TThetaWithMissP);
   if(Tree->GetBranch("TPhiWithMissP") != nullptr)
      Tree->SetBranchAddress("TPhiWithMissP", &TPhiWithMissP);
   if(Tree->GetBranch("Sphericity") != nullptr)
      Tree->SetBranchAddress("Sphericity", &Sphericity);
   if(Tree->GetBranch("STheta") != nullptr)
      Tree->SetBranchAddress("STheta", &STheta);
   if(Tree->GetBranch("SPhi") != nullptr)
      Tree->SetBranchAddress("SPhi", &SPhi);
   if(Tree->GetBranch("Aplanarity") != nullptr)
      Tree->SetBranchAddress("Aplanarity", &Aplanarity);
   if(Tree->GetBranch("Sphericity_linearized") != nullptr)
      Tree->SetBranchAddress("Sphericity_linearized", &Sphericity_linearized);
   if(Tree->GetBranch("STheta_linearized") != nullptr)
      Tree->SetBranchAddress("STheta_linearized", &STheta_linearized);
   if(Tree->GetBranch("SPhi_linearized") != nullptr)
      Tree->SetBranchAddress("SPhi_linearized", &SPhi_linearized);
   if(Tree->GetBranch("Aplanarity_linearized") != nullptr)
      Tree->SetBranchAddress("Aplanarity_linearized", &Aplanarity_linearized);
   if(Tree->GetBranch("C_linearized") != nullptr)
      Tree->SetBranchAddress("C_linearized", &C_linearized);
   if(Tree->GetBranch("D_linearized") != nullptr)
      Tree->SetBranchAddress("D_linearized", &D_linearized);
   if(Tree->GetBranch("passesLEP1TwoPC") != nullptr)
      Tree->SetBranchAddress("passesLEP1TwoPC", &passesLEP1TwoPC);
   else
      passesLEP1TwoPC = true;
   if(Tree->GetBranch("ThrustWithReco") != nullptr)
      Tree->SetBranchAddress("ThrustWithReco", &ThrustWithReco);
   if(Tree->GetBranch("TThetaWithReco") != nullptr)
      Tree->SetBranchAddress("TThetaWithReco", &TThetaWithReco);
   if(Tree->GetBranch("TPhiWithReco") != nullptr)
      Tree->SetBranchAddress("TPhiWithReco", &TPhiWithReco);
   if(Tree->GetBranch("pt_wrtThrWithReco") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrWithReco", &pt_wrtThrWithReco);
   if(Tree->GetBranch("eta_wrtThrWithReco") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrWithReco", &eta_wrtThrWithReco);
   if(Tree->GetBranch("rap_wrtThrWithReco") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrWithReco", &rap_wrtThrWithReco);
   if(Tree->GetBranch("theta_wrtThrWithReco") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrWithReco", &theta_wrtThrWithReco);
   if(Tree->GetBranch("phi_wrtThrWithReco") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrWithReco", &phi_wrtThrWithReco);
   if(Tree->GetBranch("ThrustWithGenIneff") != nullptr)
      Tree->SetBranchAddress("ThrustWithGenIneff", &ThrustWithGenIneff);
   if(Tree->GetBranch("TThetaWithGenIneff") != nullptr)
      Tree->SetBranchAddress("TThetaWithGenIneff", &TThetaWithGenIneff);
   if(Tree->GetBranch("TPhiWithGenIneff") != nullptr)
      Tree->SetBranchAddress("TPhiWithGenIneff", &TPhiWithGenIneff);
   if(Tree->GetBranch("pt_wrtThrWithGenIneff") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrWithGenIneff", &pt_wrtThrWithGenIneff);
   if(Tree->GetBranch("eta_wrtThrWithGenIneff") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrWithGenIneff", &eta_wrtThrWithGenIneff);
   if(Tree->GetBranch("rap_wrtThrWithGenIneff") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrWithGenIneff", &rap_wrtThrWithGenIneff);
   if(Tree->GetBranch("theta_wrtThrWithGenIneff") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrWithGenIneff", &theta_wrtThrWithGenIneff);
   if(Tree->GetBranch("phi_wrtThrWithGenIneff") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrWithGenIneff", &phi_wrtThrWithGenIneff);
   if(Tree->GetBranch("ThrustWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("ThrustWithGenIneffFake", &ThrustWithGenIneffFake);
   if(Tree->GetBranch("TThetaWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("TThetaWithGenIneffFake", &TThetaWithGenIneffFake);
   if(Tree->GetBranch("TPhiWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("TPhiWithGenIneffFake", &TPhiWithGenIneffFake);
   if(Tree->GetBranch("pt_wrtThrWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrWithGenIneffFake", &pt_wrtThrWithGenIneffFake);
   if(Tree->GetBranch("eta_wrtThrWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrWithGenIneffFake", &eta_wrtThrWithGenIneffFake);
   if(Tree->GetBranch("rap_wrtThrWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrWithGenIneffFake", &rap_wrtThrWithGenIneffFake);
   if(Tree->GetBranch("theta_wrtThrWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrWithGenIneffFake", &theta_wrtThrWithGenIneffFake);
   if(Tree->GetBranch("phi_wrtThrWithGenIneffFake") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrWithGenIneffFake", &phi_wrtThrWithGenIneffFake);
   if(Tree->GetBranch("ThrustWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("ThrustWithRecoCorr", &ThrustWithRecoCorr);
   if(Tree->GetBranch("TThetaWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("TThetaWithRecoCorr", &TThetaWithRecoCorr);
   if(Tree->GetBranch("TPhiWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("TPhiWithRecoCorr", &TPhiWithRecoCorr);
   if(Tree->GetBranch("pt_wrtThrWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrWithRecoCorr", &pt_wrtThrWithRecoCorr);
   if(Tree->GetBranch("eta_wrtThrWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrWithRecoCorr", &eta_wrtThrWithRecoCorr);
   if(Tree->GetBranch("rap_wrtThrWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrWithRecoCorr", &rap_wrtThrWithRecoCorr);
   if(Tree->GetBranch("theta_wrtThrWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrWithRecoCorr", &theta_wrtThrWithRecoCorr);
   if(Tree->GetBranch("phi_wrtThrWithRecoCorr") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrWithRecoCorr", &phi_wrtThrWithRecoCorr);
   if(Tree->GetBranch("ThrustWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("ThrustWithRecoCorrInverse", &ThrustWithRecoCorrInverse);
   if(Tree->GetBranch("TThetaWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("TThetaWithRecoCorrInverse", &TThetaWithRecoCorrInverse);
   if(Tree->GetBranch("TPhiWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("TPhiWithRecoCorrInverse", &TPhiWithRecoCorrInverse);
   if(Tree->GetBranch("pt_wrtThrWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrWithRecoCorrInverse", &pt_wrtThrWithRecoCorrInverse);
   if(Tree->GetBranch("eta_wrtThrWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrWithRecoCorrInverse", &eta_wrtThrWithRecoCorrInverse);
   if(Tree->GetBranch("rap_wrtThrWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrWithRecoCorrInverse", &rap_wrtThrWithRecoCorrInverse);
   if(Tree->GetBranch("theta_wrtThrWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrWithRecoCorrInverse", &theta_wrtThrWithRecoCorrInverse);
   if(Tree->GetBranch("phi_wrtThrWithRecoCorrInverse") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrWithRecoCorrInverse", &phi_wrtThrWithRecoCorrInverse);
   if(Tree->GetBranch("ThrustWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("ThrustWithRecoAndMissP", &ThrustWithRecoAndMissP);
   if(Tree->GetBranch("TThetaWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("TThetaWithRecoAndMissP", &TThetaWithRecoAndMissP);
   if(Tree->GetBranch("TPhiWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("TPhiWithRecoAndMissP", &TPhiWithRecoAndMissP);
   if(Tree->GetBranch("pt_wrtThrWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("pt_wrtThrWithRecoAndMissP", &pt_wrtThrWithRecoAndMissP);
   if(Tree->GetBranch("eta_wrtThrWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("eta_wrtThrWithRecoAndMissP", &eta_wrtThrWithRecoAndMissP);
   if(Tree->GetBranch("rap_wrtThrWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("rap_wrtThrWithRecoAndMissP", &rap_wrtThrWithRecoAndMissP);
   if(Tree->GetBranch("theta_wrtThrWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("theta_wrtThrWithRecoAndMissP", &theta_wrtThrWithRecoAndMissP);
   if(Tree->GetBranch("phi_wrtThrWithRecoAndMissP") != nullptr)
      Tree->SetBranchAddress("phi_wrtThrWithRecoAndMissP", &phi_wrtThrWithRecoAndMissP);

   return true;
}

bool ParticleTreeMessenger::GetEntry(int iEntry)
{
   P.clear();

   if(Tree == nullptr)
      return false;
   if(iEntry < 0)
      return false;
   if(iEntry >= GetEntries())
      return false;

   Tree->GetEntry(iEntry);

   P.resize(nParticle);
   for(int i = 0; i < nParticle; i++)
      P[i].SetXYZMass(px[i], py[i], pz[i], mass[i]);

   return true;
}

int ParticleTreeMessenger::GetEntries()
{
   if(Tree == nullptr)
      return 0;
   return Tree->GetEntries();
}

bool ParticleTreeMessenger::PassBaselineCut()
{
   if(passesLEP1TwoPC == false)
      return false;

   double SumP = 0;
   for(int i = 0; i < nParticle; i++)
      SumP = SumP + pmag[i];
   if(SumP > 200)
      return false;

   return true;
}

ReducedTreeMessenger::ReducedTreeMessenger()
{
   Tree = nullptr;
   Initialize();
}

ReducedTreeMessenger::ReducedTreeMessenger(TFile &file, std::string name)
{
   Tree = (TTree *)file.Get(name.c_str());
   Initialize();
}

ReducedTreeMessenger::ReducedTreeMessenger(TFile *file, std::string name)
{
   Tree = (TTree *)file->Get(name.c_str());
   Initialize();
}

ReducedTreeMessenger::ReducedTreeMessenger(TTree *tree)
{
   Tree = tree;
   Initialize();
}

bool ReducedTreeMessenger::Initialize(TTree *tree)
{
   Tree = tree;
   return Initialize();
}

bool ReducedTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   Tree->SetBranchAddress("N",        &N);
   Tree->SetBranchAddress("Momentum", &Momentum);
   Tree->SetBranchAddress("Mass",     &Mass);
   Tree->SetBranchAddress("Theta",    &Theta);
   Tree->SetBranchAddress("Phi",      &Phi);
   Tree->SetBranchAddress("Weight",   &Weight);
   Tree->SetBranchAddress("Charge",   &Charge);
   Tree->SetBranchAddress("PassCut",  &PassCut);

   return true;
}

bool ReducedTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;
   if(iEntry < 0)
      return false;
   if(iEntry >= GetEntries())
      return false;

   Tree->GetEntry(iEntry);

   P.resize(N);
   for(int i = 0; i < N; i++)
      P[i].SetSizeThetaPhiMass(Momentum[i], Theta[i], Phi[i], Mass[i]);
}

int ReducedTreeMessenger::GetEntries()
{
   if(Tree == nullptr)
      return 0;
   return Tree->GetEntries();
}





