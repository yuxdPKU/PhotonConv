#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

void hadd_kfp()
{
  TChain* chain = new TChain("tree_KFP");
  chain->Add("../Reconstructed/TrackCalo_run*_ana.root");
  TFile* outputfile = new TFile(Form("./eop_all_ana.root"),"recreate");
  TTree* outputtree = chain->CopyTree("");
  outputfile->cd();
  outputfile->Write();
  outputfile->Close();


  TChain* chain2 = new TChain("DecayTree");
  chain2->Add("../Reconstruction/TrackCalo_run*_kfp.root");
  TFile* outputfile2 = new TFile(Form("./eop_all_kfparticle.root"),"recreate");
  TTree* outputtree2 = chain2->CopyTree("");
  outputfile2->cd();
  outputfile2->Write();
  outputfile2->Close();
}
