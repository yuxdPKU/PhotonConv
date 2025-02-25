#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void hadd_track2calo(int runnumber=53744)
{
  TChain* chain = new TChain("tree_KFP");
  chain->Add(Form("../Reconstructed/%d/clusters_seeds_%d*track2calo.root",runnumber,runnumber));

  TFile* outputfile = new TFile(Form("./photonconv_track2calo.root"),"recreate");
  TTree* outputtree = chain->CopyTree("");

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();

}
