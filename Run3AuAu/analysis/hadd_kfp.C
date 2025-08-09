#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void hadd_kfp(int runnumber=53744)
{
/*
  TChain* chain = new TChain("tree_KFP");
  for (int i=0; i<nruns; i++) chain->Add(Form("./%d/TrackCalo_*_ana.root",runlist[i]));

  TFile* outputfile = new TFile(Form("./eop_all_kfp.root"),"recreate");
  TTree* outputtree = chain->CopyTree("");

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();
*/


  TChain* chain_kfp_unlikesign = new TChain("DecayTree");
  chain_kfp_unlikesign->Add(Form("../Reconstructed/%d/clusters_seeds_%d*photonconv_kfp_unlikesign.root",runnumber,runnumber));

  TFile* outputfile_kfp_unlikesign = new TFile(Form("./photonconv_kfp_unlikesign.root"),"recreate");
  TTree* outputtree_kfp_unlikesign = chain_kfp_unlikesign->CopyTree("");

  outputfile_kfp_unlikesign->cd();
  outputfile_kfp_unlikesign->Write();
  outputfile_kfp_unlikesign->Close();

}
