#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void hadd_track2calo()
{
  const int nrun=4;
  int runs[4]={53741,53742,53743,53744};

  TChain* chain_unlikesign = new TChain("tree_KFP");
  for (int i=0; i<nrun; i++)
  {
    chain_unlikesign->Add(Form("../Reconstructed/%d/clusters_seeds_%d*track2calo_unlikesign.root",runs[i],runs[i]));
  }
  
  TFile* outputfile_unlikesign = new TFile(Form("./photonconv_track2calo_unlikesign.root"),"recreate");
  TTree* outputtree_unlikesign = chain_unlikesign->CopyTree("");
  outputfile_unlikesign->cd();
  outputfile_unlikesign->Write();
  outputfile_unlikesign->Close();


  TChain* chain_likesign = new TChain("tree_KFP");
  for (int i=0; i<nrun; i++)
  {
    chain_likesign->Add(Form("../Reconstructed/%d/clusters_seeds_%d*track2calo_likesign.root",runs[i],runs[i]));
  }
  
  TFile* outputfile_likesign = new TFile(Form("./photonconv_track2calo_likesign.root"),"recreate");
  TTree* outputtree_likesign = chain_likesign->CopyTree("");
  outputfile_likesign->cd();
  outputfile_likesign->Write();
  outputfile_likesign->Close();

}
