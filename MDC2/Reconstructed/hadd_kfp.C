#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void hadd_kfp(std::string list="runlist")
{
  //const int nruns = 31;
  //int runlist[nruns] = {52401, 52406, 52412, 52421, 52424, 52435, 52438, 52441, 52444, 52448, 52456, 52402, 52410, 52416, 52422, 52431, 52436, 52439, 52442, 52446, 52452, 52404, 52411, 52420, 52423, 52434, 52437, 52440, 52443, 52447, 52454};

  std::vector<int> runlist = readNumberFromText(list);
  const int nruns = runlist.size();

  TChain* chain = new TChain("tree_KFP");
  for (int i=0; i<nruns; i++) chain->Add(Form("./%d/TrackCalo_*_ana.root",runlist[i]));

  TFile* outputfile = new TFile(Form("./eop_all_kfp.root"),"recreate");
  TTree* outputtree = chain->CopyTree("");

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();


  TChain* chain2 = new TChain("DecayTree");
  for (int i=0; i<nruns; i++) chain2->Add(Form("./%d/TrackCalo_*_kfp.root",runlist[i]));

  TFile* outputfile2 = new TFile(Form("./eop_all_kfparticle.root"),"recreate");
  TTree* outputtree2 = chain2->CopyTree("");

  outputfile2->cd();
  outputfile2->Write();
  outputfile2->Close();

}
