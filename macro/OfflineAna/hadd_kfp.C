#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void hadd_kfp()
{
  const int nruns = 67;
  int runlist[nruns] = {51452, 51560, 51605, 51718, 51771, 51874, 51490, 51562, 51606, 51719, 51777, 51886, 51491, 51563, 51607, 51732, 51824, 51493, 51564, 51613, 51733, 51831, 51505, 51565, 51614, 51735, 51838, 51508, 51569, 51615, 51736, 51840, 51509, 51570, 51616, 51738, 51841, 51510, 51571, 51619, 51740, 51842, 51511, 51572, 51710, 51742, 51843, 51512, 51575, 51712, 51762, 51854, 51518, 51576, 51713, 51764, 51855, 51520, 51603, 51714, 51765, 51856, 51559, 51604, 51715, 51768, 51860};

  TChain* chain = new TChain("tree_KFP");
  for (int i=0; i<nruns; i++) chain->Add(Form("./%d/TrackCalo_*_ana.root",runlist[i]));

  TFile* outputfile = new TFile(Form("./eop_all_kfp.root"),"recreate");
  TTree* outputtree = chain->CopyTree("");

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();

}
