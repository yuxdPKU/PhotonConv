#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void hadd(std::string list="runlist")
{
  //const int nruns = 33;
  //int runlist[nruns] = {51713, 51718,  51733,  51738,  51762,  51768,  51824,  51840,  51843,  51856,  51886,  51710,  51714, 51719,  51735,  51740,  51764,  51771,  51831,  51841,  51854,  51860,  51712,  51715, 51732,  51736,  51742,  51765,  51777,  51838,  51842,  51855,  51874};

  //const int nruns = 31;
  //int runlist[nruns] = {52401, 52406, 52412, 52421, 52424, 52435, 52438, 52441, 52444, 52448, 52456, 52402, 52410, 52416, 52422, 52431, 52436, 52439, 52442, 52446, 52452, 52404, 52411, 52420, 52423, 52434, 52437, 52440, 52443, 52447, 52454};

  std::vector<int> runlist = readNumberFromText(list);
  const int nruns = runlist.size();

  TChain* chain = new TChain("tree");
  for (int i=0; i<nruns; i++) chain->Add(Form("./%d/TrackCalo_*_ana.root",runlist[i]));

  TFile* outputfile = new TFile(Form("./eop_all.root"),"recreate");
  TTree* outputtree = chain->CopyTree("");

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();

}
