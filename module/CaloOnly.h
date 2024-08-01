// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOONLY_H
#define CALOONLY_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;

class CaloOnly : public SubsysReco
{
 public:

  CaloOnly(const std::string &name = "CaloOnly", const std::string &file = "output.root");

  ~CaloOnly() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void ResetTreeVectors();

 private:
   std::string _outfilename;
   TFile *_outfile = nullptr;
   TTree *_tree = nullptr;

   std::vector<float> _emcal_phi;
   std::vector<float> _emcal_eta;
   std::vector<float> _emcal_e;
   std::vector<float> _ihcal_phi;
   std::vector<float> _ihcal_eta;
   std::vector<float> _ihcal_e;
   std::vector<float> _ohcal_phi;
   std::vector<float> _ohcal_eta;
   std::vector<float> _ohcal_e;
   std::vector<float> _mbd_z;

};

#endif // TRACKTOCALO_H
