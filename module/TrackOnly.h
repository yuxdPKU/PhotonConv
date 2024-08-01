// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKONLY_H
#define TRACKONLY_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;

class TrackOnly : public SubsysReco
{
 public:

  TrackOnly(const std::string &name = "TrackOnly", const std::string &file = "output.root");

  ~TrackOnly() override;

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

   std::vector<int> _vertex_id;
   std::vector<float> _vertex_x;
   std::vector<float> _vertex_y;
   std::vector<float> _vertex_z;
   std::vector<int> _track_id;
   std::vector<int> _track_bc;
   std::vector<float> _track_phi;
   std::vector<float> _track_eta;
   std::vector<int> _track_nc_mvtx;
   std::vector<int> _track_nc_intt;
   std::vector<float> _track_ptq;
   std::vector<float> _track_px;
   std::vector<float> _track_py;
   std::vector<float> _track_pz;

   std::vector<int> _trClus_track_id;
   std::vector<int> _trClus_type;
   std::vector<float> _trClus_x;
   std::vector<float> _trClus_y;
   std::vector<float> _trClus_z;


};

#endif // TRACKTOCALO_H
