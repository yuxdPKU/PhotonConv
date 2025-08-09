// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               TrackOnly.h
 *  \brief              Track only info, output root file
 *  \author Antonio Silva <antonio.silva@cern.ch>, Xudong Yu <xyu3@bnl.gov>
 */

#ifndef TRACKONLY_H
#define TRACKONLY_H

#include <fun4all/SubsysReco.h>

#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;
class ActsGeometry;

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

  void fillTree();

  void createBranches();

  void setTrackPtLowCut(float pt) {m_track_pt_low_cut = pt;}

 private:
   int cnt = 0;
   std::string _outfilename;
   TFile *_outfile = nullptr;
   TTree *_tree = nullptr;

   int _runNumber = -9999;
   int _eventNumber = -9999;
   std::vector<int> _vertex_id;
   std::vector<int> _vertex_crossing;
   std::vector<int> _vertex_ntracks;
   std::vector<float> _vertex_x;
   std::vector<float> _vertex_y;
   std::vector<float> _vertex_z;
   std::vector<float> _cluster_x;
   std::vector<float> _cluster_y;
   std::vector<float> _cluster_z;
   std::vector<int> _track_id;
   std::vector<int> _track_bc;
   std::vector<float> _track_phi;
   std::vector<float> _track_eta;
   std::vector<float> _track_pcax;
   std::vector<float> _track_pcay;
   std::vector<float> _track_pcaz;
   std::vector<float> _track_crossing;
   std::vector<float> _track_vx;
   std::vector<float> _track_vy;
   std::vector<float> _track_vz;
   std::vector<float> _track_quality;
   std::vector<float> _track_dcaxy;
   std::vector<float> _track_dcaz;
   std::vector<int> _track_nc_mvtx;
   std::vector<int> _track_nc_intt;
   std::vector<int> _track_nc_tpc;
   std::vector<int> _track_ns_mvtx;
   std::vector<int> _track_ns_intt;
   std::vector<int> _track_ns_tpc;
   std::vector<float> _track_ptq;
   std::vector<float> _track_px;
   std::vector<float> _track_py;
   std::vector<float> _track_pz;

   std::vector<int> _trClus_track_id;
   std::vector<int> _trClus_type;
   std::vector<float> _trClus_x;
   std::vector<float> _trClus_y;
   std::vector<float> _trClus_z;

   SvtxVertexMap *vertexMap = nullptr;
   SvtxTrackMap *trackMap = nullptr;
   ActsGeometry *acts_Geometry = nullptr;
   TrkrClusterContainer *trkrContainer = nullptr;
   TrkrClusterCrossingAssocv1 *trkrContainerCrossing = nullptr;

   SvtxTrack *track = nullptr;
   SvtxVertex *vertex = nullptr;
   TrackSeed *seed = nullptr;
   TrackSeed *tpc_seed = nullptr;
   TrkrCluster *trkrCluster = nullptr;

   float m_track_pt_low_cut = 0.5;
   float m_vx = 0, m_vy = 0, m_vz = 0;

};

#endif // TRACKTOCALO_H
