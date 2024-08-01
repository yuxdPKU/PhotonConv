// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               TrackToCalo.h
 *  \brief              Track To Calo, output root file
 *  \author Antonio Silva <antonio.silva@cern.ch>, Xudong Yu <xyu3@bnl.gov>
 */

#ifndef TRACKTOCALO_H
#define TRACKTOCALO_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;

class TrackToCalo : public SubsysReco
{
 public:

  TrackToCalo(const std::string &name = "TrackToCalo", const std::string &file = "output.root");

  ~TrackToCalo() override;

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

  void EMcalRadiusUser(bool use) {m_use_emcal_radius = use;}
  void IHcalRadiusUser(bool use) {m_use_ihcal_radius = use;}
  void OHcalRadiusUser(bool use) {m_use_ohcal_radius = use;}
  void setEMcalRadius(float r) {m_emcal_radius_user = r;}
  void setIHcalRadius(float r) {m_ihcal_radius_user = r;}
  void setOHcalRadius(float r) {m_ohcal_radius_user = r;}

 private:
   int cnt = 0;
   bool m_use_emcal_radius = false;
   bool m_use_ihcal_radius = false;
   bool m_use_ohcal_radius = false;
   float m_emcal_radius_user = 93.5;
   float m_ihcal_radius_user = 117;
   float m_ohcal_radius_user = 177.423;
   std::string _outfilename;
   TFile *_outfile = nullptr;
   TTree *_tree = nullptr;

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
   std::vector<float> _track_x;
   std::vector<float> _track_y;
   std::vector<float> _track_z;
   std::vector<float> _track_quality;
   std::vector<float> _track_dcaxy;
   std::vector<float> _track_dcaz;
   std::vector<int> _track_nc_mvtx;
   std::vector<int> _track_nc_intt;
   std::vector<int> _track_nc_tpc;
   std::vector<float> _track_ptq;
   std::vector<float> _track_px;
   std::vector<float> _track_py;
   std::vector<float> _track_pz;
   std::vector<float> _track_phi_origin;
   std::vector<float> _track_eta_origin;
   std::vector<float> _track_x_origin;
   std::vector<float> _track_y_origin;
   std::vector<float> _track_z_origin;
   std::vector<float> _track_phi_emc;
   std::vector<float> _track_eta_emc;
   std::vector<float> _track_x_emc;
   std::vector<float> _track_y_emc;
   std::vector<float> _track_z_emc;
   std::vector<float> _track_phi_ihc;
   std::vector<float> _track_eta_ihc;
   std::vector<float> _track_x_ihc;
   std::vector<float> _track_y_ihc;
   std::vector<float> _track_z_ihc;
   std::vector<float> _track_phi_ohc;
   std::vector<float> _track_eta_ohc;
   std::vector<float> _track_x_ohc;
   std::vector<float> _track_y_ohc;
   std::vector<float> _track_z_ohc;

   std::vector<int> _trClus_track_id;
   std::vector<int> _trClus_type;
   std::vector<float> _trClus_x;
   std::vector<float> _trClus_y;
   std::vector<float> _trClus_z;

   std::vector<int> _emcal_id;
   std::vector<float> _emcal_phi;
   std::vector<float> _emcal_eta;
   std::vector<float> _emcal_x;
   std::vector<float> _emcal_y;
   std::vector<float> _emcal_z;
   std::vector<float> _emcal_e;
   std::vector<int> _emcal_tower_cluster_id;
   std::vector<float> _emcal_tower_e;
   std::vector<float> _emcal_tower_phi;
   std::vector<float> _emcal_tower_eta;
   std::vector<int> _emcal_tower_status;

   std::vector<int> _hcal_id;
   std::vector<float> _hcal_phi;
   std::vector<float> _hcal_eta;
   std::vector<float> _hcal_x;
   std::vector<float> _hcal_y;
   std::vector<float> _hcal_z;
   std::vector<float> _hcal_e;
   std::vector<int> _hcal_tower_cluster_id;
   std::vector<float> _hcal_tower_e;
   std::vector<float> _hcal_tower_phi;
   std::vector<float> _hcal_tower_eta;
   std::vector<int> _hcal_tower_status;

   std::vector<float> _mbd_x;
   std::vector<float> _mbd_y;
   std::vector<float> _mbd_z;

   std::vector<int> _triggers;

   std::vector<int> _ntracks;

};

#endif // TRACKTOCALO_H
