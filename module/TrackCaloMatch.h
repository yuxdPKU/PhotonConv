// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file   TrackCaloMatch.h
 *  \brief  Track to Calo matching, save new SvtxTrackMap, Event display
 *  \author Xudong Yu <xyu3@bnl.gov>
 */

#ifndef TRACKCALOMATCH_H
#define TRACKCALOMATCH_H

#include <fun4all/SubsysReco.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterContainer.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHNode;
class TH1;
class TH2;
class TFile;
class TTree;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class ActsGeometry;

class TrackCaloMatch : public SubsysReco
{
 public:

  TrackCaloMatch(const std::string &name = "TrackCaloMatch");

  ~TrackCaloMatch() override;

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

  const std::string& GetTrackMapName() const {return m_trackMapName;}
  void SetTrackMapName(const std::string& name) {m_trackMapName = name;}
  const std::string& GetMyTrackMapName() const {return m_trackMapName_new;}
  void SetMyTrackMapName(const std::string& name) {m_trackMapName_new = name;}

  float PiRange(float deltaPhi)
  {
    if(deltaPhi > M_PI) deltaPhi -= 2*M_PI;
    if(deltaPhi < -M_PI) deltaPhi += 2*M_PI;
    return deltaPhi;
  }

  void writeEventDisplays( bool value ) { m_write_evt_display = value; }

  void setEventDisplayPath( const std::string& path ) { m_evt_display_path = path; }
  const std::string& getEventDisplayPath() const {return m_evt_display_path;}

  void setRunDate ( const std::string& date ) { m_run_date = date; }
  const std::string& getRunDate () const {return m_run_date;}

  void event_file_start(std::ofstream &jason_file_header, const std::string& date, int runid, int evtid);

  void doSimulation(bool set) {m_is_simulation = set;}

  void EMcalRadiusUser(bool use) {m_use_emcal_radius = use;}
  void IHcalRadiusUser(bool use) {m_use_ihcal_radius = use;}
  void OHcalRadiusUser(bool use) {m_use_ohcal_radius = use;}
  void setEMcalRadius(float r) {m_emcal_radius_user = r;}
  void setIHcalRadius(float r) {m_ihcal_radius_user = r;}
  void setOHcalRadius(float r) {m_ohcal_radius_user = r;}

  void setRawClusContEMName(const std::string& name) {m_RawClusCont_EM_name = name;}
  void setRawClusContHADName(const std::string& name) {m_RawClusCont_HAD_name = name;}
  void setRawTowerGeomContName(const std::string& name) {m_RawTowerGeomCont_name = name;}

  void setTrackPtLowCut(float pt) {m_track_pt_low_cut = pt;}
  void setEmcalELowCut(float e) {m_emcal_e_low_cut = e;}
  void setnMvtxClusters(int n) {m_nmvtx_low_cut = n;}
  void setnInttClusters(int n) {m_nintt_low_cut = n;}
  void setnTpcClusters(int n) {m_ntpc_low_cut = n;}
  void setnTpotClusters(int n) {m_ntpot_low_cut = n;}
  void setTrackQuality(float q) {m_track_quality = q;}
  void setdphicut(float a) {m_dphi_cut = a;};
  void setdzcut(float a) {m_dz_cut = a;};

 private:
  bool checkTrack(SvtxTrack* track);

  int m_runNumber = 0;
  int m_evtNumber = 0;
  bool m_use_emcal_radius = false;
  bool m_use_ihcal_radius = false;
  bool m_use_ohcal_radius = false;
  float m_emcal_radius_user = 93.5;
  float m_ihcal_radius_user = 117;
  float m_ohcal_radius_user = 177.423;
  bool m_is_simulation = false;
  PHG4TruthInfoContainer *m_truthInfo = nullptr;
  PHHepMCGenEventMap *m_geneventmap = nullptr;
  PHHepMCGenEvent *m_genevt = nullptr;
  SvtxTrackMap* trackMap = nullptr;
  SvtxTrackMap_v2* trackMap_new = nullptr;
  ActsGeometry* acts_Geometry = nullptr;
  RawClusterContainer* clustersEM = nullptr;
  RawClusterContainer* clustersHAD = nullptr;
  RawClusterContainer* EMCAL_RawClusters = nullptr;
  TowerInfoContainer* EMCAL_Container = nullptr;
  TowerInfoContainer* IHCAL_Container = nullptr;
  TowerInfoContainer* OHCAL_Container = nullptr;
  TrkrHitSetContainer* trkrHitSet = nullptr;
  TrkrClusterContainer* trkrContainer = nullptr;
  RawTowerGeomContainer* EMCalGeo = nullptr;
  RawTowerGeomContainer* IHCalGeo = nullptr;
  RawTowerGeomContainer* OHCalGeo = nullptr;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_trackMapName_new = "MySvtxTrackMap";

  std::string m_RawClusCont_EM_name = "TOPOCLUSTER_EMCAL";
  std::string m_RawClusCont_HAD_name = "TOPOCLUSTER_HCAL";
  std::string m_RawTowerGeomCont_name = "TOWERGEOM_CEMC";

  bool m_write_evt_display = false;
  std::string m_evt_display_path = "./";
  std::string m_run_date = "2024-01-01";

  float m_track_pt_low_cut = 1;
  float m_emcal_e_low_cut = 0.5;
  int m_nmvtx_low_cut = 0;
  int m_nintt_low_cut = 0;
  int m_ntpc_low_cut = 20;
  int m_ntpot_low_cut = 0;
  float m_track_quality = 1000;
  float m_dphi_cut = 0.5;
  float m_dz_cut = 20;
};

#endif
// TRACKCALOMATCH_H
