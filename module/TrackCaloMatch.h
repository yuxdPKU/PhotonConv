// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               TrackCaloMatch.h
 *  \brief              Track to Calo matching, save new SvtxTrackMap, Event display
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

#include <string>
#include <vector>

class PHCompositeNode;
class PHNode;
class TH1;
class TH2;
class TFile;
class TTree;

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

  std::string GetTrackMapName() {return m_trackMapName;}
  void SetTrackMapName(std::string name) {m_trackMapName = name;}
  std::string GetMyTrackMapName() {return m_trackMapName_new;}
  void SetMyTrackMapName(std::string name) {m_trackMapName_new = name;}

  float PiRange(float deltaPhi)
  {
    if(deltaPhi > M_PI) deltaPhi -= 2*M_PI;
    if(deltaPhi < -M_PI) deltaPhi += 2*M_PI;
    return deltaPhi;
  }

  void writeEventDisplays( bool value ) { m_write_evt_display = value; }

  void setEventDisplayPath( std::string path ) { m_evt_display_path = path; }
  std::string getEventDisplayPath() {return m_evt_display_path;}

  void setRunDate ( std::string date ) { m_run_date = date; }
  std::string getRunDate () {return m_run_date;}

  void event_file_start(std::ofstream &jason_file_header, std::string date, int runid, int evtid);

  void EMcalRadiusUser(bool use) {m_use_emcal_radius = use;}
  void IHcalRadiusUser(bool use) {m_use_ihcal_radius = use;}
  void OHcalRadiusUser(bool use) {m_use_ohcal_radius = use;}
  void setEMcalRadius(float r) {m_emcal_radius_user = r;}
  void setIHcalRadius(float r) {m_ihcal_radius_user = r;}
  void setOHcalRadius(float r) {m_ohcal_radius_user = r;}

 private:
   int m_runNumber = 0;
   int m_evtNumber = 0;
   int m_event = 0;
   bool m_use_emcal_radius = false;
   bool m_use_ihcal_radius = false;
   bool m_use_ohcal_radius = false;
   float m_emcal_radius_user = 93.5;
   float m_ihcal_radius_user = 117;
   float m_ohcal_radius_user = 177.423;
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

   bool m_write_evt_display;
   std::string m_evt_display_path;
   std::string m_run_date;

};

#endif
// TRACKCALOMATCH_H
