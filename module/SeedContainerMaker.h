// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               SeedContainerMaker.h
 *  \brief              Track to Calo matching, save new SvtxTrackMap, Event display
 *  \author Xudong Yu <xyu3@bnl.gov>
 */

#ifndef SEEDCONTAINERMAKER_H
#define SEEDCONTAINERMAKER_H

#include <fun4all/SubsysReco.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHNode;
class TrackSeedContainer;

class SeedContainerMaker : public SubsysReco
{
 public:

  SeedContainerMaker(const std::string &name = "SeedContainerMaker");

  ~SeedContainerMaker() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  const std::string& GetTpcTrackMapName() const {return m_TpcTrackMapName;}
  void SetTpcTrackMapName(const std::string &name) {m_TpcTrackMapName = name;}
  const std::string& GetNewTpcTrackMapName() const {return m_NewTpcTrackMapName;}
  void SetNewTpcTrackMapName(const std::string &name) {m_NewTpcTrackMapName = name;}
  const std::string& GetSiliconTrackMapName() const {return m_SiliconTrackMapName;}
  void SetSiliconTrackMapName(const std::string &name) {m_SiliconTrackMapName = name;}
  const std::string& GetSvtxTrackMapName() const {return m_SvtxTrackMapName;}
  void SetSvtxTrackMapName(const std::string &name) {m_SvtxTrackMapName = name;}

 private:
  TrackSeedContainer *_track_map_tpc{nullptr};
  TrackSeedContainer *_track_map_tpc_new{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeedContainer *_svtx_seed_map{nullptr};

  std::string m_TpcTrackMapName = "TpcTrackSeedContainer";
  std::string m_NewTpcTrackMapName = "PureTpcTrackSeedContainer";
  std::string m_SiliconTrackMapName = "FakeSiliconTrackSeedContainer";
  std::string m_SvtxTrackMapName = "TpcSeedConvertedSvtxTrackSeedContainer";
};

#endif
// SEEDCONTAINERMAKER_H
