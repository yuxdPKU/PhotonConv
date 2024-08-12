/*!
 *  \file               TrackToCalo.cc
 *  \brief              Track To Calo, output root file
 *  \author Antonio Silva <antonio.silva@cern.ch>, Xudong Yu <xyu3@bnl.gov>
 */
#include "TrackToCalo.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffaobjects/EventHeaderv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackreco/ActsPropagator.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

//____________________________________________________________________________..
TrackToCalo::TrackToCalo(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr)
{
  std::cout << "TrackToCalo::TrackToCalo(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
TrackToCalo::~TrackToCalo()
{
  std::cout << "TrackToCalo::~TrackToCalo() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TrackToCalo::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "TrackToCalo::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  delete _outfile;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");
  delete _tree;
  _tree = new TTree("tree", "A tree with track/calo info");
  _tree->Branch("_vertex_id", &_vertex_id);
  _tree->Branch("_vertex_corssing", &_vertex_crossing);
  _tree->Branch("_vertex_ntracks", &_vertex_ntracks);
  _tree->Branch("_vertex_x", &_vertex_x);
  _tree->Branch("_vertex_y", &_vertex_y);
  _tree->Branch("_vertex_z", &_vertex_z);
  _tree->Branch("_cluster_x", &_cluster_x);
  _tree->Branch("_cluster_y", &_cluster_y);
  _tree->Branch("_cluster_z", &_cluster_z);
  _tree->Branch("_track_id", &_track_id);
  _tree->Branch("_track_bc", &_track_bc);
  _tree->Branch("_track_phi", &_track_phi);
  _tree->Branch("_track_eta", &_track_eta);
  _tree->Branch("_track_pcax", &_track_pcax);
  _tree->Branch("_track_pcay", &_track_pcay);
  _tree->Branch("_track_pcaz", &_track_pcaz);
  _tree->Branch("_track_vx", &_track_vx);
  _tree->Branch("_track_vy", &_track_vy);
  _tree->Branch("_track_vz", &_track_vz);
  _tree->Branch("_track_quality", &_track_quality);
  _tree->Branch("_track_dcaxy", &_track_dcaxy);
  _tree->Branch("_track_dcaz", &_track_dcaz);
  _tree->Branch("_track_nc_mvtx", &_track_nc_mvtx);
  _tree->Branch("_track_nc_intt", &_track_nc_intt);
  _tree->Branch("_track_nc_tpc", &_track_nc_tpc);
  _tree->Branch("_track_ptq", &_track_ptq);
  _tree->Branch("_track_px", &_track_px);
  _tree->Branch("_track_py", &_track_py);
  _tree->Branch("_track_pz", &_track_pz);
  _tree->Branch("_track_phi_origin", &_track_phi_origin);
  _tree->Branch("_track_eta_origin", &_track_eta_origin);
  _tree->Branch("_track_x_origin", &_track_x_origin);
  _tree->Branch("_track_y_origin", &_track_y_origin);
  _tree->Branch("_track_z_origin", &_track_z_origin);
  _tree->Branch("_track_phi_emc", &_track_phi_emc);
  _tree->Branch("_track_eta_emc", &_track_eta_emc);
  _tree->Branch("_track_x_emc", &_track_x_emc);
  _tree->Branch("_track_y_emc", &_track_y_emc);
  _tree->Branch("_track_z_emc", &_track_z_emc);
  _tree->Branch("_track_phi_ihc", &_track_phi_ihc);
  _tree->Branch("_track_eta_ihc", &_track_eta_ihc);
  _tree->Branch("_track_x_ihc", &_track_x_ihc);
  _tree->Branch("_track_y_ihc", &_track_y_ihc);
  _tree->Branch("_track_z_ihc", &_track_z_ihc);
  _tree->Branch("_track_phi_ohc", &_track_phi_ohc);
  _tree->Branch("_track_eta_ohc", &_track_eta_ohc);
  _tree->Branch("_track_x_ohc", &_track_x_ohc);
  _tree->Branch("_track_y_ohc", &_track_y_ohc);
  _tree->Branch("_track_z_ohc", &_track_z_ohc);
  _tree->Branch("_trClus_track_id", &_trClus_track_id);
  _tree->Branch("_trClus_type", &_trClus_type);
  _tree->Branch("_trClus_x", &_trClus_x);
  _tree->Branch("_trClus_y", &_trClus_y);
  _tree->Branch("_trClus_z", &_trClus_z);
  _tree->Branch("_emcal_id", &_emcal_id);
  _tree->Branch("_emcal_phi", &_emcal_phi);
  _tree->Branch("_emcal_eta", &_emcal_eta);
  _tree->Branch("_emcal_x", &_emcal_x);
  _tree->Branch("_emcal_y", &_emcal_y);
  _tree->Branch("_emcal_z", &_emcal_z);
  _tree->Branch("_emcal_e", &_emcal_e);
  _tree->Branch("_emcal_ecore", &_emcal_ecore);
  _tree->Branch("_emcal_chi2", &_emcal_chi2);
  _tree->Branch("_emcal_prob", &_emcal_prob);
  _tree->Branch("_emcal_tower_cluster_id", &_emcal_tower_cluster_id);
  _tree->Branch("_emcal_tower_e", &_emcal_tower_e);
  _tree->Branch("_emcal_tower_phi", &_emcal_tower_phi);
  _tree->Branch("_emcal_tower_eta", &_emcal_tower_eta);
  _tree->Branch("_emcal_tower_status", &_emcal_tower_status);
  _tree->Branch("_hcal_id", &_hcal_id);
  _tree->Branch("_hcal_phi", &_hcal_phi);
  _tree->Branch("_hcal_eta", &_hcal_eta);
  _tree->Branch("_hcal_x", &_hcal_x);
  _tree->Branch("_hcal_y", &_hcal_y);
  _tree->Branch("_hcal_z", &_hcal_z);
  _tree->Branch("_hcal_e", &_hcal_e);
  _tree->Branch("_hcal_tower_cluster_id", &_hcal_tower_cluster_id);
  _tree->Branch("_hcal_tower_e", &_hcal_tower_e);
  _tree->Branch("_hcal_tower_phi", &_hcal_tower_phi);
  _tree->Branch("_hcal_tower_eta", &_hcal_tower_eta);
  _tree->Branch("_hcal_tower_status", &_hcal_tower_status);
  _tree->Branch("_mbd_x", &_mbd_x);
  _tree->Branch("_mbd_y", &_mbd_y);
  _tree->Branch("_mbd_z", &_mbd_z);
  _tree->Branch("_triggers", &_triggers);
  _tree->Branch("_ntracks", &_ntracks);

  cnt=0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackToCalo::process_event(PHCompositeNode *topNode)
{
  std::cout<<"TrackToCalo::process_event event "<<cnt<<std::endl;
  cnt++;
  ResetTreeVectors();

  PHNodeIterator nodeIter(topNode);
  PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));
  if (evtNode)
  {
    EventHeaderv1* evtHeader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
    std::cout<<"runNumber = "<<evtHeader->get_RunNumber()<<" , m_evtNumber = "<<evtHeader->get_EvtSequence()<<std::endl;
  }

  bool has_vertex = false;
  GlobalVertex *mbd_vtx = nullptr;

  CLHEP::Hep3Vector vertex(0., 0., 0.);

  if(!vertexmap)
  {
    vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  }
  else
  {
    if(!vertexmap->empty())
    {
      mbd_vtx = vertexmap->begin()->second;
      if(mbd_vtx)
      {
        vertex.setX(mbd_vtx->get_x());
        vertex.setY(mbd_vtx->get_y());
        vertex.setZ(mbd_vtx->get_z());
        has_vertex = true;
      }
    }
  }

  if(has_vertex)
  {
    _mbd_x.push_back(vertex.x());
    _mbd_y.push_back(vertex.y());
    _mbd_z.push_back(vertex.z());
  }
  else
  {
    _mbd_x.push_back(NAN);
    _mbd_y.push_back(NAN);
    _mbd_z.push_back(NAN);
  }

  //SvtxVertex *svtx_vtx = nullptr;

  if(!vertexMap)
  {
    vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  }
  else
  {
    if(!vertexMap->empty())
    {
      //svtx_vtx = vertexMap->begin()->second;
      //if(svtx_vtx)
      //{
      //  vertex.setX(svtx_vtx->get_x());
      //  vertex.setY(svtx_vtx->get_y());
      //  vertex.setZ(svtx_vtx->get_z());
      //  _vertex_x.push_back(svtx_vtx->get_x());
      //  _vertex_y.push_back(svtx_vtx->get_y());
      //  _vertex_z.push_back(svtx_vtx->get_z());
      //}

      for (const auto& [key, svtx_vtx] : *vertexMap)
      {
        _vertex_id.push_back(svtx_vtx->get_id());
        _vertex_crossing.push_back(svtx_vtx->get_beam_crossing());
        _vertex_ntracks.push_back(svtx_vtx->size_tracks());
        _vertex_x.push_back(svtx_vtx->get_x());
        _vertex_y.push_back(svtx_vtx->get_y());
        _vertex_z.push_back(svtx_vtx->get_z());
      }

    }
  }

  if(!gl1Packet)
  {
    gl1Packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    std::cout << "gl1Packet not found" << std::endl;
  }
  else
  {
    auto scaled_vector = gl1Packet->getScaledVector();
    for(int i = 0; i < 32; i++)
    {
      if((scaled_vector & (int)std::pow(2,i)) != 0)
      {
        _triggers.push_back(i);
      }
    }
  }

  if(!trackMap)
  {
    trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    if(!trackMap)
    {
      std::cout << "trackMap not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _ntracks.push_back(trackMap->size());

  if (!acts_Geometry)
  {
    acts_Geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!acts_Geometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  //TrackSeedContainer *seedContainer = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");

  //if(!seedContainer)
  //{
  //  std::cout << "seedContainer not found! Aborting!" << std::endl;
  //  return Fun4AllReturnCodes::ABORTEVENT;
  //}

  if ( !clustersEM )
  {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_EM_name);
    if (!clustersEM)
    {
      std::cout << "TrackToCalo::process_event : FATAL ERROR, cannot find cluster container " << m_RawClusCont_EM_name << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if ( !clustersHAD )
  {
    clustersHAD = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_HAD_name);
    if (!clustersHAD)
    {
      std::cout << "TrackToCalo::process_event : FATAL ERROR, cannot find cluster container " << m_RawClusCont_HAD_name << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  //if(!EMCAL_RawClusters)
  //{
  //  EMCAL_RawClusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_POS_COR_CEMC");
  //  if(!EMCAL_RawClusters)
  //  {
  //    std::cout << "EMCAL_RawClusters not found! Aborting!" << std::endl;
  //    return Fun4AllReturnCodes::ABORTEVENT;
  //  }
  //}

  if(!EMCAL_Container)
  {
    EMCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if(!EMCAL_Container)
    {
      std::cout << "EMCAL_Container not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!IHCAL_Container)
  {
    IHCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if(!IHCAL_Container)
    {
      std::cout << "IHCAL_Container not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!OHCAL_Container)
  {
    OHCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if(!OHCAL_Container)
    {
      std::cout << "OHCAL_Container not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  //if(!trkrHitSet)
  //{
  //  trkrHitSet = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  //  if(!trkrHitSet)
  //  {
  //    std::cout << "trkrHitSet not found! Aborting!" << std::endl;
  //    return Fun4AllReturnCodes::ABORTEVENT;
  //  }
  //}

  if(!trkrContainer)
  {
    trkrContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if(!trkrContainer)
    {
      std::cout << "trkrContainer not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  TrkrClusterContainer::HitSetKeyList tpcHits = trkrContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
  for (auto &hsk : tpcHits)
  {
    auto range = trkrContainer->getClusters(hsk);
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      const auto cluskey = iter->first;
      const auto cluster = iter->second;  // auto cluster = clusters->findCluster(key);
      auto glob = acts_Geometry->getGlobalPosition(cluskey, cluster);
      auto sclusgx = glob.x();
      auto sclusgy = glob.y();
      auto sclusgz = glob.z();
      _cluster_x.push_back(sclusgx);
      _cluster_y.push_back(sclusgy);
      _cluster_z.push_back(sclusgz);
    }
  }

/*
  std::cout << "trkrContainer->size(): " << trkrContainer->size() << std::endl;

  TrkrClusterContainer::HitSetKeyList mvtxHits = trkrContainer->getHitSetKeys(TrkrDefs::TrkrId::mvtxId);
  TrkrClusterContainer::HitSetKeyList inttHits = trkrContainer->getHitSetKeys(TrkrDefs::TrkrId::inttId);
  TrkrClusterContainer::HitSetKeyList tpcHits = trkrContainer->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
  TrkrClusterContainer::HitSetKeyList tpotHits = trkrContainer->getHitSetKeys(TrkrDefs::TrkrId::micromegasId);

  for(auto mvtx_hit : mvtxHits)
  {
    TrkrClusterContainer::ConstRange cluster_range = trkrContainer->getClusters(mvtx_hit);
    for(TrkrClusterContainer::ConstIterator cIter = cluster_range.first; cIter != cluster_range.second; ++cIter)
    {
      auto cluster_key = cIter->first;

      if(TrkrDefs::getTrkrId(cluster_key)  == TrkrDefs::TrkrId::mvtxId)
      {
        int id = TrkrDefs::getTrkrId(cluster_key);
        std::cout << "MVTX cluster: " << id << std::endl;
      }
    }
  }


  for(auto tpc_hit : tpcHits)
  {
    TrkrClusterContainer::ConstRange cluster_range = trkrContainer->getClusters(tpc_hit);
    for(TrkrClusterContainer::ConstIterator cIter = cluster_range.first; cIter != cluster_range.second; ++cIter)
    {
      auto cluster_key = cIter->first;

      if(TrkrDefs::getTrkrId(cluster_key)  == TrkrDefs::TrkrId::tpcId)
      {
        int id = TrkrDefs::getTrkrId(cluster_key);
        std::cout << "TPC cluster: " << id << std::endl;
      }
    }
  }


  for(auto intt_hit : inttHits)
  {
    TrkrClusterContainer::ConstRange cluster_range = trkrContainer->getClusters(intt_hit);
    for(TrkrClusterContainer::ConstIterator cIter = cluster_range.first; cIter != cluster_range.second; ++cIter)
    {
      auto cluster_key = cIter->first;

      if(TrkrDefs::getTrkrId(cluster_key)  == TrkrDefs::TrkrId::inttId)
      {
        int id = TrkrDefs::getTrkrId(cluster_key);
        std::cout << "INTT cluster: " << id << std::endl;
      }
    }
  }

  for(auto tpot_hit : tpotHits)
  {
    TrkrClusterContainer::ConstRange cluster_range = trkrContainer->getClusters(tpot_hit);
    for(TrkrClusterContainer::ConstIterator cIter = cluster_range.first; cIter != cluster_range.second; ++cIter)
    {
      auto cluster_key = cIter->first;

      if(TrkrDefs::getTrkrId(cluster_key)  == TrkrDefs::TrkrId::micromegasId)
      {
        int id = TrkrDefs::getTrkrId(cluster_key);
        std::cout << "TPOT cluster: " << id << std::endl;
      }
    }
  }

  for(TrkrClusterContainer::ConstIterator cIter = cluster_range.first; cIter != cluster_range.second; ++cIter)
  {
    auto cluster_key = cIter->first;

    if(TrkrDefs::getTrkrId(cluster_key)  == TrkrDefs::TrkrId::tpcId)
    {
      std::cout << "TPC cluster: " << TrkrDefs::getTrkrId(cluster_key) << std::endl;
    }
    else
    {
      std::cout << "Silicon cluster: " << TrkrDefs::getTrkrId(cluster_key) << std::endl;
    }
  }
*/



  //TrkrClusterCrossingAssocv1 *trkrContainerCrossing = findNode::getClass<TrkrClusterCrossingAssocv1>(topNode, "TRKR_CLUSTERCROSSINGASSOC");

  //if(!trkrContainerCrossing)
  //{
  //  std::cout << "trkrContainerCrossing not found! Aborting!" << std::endl;
  //  return Fun4AllReturnCodes::ABORTEVENT;
  //}
/*
  TrkrClusterCrossingAssocv1::ConstRange crange = trkrContainerCrossing->getAll();

  for(TrkrClusterCrossingAssocv1::ConstIterator citer = crange.first; citer != crange.second; ++citer)
  {
    TrkrDefs::cluskey mykey = citer.first;
    short int bunch_crossing_number = citer.second;
    std::cout << "mykey: " << mykey << " bunch_crossing_number: " << bunch_crossing_number << std::endl;

  }
  */

  if(!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(!EMCalGeo)
    {
      std::cout << "EMCalGeo not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!IHCalGeo)
  {
    IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if(!IHCalGeo)
    {
      std::cout << "IHCalGeo not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!OHCalGeo)
  {
    OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if(!OHCalGeo)
    {
      std::cout << "OHCalGeo not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  double caloRadiusEMCal;
  double caloRadiusIHCal;
  double caloRadiusOHCal;
  if (m_use_emcal_radius)
  {
    caloRadiusEMCal = m_emcal_radius_user;
  }
  else
  {
    caloRadiusEMCal = EMCalGeo->get_radius();
  }
  if (m_use_ihcal_radius)
  {
    caloRadiusIHCal = m_ihcal_radius_user;
  }
  else
  {
    caloRadiusIHCal = IHCalGeo->get_radius();
  }
  if (m_use_ohcal_radius)
  {
    caloRadiusOHCal = m_ohcal_radius_user;
  }
  else
  {
    caloRadiusOHCal = OHCalGeo->get_radius();
  }

  SvtxTrackState *thisState = nullptr;
  SvtxTrack *track = nullptr;
  TrackSeed *seed = nullptr;
  TrackSeed *tpc_seed = nullptr;
  TrkrCluster *trkrCluster = nullptr;

  //Acts::Vector3 acts_vertex(vertex.x(), vertex.y(), vertex.z());

  for (auto &iter : *trackMap)
  {
    track = iter.second;

    if(!track) continue;

    if(track->get_pt() < m_track_pt_low_cut) continue;

    seed = track->get_silicon_seed();

    int n_mvtx_clusters = 0;
    int n_intt_clusters = 0;
    int n_tpc_clusters = 0;
    short int bunch_crossing_number = -1;

    if(!seed)
    {
      _track_nc_mvtx.push_back(0);
      _track_nc_intt.push_back(0);
    }
    else
    {
      for(auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        trkrCluster = trkrContainer->findCluster(cluster_key);
        if(!trkrCluster)
        {
          continue;
        }
        //unsigned int cluster_detector = TrkrDefs::getTrkrId(cluster_key);
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::mvtxId)
        {
          n_mvtx_clusters++;
        }
        //std::cout << "TrkrDefs::getTrkrId(cluster_key): " << cluster_detector << std::endl;
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::inttId)
        {
          n_intt_clusters++;
          //TrkrClusterCrossingAssocv1::ConstRange bc_range = trkrContainerCrossing->getCrossings(cluster_key);
          //for(TrkrClusterCrossingAssocv1::ConstIterator bcIter = bc_range.first; bcIter != bc_range.second; ++bcIter)
          //{
          //  if(bunch_crossing_number < 0) bunch_crossing_number = bcIter->second;
          //  if(bunch_crossing_number != bcIter->second)
          //  {
          //    bunch_crossing_number = -1;
          //    break;
          //  }
          //}
        }
        Acts::Vector3 global(0., 0., 0.);
        global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);
        _trClus_track_id.push_back(track->get_id());
        _trClus_type.push_back(TrkrDefs::getTrkrId(cluster_key));
        _trClus_x.push_back(global[0]);
        _trClus_y.push_back(global[1]);
        _trClus_z.push_back(global[2]);
      }
      _track_nc_mvtx.push_back(n_mvtx_clusters);
      _track_nc_intt.push_back(n_intt_clusters);
    }

    tpc_seed = track->get_tpc_seed();

    if(tpc_seed)
    {
      for(auto key_iter = tpc_seed->begin_cluster_keys(); key_iter != tpc_seed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        trkrCluster = trkrContainer->findCluster(cluster_key);
        if(!trkrCluster)
        {
          continue;
        }
        if(TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::TrkrId::tpcId)
        {
          n_tpc_clusters++;
        }
        Acts::Vector3 global(0., 0., 0.);
        global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);
        _trClus_track_id.push_back(track->get_id());
        _trClus_type.push_back(TrkrDefs::getTrkrId(cluster_key));
        _trClus_x.push_back(global[0]);
        _trClus_y.push_back(global[1]);
        _trClus_z.push_back(global[2]);
      }
    }

    // project to R=0
    thisState = track->get_state(0);

    if(!thisState)
    {
      _track_phi_origin.push_back(NAN);
      _track_eta_origin.push_back(NAN);
      _track_x_origin.push_back(NAN);
      _track_y_origin.push_back(NAN);
      _track_z_origin.push_back(NAN);
    }
    else
    {
      _track_phi_origin.push_back(atan2(thisState->get_y(), thisState->get_x()));
      _track_eta_origin.push_back(asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y())));
      _track_x_origin.push_back(thisState->get_x());
      _track_y_origin.push_back(thisState->get_y());
      _track_z_origin.push_back(thisState->get_z());
    }

    // project to R_EMCAL
    thisState = track->get_state(caloRadiusEMCal);

    if(!thisState)
    {
      _track_phi_emc.push_back(NAN);
      _track_eta_emc.push_back(NAN);
      _track_x_emc.push_back(NAN);
      _track_y_emc.push_back(NAN);
      _track_z_emc.push_back(NAN);
    }
    else
    {
      _track_phi_emc.push_back(atan2(thisState->get_y(), thisState->get_x()));
      _track_eta_emc.push_back(asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y())));
      _track_x_emc.push_back(thisState->get_x());
      _track_y_emc.push_back(thisState->get_y());
      _track_z_emc.push_back(thisState->get_z());
    }

    // project to R_IHCAL
    thisState = track->get_state(caloRadiusIHCal);

    if(!thisState)
    {
      _track_phi_ihc.push_back(NAN);
      _track_eta_ihc.push_back(NAN);
      _track_x_ihc.push_back(NAN);
      _track_y_ihc.push_back(NAN);
      _track_z_ihc.push_back(NAN);
    }
    else
    {
      _track_phi_ihc.push_back(atan2(thisState->get_y(), thisState->get_x()));
      _track_eta_ihc.push_back(asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y())));
      _track_x_ihc.push_back(thisState->get_x());
      _track_y_ihc.push_back(thisState->get_y());
      _track_z_ihc.push_back(thisState->get_z());
    }

    // project to R_OHCAL
    thisState = track->get_state(caloRadiusOHCal);

    if(!thisState)
    {
      _track_phi_ohc.push_back(NAN);
      _track_eta_ohc.push_back(NAN);
      _track_x_ohc.push_back(NAN);
      _track_y_ohc.push_back(NAN);
      _track_z_ohc.push_back(NAN);
    }
    else
    {
      _track_phi_ohc.push_back(atan2(thisState->get_y(), thisState->get_x()));
      _track_eta_ohc.push_back(asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y())));
      _track_x_ohc.push_back(thisState->get_x());
      _track_y_ohc.push_back(thisState->get_y());
      _track_z_ohc.push_back(thisState->get_z());
    }

    unsigned int m_vertexid = track->get_vertex_id();
    bool track_have_vertex = false;
    if (vertexMap)
    {
      auto vertexit = vertexMap->find(m_vertexid);
      if (vertexit != vertexMap->end())
      {
        auto svtxvertex = vertexit->second;
        m_vx = svtxvertex->get_x();
        m_vy = svtxvertex->get_y();
        m_vz = svtxvertex->get_z();
        track_have_vertex = true;
      }
    }

    if (track_have_vertex)
    {
      _track_vx.push_back(m_vx);
      _track_vy.push_back(m_vy);
      _track_vz.push_back(m_vz);
    }
    else
    {
      _track_vx.push_back(NAN);
      _track_vy.push_back(NAN);
      _track_vz.push_back(NAN);
    }

    _track_id.push_back(track->get_id());
    _track_quality.push_back(track->get_quality());
    //auto dcapair = TrackAnalysisUtils::get_dca(track, acts_vertex);
    Acts::Vector3 zero = Acts::Vector3::Zero();
    auto dcapair = TrackAnalysisUtils::get_dca(track, zero);
    _track_dcaxy.push_back(dcapair.first.first);
    _track_dcaz.push_back(dcapair.second.first);
    _track_nc_tpc.push_back(n_tpc_clusters);
    _track_bc.push_back(bunch_crossing_number);
    _track_ptq.push_back(track->get_charge()*track->get_pt());
    _track_px.push_back(track->get_px());
    _track_py.push_back(track->get_py());
    _track_pz.push_back(track->get_pz());
    _track_phi.push_back(track->get_phi());
    _track_eta.push_back(track->get_eta());
    _track_pcax.push_back(track->get_x());
    _track_pcay.push_back(track->get_y());
    _track_pcaz.push_back(track->get_z());

  }

  /*
  caloRadiusEMCal *= Acts::UnitConstants::cm;

  const auto eta = 2.5;
  const auto theta = 2. * atan(exp(-eta));
  const auto halfZ = caloRadiusEMCal / tan(theta) * Acts::UnitConstants::cm;

  auto transform = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> emcal_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(transform, caloRadiusEMCal, halfZ);


  caloRadiusIHCal *= Acts::UnitConstants::cm;

  const auto halfZ_IHCal = caloRadiusIHCal / tan(theta) * Acts::UnitConstants::cm;

  auto transform_IHCal = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> ihcal_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(transform_IHCal, caloRadiusIHCal, halfZ_IHCal);

  caloRadiusOHCal *= Acts::UnitConstants::cm;

  const auto halfZ_OHCal = caloRadiusOHCal / tan(theta) * Acts::UnitConstants::cm;

  auto transform_OHCal = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> ohcal_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(transform_OHCal, caloRadiusOHCal, halfZ_OHCal);


  TrackSeed *seed = nullptr;

  RawCluster *cluster = nullptr;

  ActsPropagator prop(acts_Geometry);

  TowerInfo *tInfo = nullptr;

  for(std::size_t i = 0; i < seedContainer->size(); i++)
  {
    seed = seedContainer->get(i);
    if(!seed)
    {
      continue;
    }

    if(seed->get_pt() < 0.5)
    {
      continue;
    }

    TrkrCluster *trkrCluster = nullptr;

    Acts::Vector3 global(0., 0., 0.);

    TrkrDefs::cluskey c_key = UINT64_MAX;

    for (auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter)
    {
      const auto& cluster_key = *key_iter;

      trkrCluster = trkrContainer->findCluster(cluster_key);

      c_key = cluster_key;

      global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);

    }

    std::shared_ptr<const Acts::Surface> surf_mvtx = acts_Geometry->maps().getSurface(c_key, trkrCluster);


    SvtxTrackState_v1 *state = new SvtxTrackState_v1(0.);
    state->set_x(global[0]);
    state->set_y(global[1]);
    state->set_z(global[2]);
    float c_phi = std::atan2(global[1], global[0]);
    float c_r = sqrt(global[0]*global[0] + global[1]*global[1]);
    float c_theta = std::atan2(c_r, global[2]);
    float c_eta = -std::log(std::tan(c_theta / 2.));

    float c_px = seed->get_pt()*std::cos(c_phi);
    float c_py = seed->get_pt()*std::sin(c_phi);
    float c_pz = seed->get_pt()*std::sinh(c_eta);

    state->set_px(c_px);
    state->set_py(c_py);
    state->set_pz(c_pz);

    auto params = prop.makeTrackParams(state, seed->get_charge(), surf_mvtx);

    if(!params.ok())
    {
      std::cout << "NOT OK!" << std::endl;
      continue;
    }

    ActsPropagator propagator(acts_Geometry);
    propagator.constField();
    propagator.verbosity(Verbosity());
    propagator.setConstFieldValue(1.4 * Acts::UnitConstants::T);

    auto result_EMCal = propagator.propagateTrackFast((const Acts::BoundTrackParameters&)params, emcal_surf);

    auto result_IHCal = propagator.propagateTrackFast((const Acts::BoundTrackParameters&)params, ihcal_surf);

    auto result_OHCal = propagator.propagateTrackFast((const Acts::BoundTrackParameters&)params, ohcal_surf);

    float EMCal_proj_phi = 9999.;
    float IHCal_proj_phi = 9999.;
    float OHCal_proj_phi = 9999.;

    if (result_EMCal.ok())
    {
      auto projectionPos = result_EMCal.value().second.position(acts_Geometry->geometry().getGeoContext());
      EMCal_proj_phi = std::atan2(projectionPos.y(),projectionPos.x());
    }
    else
    {
      std::cout << "result_EMCal failed!" << std::endl;
      continue;
    }

    if (result_IHCal.ok())
    {
      auto projectionPos = result_IHCal.value().second.position(acts_Geometry->geometry().getGeoContext());
      IHCal_proj_phi = std::atan2(projectionPos.y(),projectionPos.x());
    }
    else
    {
      std::cout << "result_IHCal failed!" << std::endl;
      continue;
    }

    if (result_OHCal.ok())
    {
      auto projectionPos = result_OHCal.value().second.position(acts_Geometry->geometry().getGeoContext());
      OHCal_proj_phi = std::atan2(projectionPos.y(),projectionPos.x());
    }
    else
    {
      std::cout << "result_OHCal failed!" << std::endl;
      continue;
    }

    if(seed->get_charge() > 0) _EMCalProjShiftPhiPos->Fill(seed->get_phi()-EMCal_proj_phi);
    if(seed->get_charge() < 0) _EMCalProjShiftPhiNeg->Fill(seed->get_phi()-EMCal_proj_phi);

    _EMCalProjShiftPhiVsSeedPt->Fill(seed->get_pt(), seed->get_phi()-EMCal_proj_phi);

    _SeedPt->Fill(seed->get_pt());



    track_mult += 1.;

    _SeedPhi->Fill(seed->get_phi());

    bool EMCal_match = false;
    bool IHCal_match = false;
    bool OHCal_match = false;

    RawClusterContainer::Range begin_end_EMC = EMCAL_RawClusters->getClusters();
    RawClusterContainer::Iterator clusIter_EMC;

    /// Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
    {
      cluster = clusIter_EMC->second;
      if(cluster->get_energy() < 0.2) continue;



      float deltaPhi = EMCal_proj_phi - RawClusterUtility::GetAzimuthAngle(*cluster, vertex);
      if(deltaPhi > M_PI)
      {
        deltaPhi = deltaPhi - 2*M_PI;
      }
      if(deltaPhi < -M_PI)
      {
        deltaPhi = deltaPhi + 2*M_PI;
      }

      _TrackEMCalPhi->Fill(deltaPhi);

      if(std::fabs(deltaPhi) < 0.2)
      {
        _TrackEMCalMatchPhi->Fill(deltaPhi);
        EMCal_match = true;
      }
    }

    //IHCal Towers
    for(unsigned int ihcal = 0; ihcal < IHCAL_Container->size(); ihcal++)
    {
      tInfo = IHCAL_Container->get_tower_at_channel(ihcal);

      if(!tInfo)
      {
        continue;
      }

      if(tInfo->get_energy() < 0.2) continue;

      unsigned int towerinfo_key = IHCAL_Container->encode_key(ihcal);
      int ti_ieta = IHCAL_Container->getTowerEtaBin(towerinfo_key);
      int ti_iphi = IHCAL_Container->getTowerPhiBin(towerinfo_key);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ti_ieta, ti_iphi);

      RawTowerGeom *tower_geom = IHCalGeo->get_tower_geometry(key);

      float deltaPhi = IHCal_proj_phi - tower_geom->get_phi();

      if(deltaPhi > M_PI)
      {
        deltaPhi = deltaPhi - 2*M_PI;
      }
      if(deltaPhi < -M_PI)
      {
        deltaPhi = deltaPhi + 2*M_PI;
      }

      _TrackIHCalPhi->Fill(deltaPhi);

      if(fabs(deltaPhi) < 0.2)
      {
        IHCal_match = true;
      }

    }

    //OHCal Towers
    for(unsigned int ohcal = 0; ohcal < OHCAL_Container->size(); ohcal++)
    {
      tInfo = OHCAL_Container->get_tower_at_channel(ohcal);

      if(!tInfo)
      {
        continue;
      }

      if(tInfo->get_energy() < 0.2) continue;

      unsigned int towerinfo_key = OHCAL_Container->encode_key(ohcal);
      int ti_ieta = OHCAL_Container->getTowerEtaBin(towerinfo_key);
      int ti_iphi = OHCAL_Container->getTowerPhiBin(towerinfo_key);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ti_ieta, ti_iphi);

      RawTowerGeom *tower_geom = OHCalGeo->get_tower_geometry(key);

      float deltaPhi = OHCal_proj_phi - tower_geom->get_phi();

      if(deltaPhi > M_PI)
      {
        deltaPhi = deltaPhi - 2*M_PI;
      }
      if(deltaPhi < -M_PI)
      {
        deltaPhi = deltaPhi + 2*M_PI;
      }

      _TrackOHCalPhi->Fill(deltaPhi);

      if(fabs(deltaPhi) < 0.2)
      {
        OHCal_match = true;
      }


    }

  }
  */

  RawCluster *cluster = nullptr;

  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < m_emcal_e_low_cut) continue;

    _emcal_id.push_back(clusIter_EMC->first);
    _emcal_e.push_back(cluster->get_energy());
    _emcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _emcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));
    _emcal_x.push_back(cluster->get_x());
    _emcal_y.push_back(cluster->get_y());
    _emcal_z.push_back(cluster->get_z());
    _emcal_ecore.push_back(cluster->get_ecore());
    _emcal_chi2.push_back(cluster->get_chi2());
    _emcal_prob.push_back(cluster->get_prob());
//std::cout<<"cluster->get_ecore() = "<<cluster->get_ecore()<<std::endl;
//std::cout<<"cluster->get_chi2() = "<<cluster->get_chi2()<<std::endl;
//std::cout<<"cluster->get_prob() = "<<cluster->get_prob()<<std::endl;

    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;

    TowerInfo *towerInfo = nullptr;

    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
    {
      _emcal_tower_cluster_id.push_back(clusIter_EMC->first);
      RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(toweriter->first);
      _emcal_tower_phi.push_back(tower_geom->get_phi());
      _emcal_tower_eta.push_back(tower_geom->get_eta());
      _emcal_tower_e.push_back(toweriter->second);
      unsigned int key = TowerInfoDefs::encode_emcal(tower_geom->get_bineta(), tower_geom->get_binphi());
      towerInfo = EMCAL_Container->get_tower_at_key(key);
      _emcal_tower_status.push_back(towerInfo->get_status());

    }
  }

  RawClusterContainer::Range begin_end_HAD = clustersHAD->getClusters();
  RawClusterContainer::Iterator clusIter_HAD;

  /// Loop over the HCal clusters
  for (clusIter_HAD = begin_end_HAD.first; clusIter_HAD != begin_end_HAD.second; ++clusIter_HAD)
  {
    cluster = clusIter_HAD->second;
    if(cluster->get_energy() < 0.2) continue;

    _hcal_id.push_back(clusIter_HAD->first);
    _hcal_e.push_back(cluster->get_energy());
    _hcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _hcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));
    _hcal_x.push_back(cluster->get_x());
    _hcal_y.push_back(cluster->get_y());
    _hcal_z.push_back(cluster->get_z());

    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;

    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
    {
      _hcal_tower_cluster_id.push_back(clusIter_HAD->first);
      RawTowerGeom *tower_geom = nullptr;
      TowerInfo *towerInfo = nullptr;

      if(RawTowerDefs::decode_caloid(toweriter->first) == RawTowerDefs::CalorimeterId::HCALOUT)
      {
        tower_geom = OHCalGeo->get_tower_geometry(toweriter->first);
        unsigned int key = TowerInfoDefs::encode_hcal(tower_geom->get_bineta(), tower_geom->get_binphi());
        towerInfo = OHCAL_Container->get_tower_at_key(key);
      }
      else if(RawTowerDefs::decode_caloid(toweriter->first) == RawTowerDefs::CalorimeterId::HCALIN)
      {
        tower_geom = IHCalGeo->get_tower_geometry(toweriter->first);
        unsigned int key = TowerInfoDefs::encode_hcal(tower_geom->get_bineta(), tower_geom->get_binphi());
        towerInfo = IHCAL_Container->get_tower_at_key(key);
      }

      _hcal_tower_phi.push_back(tower_geom->get_phi());
      _hcal_tower_eta.push_back(tower_geom->get_eta());
      _hcal_tower_e.push_back(toweriter->second);
      _hcal_tower_status.push_back(towerInfo->get_status());
    }
  }

  /*

  RawClusterContainer::Range begin_end_EMC = EMCAL_RawClusters->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < 0.2) continue;

    _emcal_e.push_back(cluster->get_energy());
    _emcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _emcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));

  }

  for(unsigned int ihcal = 0; ihcal < IHCAL_Container->size(); ihcal++)
  {
    tInfo = IHCAL_Container->get_tower_at_channel(ihcal);

    if(!tInfo)
    {
      continue;
    }

    if(tInfo->get_energy() < 0.2) continue;

    unsigned int towerinfo_key = IHCAL_Container->encode_key(ihcal);
    int ti_ieta = IHCAL_Container->getTowerEtaBin(towerinfo_key);
    int ti_iphi = IHCAL_Container->getTowerPhiBin(towerinfo_key);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ti_ieta, ti_iphi);

    RawTowerGeom *tower_geom = IHCalGeo->get_tower_geometry(key);

    _ihcal_e.push_back(tInfo->get_energy());
    _ihcal_phi.push_back(tower_geom->get_phi());
    _ihcal_eta.push_back(tower_geom->get_eta());

  }

  for(unsigned int ohcal = 0; ohcal < OHCAL_Container->size(); ohcal++)
  {
    tInfo = OHCAL_Container->get_tower_at_channel(ohcal);

    if(!tInfo)
    {
      continue;
    }

    if(tInfo->get_energy() < 0.2) continue;

    unsigned int towerinfo_key = OHCAL_Container->encode_key(ohcal);
    int ti_ieta = OHCAL_Container->getTowerEtaBin(towerinfo_key);
    int ti_iphi = OHCAL_Container->getTowerPhiBin(towerinfo_key);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ti_ieta, ti_iphi);

    RawTowerGeom *tower_geom = OHCalGeo->get_tower_geometry(key);

    _ohcal_e.push_back(tInfo->get_energy());
    _ohcal_phi.push_back(tower_geom->get_phi());
    _ohcal_eta.push_back(tower_geom->get_eta());

  }

  */

  _tree->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackToCalo::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackToCalo::ResetTreeVectors()
{
  _vertex_id.clear();
  _vertex_crossing.clear();
  _vertex_ntracks.clear();
  _vertex_x.clear();
  _vertex_y.clear();
  _vertex_z.clear();
  _cluster_x.clear();
  _cluster_y.clear();
  _cluster_z.clear();
  _track_id.clear();
  _track_bc.clear();
  _track_phi.clear();
  _track_eta.clear();
  _track_pcax.clear();
  _track_pcay.clear();
  _track_pcaz.clear();
  _track_vx.clear();
  _track_vy.clear();
  _track_vz.clear();
  _track_quality.clear();
  _track_dcaxy.clear();
  _track_dcaz.clear();
  _track_nc_mvtx.clear();
  _track_nc_intt.clear();
  _track_nc_tpc.clear();
  _track_ptq.clear();
  _track_px.clear();
  _track_py.clear();
  _track_pz.clear();
  _track_phi_origin.clear();
  _track_eta_origin.clear();
  _track_x_origin.clear();
  _track_y_origin.clear();
  _track_z_origin.clear();
  _track_phi_emc.clear();
  _track_eta_emc.clear();
  _track_x_emc.clear();
  _track_y_emc.clear();
  _track_z_emc.clear();
  _track_phi_ihc.clear();
  _track_eta_ihc.clear();
  _track_x_ihc.clear();
  _track_y_ihc.clear();
  _track_z_ihc.clear();
  _track_phi_ohc.clear();
  _track_eta_ohc.clear();
  _track_x_ohc.clear();
  _track_y_ohc.clear();
  _track_z_ohc.clear();
  _trClus_track_id.clear();
  _trClus_type.clear();
  _trClus_x.clear();
  _trClus_y.clear();
  _trClus_z.clear();
  _emcal_id.clear();
  _emcal_phi.clear();
  _emcal_eta.clear();
  _emcal_x.clear();
  _emcal_y.clear();
  _emcal_z.clear();
  _emcal_e.clear();
  _emcal_ecore.clear();
  _emcal_prob.clear();
  _emcal_chi2.clear();
  _emcal_tower_cluster_id.clear();
  _emcal_tower_e.clear();
  _emcal_tower_phi.clear();
  _emcal_tower_eta.clear();
  _emcal_tower_status.clear();
  _hcal_id.clear();
  _hcal_phi.clear();
  _hcal_eta.clear();
  _hcal_x.clear();
  _hcal_y.clear();
  _hcal_z.clear();
  _hcal_e.clear();
  _hcal_tower_cluster_id.clear();
  _hcal_tower_e.clear();
  _hcal_tower_phi.clear();
  _hcal_tower_eta.clear();
  _hcal_tower_status.clear();
  _mbd_x.clear();
  _mbd_y.clear();
  _mbd_z.clear();
  _triggers.clear();
  _ntracks.clear();
}
