/*!
 *  \file               TrackOnly.cc
 *  \brief              Track only info, output root file
 *  \author Antonio Silva <antonio.silva@cern.ch>, Xudong Yu <xyu3@bnl.gov>
 */
#include "TrackOnly.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

#include <ffaobjects/EventHeaderv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrDefs.h>
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
TrackOnly::TrackOnly(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr)
{
  std::cout << "TrackOnly::TrackOnly(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
TrackOnly::~TrackOnly()
{
  std::cout << "TrackOnly::~TrackOnly() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TrackOnly::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "TrackOnly::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  delete _outfile;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");

  createBranches();

  cnt=0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void TrackOnly::createBranches()
{
  delete _tree;
  _tree = new TTree("tree", "A tree with track info");
  _tree->Branch("_runNumber", &_runNumber);
  _tree->Branch("_eventNumber", &_eventNumber);
  _tree->Branch("_vertex_id", &_vertex_id);
  _tree->Branch("_vertex_crossing", &_vertex_crossing);
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
  _tree->Branch("_track_crossing", &_track_crossing);
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
  _tree->Branch("_trClus_track_id", &_trClus_track_id);
  _tree->Branch("_trClus_type", &_trClus_type);
  _tree->Branch("_trClus_x", &_trClus_x);
  _tree->Branch("_trClus_y", &_trClus_y);
  _tree->Branch("_trClus_z", &_trClus_z);
}

//____________________________________________________________________________..
int TrackOnly::process_event(PHCompositeNode *topNode)
{
  std::cout<<"TrackOnly::process_event event "<<cnt<<std::endl;
  cnt++;
  PHNodeIterator nodeIter(topNode);
  PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));
  if (evtNode)
  {
    EventHeaderv1* evtHeader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
    std::cout<<"runNumber = "<<evtHeader->get_RunNumber()<<" , m_evtNumber = "<<evtHeader->get_EvtSequence()<<std::endl;
    _runNumber = evtHeader->get_RunNumber();
    _eventNumber = evtHeader->get_EvtSequence();
  }
  else
  {
    _runNumber = 0;
    _eventNumber = -1;
  }

  if(!vertexMap)
  {
    vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
    if(!vertexMap)
    {
      std::cout << "TrackOnly::process_event: SvtxVertexMap not found!!!" << std::endl;
    }
  }

  if(!trackMap)
  {
    trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    if(!trackMap)
    {
      std::cout << "TrackOnly::process_event: SvtxTrackMap not found!!!" << std::endl;
    }
  }

  if (!acts_Geometry)
  {
    acts_Geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!acts_Geometry)
    {
      std::cout << "TrackOnly::process_event: ActsGeometry not found!!!" << std::endl;
    }
  }

  if(!trkrContainer)
  {
    trkrContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if(!trkrContainer)
    {
      std::cout << "TrackOnly::process_event: TRKR_CLUSTER not found!!!" << std::endl;
    }
  }

/*
  if(!trkrContainerCrossing)
  {
    trkrContainerCrossing = findNode::getClass<TrkrClusterCrossingAssocv1>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
    if(!trkrContainerCrossing)
    {
      std::cout << "TrackOnly::process_event: TRKR_CLUSTERCROSSINGASSOC not found!!!" << std::endl;
    }
  }
*/

  ResetTreeVectors();

  fillTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackOnly::fillTree()
{
  if (!trackMap || !acts_Geometry || !trkrContainer)
  {
    std::cout << PHWHERE << "missing node trees, can't continue with track calo matching"
              << std::endl;
    return;
  }

  CLHEP::Hep3Vector vertex(0., 0., 0.);

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
//std::cout<<"Track id "<<track->get_id()<<" , charge = "<<track->get_charge()<<" , quality = "<<track->get_quality()<<" , ntpc = "<<n_tpc_clusters<<" , track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , tpc seed px = "<<tpc_seed->get_px()<<" , py = "<<tpc_seed->get_py()<<" , pz = "<<tpc_seed->get_pz()<<std::endl;
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
    _track_crossing.push_back(track->get_crossing());

  }

  _tree->Fill();
}

//____________________________________________________________________________..
int TrackOnly::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackOnly::ResetTreeVectors()
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
  _track_crossing.clear();
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
  _trClus_track_id.clear();
  _trClus_type.clear();
  _trClus_x.clear();
  _trClus_y.clear();
  _trClus_z.clear();
}
