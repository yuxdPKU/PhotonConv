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
  delete _tree;
  _tree = new TTree("tree", "A tree with track/calo info");
  _tree->Branch("_vertex_id", &_vertex_id);
  _tree->Branch("_vertex_x", &_vertex_x);
  _tree->Branch("_vertex_y", &_vertex_y);
  _tree->Branch("_vertex_z", &_vertex_z);
  _tree->Branch("_track_id", &_track_id);
  _tree->Branch("_track_bc", &_track_bc);
  _tree->Branch("_track_phi", &_track_phi);
  _tree->Branch("_track_eta", &_track_eta);
  _tree->Branch("_track_nc_mvtx", &_track_nc_mvtx);
  _tree->Branch("_track_nc_intt", &_track_nc_intt);
  _tree->Branch("_track_ptq", &_track_ptq);
  _tree->Branch("_track_px", &_track_px);
  _tree->Branch("_track_py", &_track_py);
  _tree->Branch("_track_pz", &_track_pz);

  _tree->Branch("_trClus_track_id", &_trClus_track_id);
  _tree->Branch("_trClus_type", &_trClus_type);
  _tree->Branch("_trClus_x", &_trClus_x);
  _tree->Branch("_trClus_y", &_trClus_y);
  _tree->Branch("_trClus_z", &_trClus_z);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackOnly::process_event(PHCompositeNode *topNode)
{
  ResetTreeVectors();

  /*
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "TrackOnly::process_event - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    assert(vertexmap);  // force quit

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (vertexmap->empty())
  {
    std::cout << "TrackOnly::process_event - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  GlobalVertex *vtx = vertexmap->begin()->second;
  if (vtx == nullptr)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */

  SvtxVertexMap *vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

  if (!vertexMap)
  {
    std::cout << "Fatal Error - vertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  CLHEP::Hep3Vector vertex(0., 0., 0.);


/*
  SvtxVertex *mvtxVertex = vertexMap->begin()->second;

  if(mvtxVertex)
  {
    std::cout << "vx: " << mvtxVertex->get_x() << " vy: " << mvtxVertex->get_y() << " vz: " << mvtxVertex->get_z() << std::endl;
  }
  else
  {
    std::cout << "no vertex object" << std::endl;
  }
*/



  SvtxTrackMap *trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if(!trackMap)
  {
    std::cout << "trackMap not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
std::cout<<"nTracks = "<<trackMap->size()<<std::endl;

  ActsGeometry *acts_Geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!acts_Geometry)
  {
    std::cout << "ActsTrackingGeometry not on node tree. Exiting." << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TrackSeedContainer *seedContainer = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");

  if(!seedContainer)
  {
    std::cout << "seedContainer not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TrkrClusterContainer *trkrContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  if(!trkrContainer)
  {
    std::cout << "trkrContainer not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TrkrClusterCrossingAssocv1 *trkrContainerCrossing = findNode::getClass<TrkrClusterCrossingAssocv1>(topNode, "TRKR_CLUSTERCROSSINGASSOC");

  if(!trkrContainerCrossing)
  {
    std::cout << "trkrContainerCrossing not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
/*
  TrkrClusterCrossingAssocv1::ConstRange crange = trkrContainerCrossing->getAll();

  for(TrkrClusterCrossingAssocv1::ConstIterator citer = crange.first; citer != crange.second; ++citer)
  {
    TrkrDefs::cluskey mykey = citer.first;
    short int bunch_crossing_number = citer.second;
    std::cout << "mykey: " << mykey << " bunch_crossing_number: " << bunch_crossing_number << std::endl;

  }
  */


  SvtxTrack *track = nullptr;

  TrackSeed *seed = nullptr;
  TrkrCluster *trkrCluster = nullptr;
  SvtxVertex *mvtxVertex = nullptr;

  for (auto &iter : *trackMap)
  {
    track = iter.second;

    if(!track) continue;

    if(track->get_pt() < 0.5) continue;

    seed = track->get_silicon_seed();

    int n_mvtx_clusters = 0;
    int n_intt_clusters = 0;
    short int bunch_crossing_number = -1;

    if(!seed)
    {
      _track_nc_mvtx.push_back(0);
      _track_nc_intt.push_back(0);
    }
    else
    {
      //std::cout << "seed->size_cluster_keys(): " << seed->size_cluster_keys() << std::endl;
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
          TrkrClusterCrossingAssocv1::ConstRange bc_range = trkrContainerCrossing->getCrossings(cluster_key);
          for(TrkrClusterCrossingAssocv1::ConstIterator bcIter = bc_range.first; bcIter != bc_range.second; ++bcIter)
          {
            if(bunch_crossing_number < 0) bunch_crossing_number = bcIter->second;
            if(bunch_crossing_number != bcIter->second)
            {
              bunch_crossing_number = -1;
              break;
            }
          }
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

    mvtxVertex = vertexMap->get(track->get_vertex_id());

    if(mvtxVertex)
    {
      _vertex_id.push_back(track->get_vertex_id());
      _vertex_x.push_back(mvtxVertex->get_x());
      _vertex_y.push_back(mvtxVertex->get_y());
      _vertex_z.push_back(mvtxVertex->get_z());
    }
    else
    {
      _vertex_id.push_back(-1);
      _vertex_x.push_back(NAN);
      _vertex_y.push_back(NAN);
      _vertex_z.push_back(NAN);
    }


    _track_id.push_back(track->get_id());
    _track_bc.push_back(bunch_crossing_number);
    _track_ptq.push_back(track->get_charge()*track->get_pt());
    _track_px.push_back(track->get_px());
    _track_py.push_back(track->get_py());
    _track_pz.push_back(track->get_pz());
    _track_phi.push_back(track->get_phi());
    _track_eta.push_back(track->get_eta());

  }

  _tree->Fill();


  return Fun4AllReturnCodes::EVENT_OK;
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
  _vertex_x.clear();
  _vertex_y.clear();
  _vertex_z.clear();
  _track_id.clear();
  _track_bc.clear();
  _track_phi.clear();
  _track_eta.clear();
  _track_nc_mvtx.clear();
  _track_nc_intt.clear();
  _track_ptq.clear();
  _track_px.clear();
  _track_py.clear();
  _track_pz.clear();

  _trClus_track_id.clear();
  _trClus_x.clear();
  _trClus_y.clear();
  _trClus_z.clear();
  _trClus_type.clear();

}
