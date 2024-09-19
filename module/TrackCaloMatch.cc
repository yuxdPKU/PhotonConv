/*!
 *  \file               TrackCaloMatch.h
 *  \brief              Track to Calo matching, save new SvtxTrackMap, Event display
 *  \author Xudong Yu <xyu3@bnl.gov>
 */

#include "TrackCaloMatch.h"

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
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
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

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/recoConsts.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include <boost/format.hpp>
#include <boost/math/special_functions/sign.hpp>

//____________________________________________________________________________..
TrackCaloMatch::TrackCaloMatch(const std::string &name):
 SubsysReco(name)
{
  std::cout << "TrackCaloMatch::TrackCaloMatch(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
TrackCaloMatch::~TrackCaloMatch()
{
  std::cout << "TrackCaloMatch::~TrackCaloMatch() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TrackCaloMatch::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "TrackCaloMatch::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_event=0;

  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsTrkFitter::createNodes");
  }

  PHNodeIterator dstIter(topNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  trackMap_new = findNode::getClass<SvtxTrackMap_v2>(topNode, m_trackMapName_new);
  if(!trackMap_new)
  {
    trackMap_new = new SvtxTrackMap_v2;
    PHIODataNode<PHObject>* trackNode =
        new PHIODataNode<PHObject>(trackMap_new, m_trackMapName_new, "PHObject");
    svtxNode->addNode(trackNode);
    trackMap_new->clear();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackCaloMatch::process_event(PHCompositeNode* topNode)
{
  std::cout << "TrackCaloMatch::process_event event " << m_event << std::endl;

  recoConsts *rc = recoConsts::instance();
  //m_runNumber = rc->get_IntFlag("RUNNUMBER");

  if(!trackMap)
  {
    trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    if(!trackMap)
    {
      std::cout << "trackMap not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (!acts_Geometry)
  {
    acts_Geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!acts_Geometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if ( !clustersEM ) {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_EM_name);
    if (!clustersEM)
    {
      std::cout << "TrackCaloMatch::process_event : FATAL ERROR, cannot find cluster container " << m_RawClusCont_EM_name << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!trkrContainer)
  {
    trkrContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if(!trkrContainer)
    {
      std::cout << "trkrContainer not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(!EMCalGeo)
    {
      std::cout << "EMCalGeo not found! Aborting!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if(m_rejectLaserEvent)
  {
    if(!laserEventInfo)
    {
      laserEventInfo = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
      if(!laserEventInfo)
      {
        std::cout << "LaserEventInfo not found! Aborting!" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }

    if (laserEventInfo->isLaserEvent())
    {
      std::cout << "This is a laser event!" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  double caloRadiusEMCal;
  if (m_use_emcal_radius)
  {
    caloRadiusEMCal = m_emcal_radius_user;
  }
  else
  {
    caloRadiusEMCal = EMCalGeo->get_radius();
  }

  SvtxTrackState *thisState = nullptr;
  SvtxTrack *track = nullptr;
  TrackSeed *tpc_seed = nullptr;
  TrkrCluster *trkrCluster = nullptr;

  //for (auto &iter : *trackMap)
  //{
  //  SvtxTrack* kfp = iter.second;
  //  std::cout<<"iter.first = "<<iter.first<<std::endl;
  //  std::cout<<"svtxtrack id = "<<kfp->get_id()<<" px = "<<kfp->get_px()<<" py = "<<kfp->get_py()<<" pz = "<<kfp->get_pz()<<std::endl;
  //}

  int num_matched_pair = 0;
  for (auto &iter : *trackMap)
  {
    track = iter.second;

    if(!track) continue;

    if(track->get_pt() < m_track_pt_low_cut)
    {
      continue;
    }

    int n_tpc_clusters = 0;

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
      }
    }

    if(n_tpc_clusters < m_ntpc_low_cut) 
    {
      continue;
    }

    thisState = track->get_state(caloRadiusEMCal);
    float _track_phi_emc = NAN;
    //float _track_eta_emc = NAN;
    //float _track_x_emc = NAN;
    //float _track_y_emc = NAN;
    float _track_z_emc = NAN;

    if(!thisState)
    {
      continue;
    }
    else
    {
      _track_phi_emc = atan2(thisState->get_y(), thisState->get_x());
      //_track_eta_emc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
      //_track_x_emc = thisState->get_x();
      //_track_y_emc = thisState->get_y();
      _track_z_emc = thisState->get_z();
    }

    bool is_match = false;

    RawCluster *cluster = nullptr;

    RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
    RawClusterContainer::Iterator clusIter_EMC;

    /// Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
    {
      cluster = clusIter_EMC->second;
      if(cluster->get_energy() < m_emcal_e_low_cut)
      {
        continue;
      }

      float _emcal_phi = atan2(cluster->get_y(), cluster->get_x());
      //float _emcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
      float _emcal_x = cluster->get_x();
      float _emcal_y = cluster->get_y();
      float radius_scale = caloRadiusEMCal / sqrt(_emcal_x*_emcal_x+_emcal_y*_emcal_y);
      float _emcal_z = radius_scale*cluster->get_z();

      float dphi = PiRange(_track_phi_emc - _emcal_phi);
      float dz = _track_z_emc - _emcal_z;

      if(fabs(dphi)<m_dphi_cut && fabs(dz)<m_dz_cut)
      {
        //std::cout<<"matched tracks!!!"<<std::endl;
        //std::cout<<"emcal x = "<<_emcal_x<<" , y = "<<_emcal_y<<" , z = "<<_emcal_z<<" , phi = "<<_emcal_phi<<" , eta = "<<_emcal_eta<<std::endl;
        //std::cout<<"track projected x = "<<_track_x_emc<<" , y = "<<_track_y_emc<<" , z = "<<_track_z_emc<<" , phi = "<<_track_phi_emc<<" , eta = "<<_track_eta_emc<<std::endl;
        //std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;

        is_match = true;
      }
    }

    if(is_match)
    {
      //trackMap_new->insert(iter.second);
      trackMap_new->insertWithKey(iter.second,iter.first);
//std::cout<<"Hey!!!!!! insertWithKey iter.first = "<<iter.first<<" , track->get_id() = "<<track->get_id()<<std::endl;
      num_matched_pair++;
    }
  }

  //Set up the event display writer
  std::ofstream outFile;
  if (num_matched_pair>=2 && m_write_evt_display)
  {

    PHNodeIterator nodeIter(topNode);

    PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));

    if (evtNode)
    {
      EventHeaderv1* evtHeader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
      m_runNumber = evtHeader->get_RunNumber();
      m_evtNumber = evtHeader->get_EvtSequence();
    }
    else
    {
      m_runNumber = m_evtNumber = -1;
    }

    outFile.open(m_evt_display_path + "EvtDisplay_" + m_runNumber + "_" + m_evtNumber + ".json");
    event_file_start(outFile, m_run_date, m_runNumber, m_evtNumber);
    outFile << "     \"HITS\": {\n   \n     \"CEMC\":  [";
    bool firstEMCALHits = true;
    bool firstHCALHits = true;
    bool firstHits = true;

    RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
    RawClusterContainer::Iterator clusIter_EMC;

    // draw EMC clusters towers
    // without energy cuts
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
    {
      const auto cluster = clusIter_EMC->second;

      RawCluster::TowerConstRange towers = cluster->get_towers();
      RawCluster::TowerConstIterator toweriter;
      TowerInfo *towerInfo = nullptr;

      for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
      {
        RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(toweriter->first);
        std::ostringstream spts;
        if (firstEMCALHits)
        {
          firstEMCALHits = false;
        }
        else
        {
          spts << ",";
        }
        spts << "{ \"eta\": ";
        spts << tower_geom->get_eta();
        spts << ", \"phi\": ";
        spts << tower_geom->get_phi();
        spts << ", \"e\": ";
        spts << toweriter->second;
        spts << ", \"event\": ";
        spts << 0; // adcMax?
        spts << "}";

        outFile << (boost::format("%1%") % spts.str());
        spts.clear();
        spts.str("");
      }
    }
    outFile << "],\n " << std::endl;

    outFile <<  "\"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n" << std::endl;
    outFile << "    \"TRACKHITS\": [\n\n ";

    // draw all tracker clusters
    // without momentum cuts
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

        std::ostringstream spts;
        if (firstHits)
        {
          firstHits = false;
        }
        else
        {
          spts << ",";
        }
        spts << "{ \"x\": ";
        spts << sclusgx;
        spts << ", \"y\": ";
        spts << sclusgy;
        spts << ", \"z\": ";
        spts << sclusgz;
        spts << ", \"e\": ";
        spts << 0; // adcMax?
        spts << "}";

        outFile << (boost::format("%1%") % spts.str());
        spts.clear();
        spts.str("");
      }
    }

    outFile << "],\n    \"JETS\": [\n         ],\n\n" << std::endl;
    outFile << "    \"INNERTRACKER\": [\n\n ";

    firstHits = true;
    for (auto &iter : *trackMap_new)
    {
      track = iter.second;
      if(!track) continue;
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
          Acts::Vector3 global(0., 0., 0.);
          global = acts_Geometry->getGlobalPosition(cluster_key, trkrCluster);

          std::ostringstream spts;
          if (firstHits)
          {
            firstHits = false;
          }
          else
          {
            spts << ",";
          }
          spts << "{ \"x\": ";
          spts << global[0];
          spts << ", \"y\": ";
          spts << global[1];
          spts << ", \"z\": ";
          spts << global[2];
          spts << ", \"e\": ";
          spts << 0; // adcMax?
          spts << "}";

          outFile << (boost::format("%1%") % spts.str());
          spts.clear();
          spts.str("");
        }
      }
    }

    // information of fitted track by KFP
    // two tracks for dielectron
    // x,y,z is the position of secondary vertex
    // px,py,pz is the 4-momentum of track in SV
    // track bending in B-field
    outFile << "]\n},\n\"TRACKS\": {" << std::endl;
    outFile << "   \"B\": 0.000014," << std::endl;
    outFile << "   \"TRACKHITS\": [" << std::endl;
    outFile << "     {" << std::endl;
    outFile << "       \"color\": 16776960," << std::endl;
    outFile << "       \"l\": 100.8824," << std::endl;
    outFile << "       \"nh\": 60," << std::endl;
    outFile << "       \"pxyz\": [1.43648, 4.50954, -6.25811]," << std::endl;
    outFile << "       \"q\": -1," << std::endl;
    outFile << "       \"xyz\": [14.1209, 26.6029, -8.50071]" << std::endl;
    outFile << "     },\n     {" << std::endl;
    outFile << "       \"color\": 16776960," << std::endl;
    outFile << "       \"l\": 100.8824," << std::endl;
    outFile << "       \"nh\": 60," << std::endl;
    outFile << "       \"pxyz\": [0.697716, 1.20626, -1.85189]," << std::endl;
    outFile << "       \"q\": 1," << std::endl;
    outFile << "       \"xyz\": [14.1177,26.5955,-8.5073]" << std::endl;
    outFile << "     }" << std::endl;

    outFile << "]" << std::endl;
    outFile << "}" << std::endl;
    outFile << "}" << std::endl;
  }

  outFile.close();
  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackCaloMatch::event_file_start(std::ofstream &jason_file_header, std::string date, int runid, int evtid)
{
    jason_file_header << "{\n    \"EVENT\": {\n        \"runid\": " << runid << ", \n        \"evtid\": " << evtid << ", \n        \"time\": 0, \n        \"type\": \"Collision\", \n        \"s_nn\": 0, \n        \"B\": 3.0,\n        \"pv\": [0,0,0],\n        \"runstats\": [\"sPHENIX Internal\",        \n        \"200 GeV pp\",        \n        \"" << date << ", Run " << runid << "\",        \n        \"Event #" << evtid << "\"]  \n    },\n" << std::endl;

    jason_file_header << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 6.0,\n              \"color\": 16711680\n              } \n          },\n" << std::endl;
    jason_file_header << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2.0,\n              \"transparent\": 0.6,\n              \"color\": 16777215\n              } \n          },\n" << std::endl;
    jason_file_header << "          \"CEMC\": {\n              \"type\": \"PROJECTIVE\",\n              \"options\": {\n                  \"rmin\": 90,\n                  \"rmax\": 136.1,\n                  \"deta\": 0.025,\n                  \"dphi\": 0.025,\n                  \"color\": 16766464,\n                  \"transparent\": 0.6,\n                  \"scaleminmax\": true\n              }\n          },\n" << std::endl;
    jason_file_header << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << std::endl;
}

//____________________________________________________________________________..
int TrackCaloMatch::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "TrackCaloMatch::End(PHCompositeNode *topNode) Endding" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
