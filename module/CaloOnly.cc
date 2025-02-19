#include "CaloOnly.h"

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
CaloOnly::CaloOnly(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr)
{
  std::cout << "CaloOnly::CaloOnly(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloOnly::~CaloOnly()
{
  std::cout << "CaloOnly::~CaloOnly() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloOnly::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "CaloOnly::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");
  delete _tree;
  _tree = new TTree("tree", "A tree with track/calo info");
  _tree->Branch("_emcal_phi", &_emcal_phi);
  _tree->Branch("_emcal_eta", &_emcal_eta);
  _tree->Branch("_emcal_e", &_emcal_e);
  _tree->Branch("_ihcal_phi", &_ihcal_phi);
  _tree->Branch("_ihcal_eta", &_ihcal_eta);
  _tree->Branch("_ihcal_e", &_ihcal_e);
  _tree->Branch("_ohcal_phi", &_ohcal_phi);
  _tree->Branch("_ohcal_eta", &_ohcal_eta);
  _tree->Branch("_ohcal_e", &_ohcal_e);
  _tree->Branch("_mbd_z", &_mbd_z);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloOnly::process_event(PHCompositeNode *topNode)
{
  ResetTreeVectors();

  bool hasMBDvertex = true;


  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    //std::cout << "CaloOnly::process_event - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    //assert(vertexmap);  // force quit

    hasMBDvertex = false;
  }
  else if (vertexmap->empty())
  {
    //std::cout << "CaloOnly::process_event - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    hasMBDvertex = false;
  }

  if (vertexmap)
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx == nullptr)
    {
      hasMBDvertex = false;
    }

    if(hasMBDvertex == true)
    {
      _mbd_z.push_back(vtx->get_z());
    }
  }



  /*

  SvtxVertexMap *vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

  if (!vertexMap)
  {
    std::cout << "Fatal Error - vertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  CLHEP::Hep3Vector vertex(0., 0., 0.);

  */


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


/*
  RawClusterContainer *EMCAL_RawClusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_POS_COR_CEMC");

  if(!EMCAL_RawClusters)
  {
    std::cout << "EMCAL_RawClusters not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
*/
  TowerInfoContainer *EMCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *IHCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *OHCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  if(!EMCAL_Container)
  {
    std::cout << "EMCAL_Container not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if(!IHCAL_Container)
  {
    std::cout << "IHCAL_Container not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if(!OHCAL_Container)
  {
    std::cout << "OHCAL_Container not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
/*
  RawClusterContainer *clustersEM = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_EMCAL");
  RawClusterContainer *clustersHAD = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_HCAL");

  if ( !clustersEM ) {
    std::cout << "ParticleFlowReco_x::process_event : FATAL ERROR, cannot find cluster container TOPOCLUSTER_EMCAL" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if ( !clustersHAD ) {
    std::cout << "ParticleFlowReco_x::process_event : FATAL ERROR, cannot find cluster container TOPOCLUSTER_HCAL" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
*/
  RawTowerGeomContainer *EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  if(!EMCalGeo)
  {
    std::cout << "EMCalGeo not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawTowerGeomContainer *IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  if(!IHCalGeo)
  {
    std::cout << "IHCalGeo not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawTowerGeomContainer *OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if(!OHCalGeo)
  {
    std::cout << "OHCalGeo not found! Aborting!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  TowerInfo *tInfo = nullptr;
/*
  RawCluster *cluster = nullptr;
  CLHEP::Hep3Vector vertex(0., 0., 0.);

  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < 0.2) continue;

    float clusPhi = RawClusterUtility::GetAzimuthAngle(*cluster, vertex);
    float clusEta = RawClusterUtility::GetPseudorapidity(*cluster, vertex);

    _emcal_phi.push_back(clusPhi);
    _emcal_eta.push_back(clusEta);
    _emcal_e.push_back(cluster->get_energy());

  }
  */

  /*

  RawClusterContainer::Range begin_end_HAD = clustersHAD->getClusters();
  RawClusterContainer::Iterator clusIter_HAD;

  /// Loop over the EMCal clusters
  for (clusIter_HAD = begin_end_HAD.first; clusIter_HAD != begin_end_HAD.second; ++clusIter_HAD)
  {
    cluster = clusIter_HAD->second;
    if(cluster->get_energy() < 0.2) continue;

    float clusPhi = RawClusterUtility::GetAzimuthAngle(*cluster, vertex);
    float clusEta = RawClusterUtility::GetPseudorapidity(*cluster, vertex);

    _hcal_phi.push_back(clusPhi);
    _hcal_eta.push_back(clusEta);
    _hcal_e.push_back(cluster->get_energy());

  }
  */


  for(unsigned int iem = 0; iem < EMCAL_Container->size(); iem++)
  {
    tInfo = EMCAL_Container->get_tower_at_channel(iem);

    if(!tInfo)
    {
      continue;
    }

    if(!tInfo->get_isGood())
    {
      continue;
    }

    //if(tInfo->get_energy() < 0.2) continue;

    unsigned int towerinfo_key = EMCAL_Container->encode_key(iem);
    int ti_ieta = EMCAL_Container->getTowerEtaBin(towerinfo_key);
    int ti_iphi = EMCAL_Container->getTowerPhiBin(towerinfo_key);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ti_ieta, ti_iphi);

    RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(key);

    _emcal_e.push_back(tInfo->get_energy());
    _emcal_phi.push_back(tower_geom->get_phi());
    _emcal_eta.push_back(tower_geom->get_eta());

  }

  for(unsigned int ihcal = 0; ihcal < IHCAL_Container->size(); ihcal++)
  {
    tInfo = IHCAL_Container->get_tower_at_channel(ihcal);

    if(!tInfo)
    {
      continue;
    }

    if(!tInfo->get_isGood())
    {
      continue;
    }

    //if(tInfo->get_energy() < 0.2) continue;

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

    if(!tInfo->get_isGood())
    {
      continue;
    }

    //if(tInfo->get_energy() < 0.2) continue;

    unsigned int towerinfo_key = OHCAL_Container->encode_key(ohcal);
    int ti_ieta = OHCAL_Container->getTowerEtaBin(towerinfo_key);
    int ti_iphi = OHCAL_Container->getTowerPhiBin(towerinfo_key);
    const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ti_ieta, ti_iphi);

    RawTowerGeom *tower_geom = OHCalGeo->get_tower_geometry(key);

    _ohcal_e.push_back(tInfo->get_energy());
    _ohcal_phi.push_back(tower_geom->get_phi());
    _ohcal_eta.push_back(tower_geom->get_eta());

  }


  _tree->Fill();

  //_TrackMultVsEMCalTotalE->Fill(track_mult, EMCal_E);

  //std::cout << "Event #" << _e_counter << std::endl;
  //_e_counter++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloOnly::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloOnly::ResetTreeVectors()
{
  _emcal_phi.clear();
  _emcal_eta.clear();
  _emcal_e.clear();
  _ihcal_phi.clear();
  _ihcal_eta.clear();
  _ihcal_e.clear();
  _ohcal_phi.clear();
  _ohcal_eta.clear();
  _ohcal_e.clear();
  _mbd_z.clear();
}
