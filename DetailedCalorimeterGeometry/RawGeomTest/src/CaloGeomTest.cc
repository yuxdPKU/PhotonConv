#include "CaloGeomTest.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

CaloGeomTest::CaloGeomTest(const std::string& name)
  : SubsysReco(name)
{ 
}

CaloGeomTest::~CaloGeomTest()
{ 
}


int CaloGeomTest::InitRun(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloGeomTest::process_event(PHCompositeNode *topNode)
{
  //std::cout << "CaloGeomTest::process_event()\n";
  // Load TowerGeom object
  if (m_TowerGeomNodeName.empty())
  {
    m_TowerGeomNodeName = "TOWERGEOM_CEMC";
  }
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << m_TowerGeomNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  // Loop over the towergeom values
  int ixmin = 999999;
  int iymin = 999999;
  RawTowerGeomContainer::ConstRange begin_end_geom = towergeom->get_tower_geometries();
  RawTowerGeomContainer::ConstIterator itr_geom = begin_end_geom.first;
  for (; itr_geom != begin_end_geom.second; ++itr_geom)
  {
    //std::cout << "first geom loop\n";
    RawTowerGeom *towerg = itr_geom->second;
    RawTowerDefs::keytype towerid = towerg->get_id();
    int ix = RawTowerDefs::decode_index2(towerid);  // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(towerid);  // index1 is eta in CYL
    if (ixmin > ix)
    {
      ixmin = ix;
    }
    if (iymin > iy)
    {
      iymin = iy;
    }
  }

  // Show the values 
  itr_geom = begin_end_geom.first;
  for (; itr_geom != begin_end_geom.second; ++itr_geom)
  {
    RawTowerGeom *towerg = itr_geom->second;
    RawTowerDefs::keytype towerid = towerg->get_id();
    //    int itype = towerg->get_tower_type();
    //    if( itype==2 ) { // PbSc
    int ix = RawTowerDefs::decode_index2(towerid);  // index2 is phi in CYL
    int iy = RawTowerDefs::decode_index1(towerid);  // index1 is eta in CYL
    ix -= ixmin;
    iy -= iymin;

    if ((ix == _ix) && (iy == _iy)) {
      std::cout << "Geometry tower (" << ix << ", " << iy << "):\n";
      std::cout << "center of tower: " << "(" 
                << towerg->get_center_x() << "," 
                << towerg->get_center_y() << "," 
                << towerg->get_center_z() << ")\n"; 
      std::cout << "center of tower inner face: " << "(" 
                << towerg->get_center_int_x() << "," 
                << towerg->get_center_int_y() << "," 
                << towerg->get_center_int_z() << ")\n"; 
      std::cout << "center of tower outer face: " << "(" 
                << towerg->get_center_ext_x() << "," 
                << towerg->get_center_ext_y() << "," 
                << towerg->get_center_ext_z() << ")\n"; 
      std::cout << "center of tower lateral face (low phi): " << "(" 
                << towerg->get_center_low_phi_x() << "," 
                << towerg->get_center_low_phi_y() << "," 
                << towerg->get_center_low_phi_z() << ")\n"; 
      std::cout << "center of tower lateral face (high phi): " << "(" 
                << towerg->get_center_high_phi_x() << "," 
                << towerg->get_center_high_phi_y() << "," 
                << towerg->get_center_high_phi_z() << ")\n"; 
      std::cout << "center of tower lateral face (low eta): " << "(" 
                << towerg->get_center_low_eta_x() << "," 
                << towerg->get_center_low_eta_y() << "," 
                << towerg->get_center_low_eta_z() << ")\n"; 
      std::cout << "center of tower lateral face (high eta): " << "(" 
                << towerg->get_center_high_eta_x() << "," 
                << towerg->get_center_high_eta_y() << "," 
                << towerg->get_center_high_eta_z() << ")\n"; 
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

