#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloGeomMappingv2.h>
#include <calogeomtest/CaloGeomTest.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libcalogeomtest.so)

void Fun4AllCaloGeomTest()
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", 24);

  // Load the modified geometry
  CaloGeomMappingv2 *cgm = new CaloGeomMappingv2();
  cgm->set_detector_name("CEMC");
  cgm->setTowerGeomNodeName("TOWERGEOM_CEMCv3");
  se->registerSubsystem(cgm);
  
  //gSystem->ListLibraries();

  // Show the geometry of a given tower
  CaloGeomTest* cgt = new CaloGeomTest();
  cgt->setTowerGeomNodeName("TOWERGEOM_CEMCv3");
  cgt->setIndexX(2);
  cgt->setIndexY(10);
  se->registerSubsystem(cgt);

  se->run(1);
  se->End();
  delete se;
  
  gSystem->Exit(0);
 }
