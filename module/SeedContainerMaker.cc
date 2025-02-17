/*!
 *  \file               SeedContainerMaker.h
 *  \brief              Make dummpy SiliconTrackSeedContainer, and copy TpcTrackSeedContainer to SvtxTrackSeedContainer
 *  \author Xudong Yu <xyu3@bnl.gov>
 */

#include "SeedContainerMaker.h"

#include <ffaobjects/EventHeaderv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterCrossingAssocv1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/recoConsts.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>
#include <vector>

//____________________________________________________________________________..
SeedContainerMaker::SeedContainerMaker(const std::string &name):
 SubsysReco(name)
{
  std::cout << "SeedContainerMaker::SeedContainerMaker(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
SeedContainerMaker::~SeedContainerMaker()
{
  std::cout << "SeedContainerMaker::~SeedContainerMaker() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int SeedContainerMaker::InitRun(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "SeedContainerMaker::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;

  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in SeedContainerMaker::Init");
  }

  PHNodeIterator dstIter(topNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  _track_map_tpc = findNode::getClass<TrackSeedContainer>(topNode, m_TpcTrackMapName);
  if (!_track_map_tpc)
  {
    std::cerr << PHWHERE << " ERROR: Can't find " << m_TpcTrackMapName.c_str() << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_tpc_new = findNode::getClass<TrackSeedContainer>(topNode, m_NewTpcTrackMapName);
  if (!_track_map_tpc_new)
  {
    _track_map_tpc_new = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(_track_map_tpc_new, m_NewTpcTrackMapName, "PHObject");
    svtxNode->addNode(node);
  }

  _track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, m_SiliconTrackMapName);
  if(!_track_map_silicon)
  {
    _track_map_silicon = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(_track_map_silicon, m_SiliconTrackMapName, "PHObject");
    svtxNode->addNode(node);
  }

  _svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, m_SvtxTrackMapName);
  if(!_svtx_seed_map)
  {
    _svtx_seed_map = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(_svtx_seed_map, m_SvtxTrackMapName, "PHObject");
    svtxNode->addNode(node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SeedContainerMaker::process_event(PHCompositeNode* topNode)
{

  // _track_map_tpc contains the TPC seed track stubs
  // _track_map_silicon is a dummy container
  // _svtx_seed_map contains only tpc track seeds

  // loop over the TPC track seeds
  for (unsigned int tpcid = 0; tpcid < _track_map_tpc->size(); ++tpcid)
  {
    auto svtxseed = std::make_unique<SvtxTrackSeed_v2>();
    svtxseed->set_tpc_seed_index(tpcid);
    _svtx_seed_map->insert(svtxseed.get());

    if (Verbosity() > 1)
    {
      std::cout << "  converted TPC seed id " << _svtx_seed_map->size() - 1 << " tpc id " << tpcid << std::endl;
    }
  }

  TrackSeed *tpc_seed = nullptr;

  for (auto &iter : *_track_map_tpc)
  {
    _track_map_tpc_new->insert(iter);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SeedContainerMaker::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "SeedContainerMaker::End(PHCompositeNode *topNode) Endding" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
