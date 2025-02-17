#ifndef __CALOGEOMTEST_H__
#define __CALOGEOMTEST_H__

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class RawTowerGeomContainer;
class RawTower;

class CaloGeomTest : public SubsysReco
{
 public:
  CaloGeomTest(const std::string& name = "CaloGeomTest");
  ~CaloGeomTest() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode* topNode) override;

  void setTowerGeomNodeName(const std::string& name)
  {
    m_TowerGeomNodeName = name;
  }

  void setIndexX(int ix) { _ix = ix; }
  void setIndexY(int iy) { _iy = iy; }

 protected:
  int _ix = 0;
  int _iy = 0;
  std::string m_TowerGeomNodeName = "";
};

#endif
