#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

void EvtDisplay2D(){

  //SetsPhenixStyle();

  TChain* chain = new TChain("tree");
  chain->Add(Form("eop_kfp_unlikesign.root"));

  int _runNumber, _eventNumber;
  std::vector<float> *_ep_clus_x=0, *_ep_clus_y=0, *_ep_clus_z=0;
  std::vector<float> *_em_clus_x=0, *_em_clus_y=0, *_em_clus_z=0;
  float _ep_x_projemc, _ep_y_projemc, _ep_z_projemc;
  float _em_x_projemc, _em_y_projemc, _em_z_projemc;
  float _ep_emcal_x, _ep_emcal_y, _ep_emcal_z;
  float _em_emcal_x, _em_emcal_y, _em_emcal_z;

  chain->SetBranchAddress("_runNumber", &_runNumber);
  chain->SetBranchAddress("_eventNumber", &_eventNumber);
  chain->SetBranchAddress("_ep_clus_x", &_ep_clus_x);
  chain->SetBranchAddress("_em_clus_x", &_em_clus_x);
  chain->SetBranchAddress("_ep_clus_y", &_ep_clus_y);
  chain->SetBranchAddress("_em_clus_y", &_em_clus_y);
  chain->SetBranchAddress("_ep_clus_z", &_ep_clus_z);
  chain->SetBranchAddress("_em_clus_z", &_em_clus_z);
  chain->SetBranchAddress("_ep_x_projemc", &_ep_x_projemc);
  chain->SetBranchAddress("_em_x_projemc", &_em_x_projemc);
  chain->SetBranchAddress("_ep_y_projemc", &_ep_y_projemc);
  chain->SetBranchAddress("_em_y_projemc", &_em_y_projemc);
  chain->SetBranchAddress("_ep_z_projemc", &_ep_z_projemc);
  chain->SetBranchAddress("_em_z_projemc", &_em_z_projemc);
  chain->SetBranchAddress("_ep_emcal_x", &_ep_emcal_x);
  chain->SetBranchAddress("_em_emcal_x", &_em_emcal_x);
  chain->SetBranchAddress("_ep_emcal_y", &_ep_emcal_y);
  chain->SetBranchAddress("_em_emcal_y", &_em_emcal_y);
  chain->SetBranchAddress("_ep_emcal_z", &_ep_emcal_z);
  chain->SetBranchAddress("_em_emcal_z", &_em_emcal_z);

  int nevent  = chain->GetEntries();
  cout<<"total nevent = "<<nevent<<endl;

  std::vector<float> vec_cluster_x;
  std::vector<float> vec_cluster_y;
  std::vector<float> vec_cluster_r;
  std::vector<float> vec_cluster_z;
  std::vector<float> vec_track_proj_x;
  std::vector<float> vec_track_proj_y;
  std::vector<float> vec_track_proj_z;
  std::vector<float> vec_track_proj_r;
  std::vector<float> vec_emcal_x;
  std::vector<float> vec_emcal_y;
  std::vector<float> vec_emcal_z;
  std::vector<float> vec_emcal_r;

  vec_cluster_x.clear();
  vec_cluster_y.clear();
  vec_cluster_r.clear();
  vec_cluster_z.clear();
  vec_track_proj_x.clear();
  vec_track_proj_y.clear();
  vec_track_proj_z.clear();
  vec_track_proj_r.clear();
  vec_emcal_x.clear();
  vec_emcal_y.clear();
  vec_emcal_z.clear();
  vec_emcal_r.clear();

  int runnumber = 53741, eventid = 1035318, irow = 342;
  //int runnumber = 53741, eventid = 1291503, irow = 1533;
  //int runnumber = 53741, eventid = 1553134, irow = 2983;

  chain->GetEntry(irow);

  //if (_runNumber!=runnumber || _eventNumber!=eventid) continue;

  for (int j = 0; j < (_ep_clus_x->size()); j++)
  {
    vec_cluster_x.push_back(_ep_clus_x->at(j));
    vec_cluster_y.push_back(_ep_clus_y->at(j));
    vec_cluster_r.push_back( sqrt(pow(_ep_clus_x->at(j),2) + pow(_ep_clus_y->at(j),2)) );
    vec_cluster_z.push_back(_ep_clus_z->at(j));
    //std::cout<<"ep x "<<_ep_clus_x->at(j)<<" y "<<_ep_clus_y->at(j)<<" z "<<_ep_clus_z->at(j)<<std::endl;
  }

  for (int j = 0; j < (_em_clus_x->size()); j++)
  {
    vec_cluster_x.push_back(_em_clus_x->at(j));
    vec_cluster_y.push_back(_em_clus_y->at(j));
    vec_cluster_r.push_back( sqrt(pow(_em_clus_x->at(j),2) + pow(_em_clus_y->at(j),2)) );
    vec_cluster_z.push_back(_em_clus_z->at(j));
    //std::cout<<"em x "<<_em_clus_x->at(j)<<" y "<<_em_clus_y->at(j)<<" z "<<_em_clus_z->at(j)<<std::endl;
  }

  vec_track_proj_x.push_back(_ep_x_projemc);
  vec_track_proj_x.push_back(_em_x_projemc);
  vec_track_proj_y.push_back(_ep_y_projemc);
  vec_track_proj_y.push_back(_em_y_projemc);
  vec_track_proj_r.push_back( sqrt(pow(_ep_x_projemc,2) + pow(_ep_y_projemc,2)) );
  vec_track_proj_r.push_back( sqrt(pow(_em_x_projemc,2) + pow(_em_y_projemc,2)) );
  vec_track_proj_z.push_back(_ep_z_projemc);
  vec_track_proj_z.push_back(_em_z_projemc);

  vec_emcal_x.push_back(_ep_emcal_x);
  vec_emcal_x.push_back(_em_emcal_x);
  vec_emcal_y.push_back(_ep_emcal_y);
  vec_emcal_y.push_back(_em_emcal_y);
  vec_emcal_z.push_back(_ep_emcal_z);
  vec_emcal_z.push_back(_em_emcal_z);
  vec_emcal_r.push_back( sqrt(pow(_ep_emcal_x,2) + pow(_ep_emcal_y,2)) );
  vec_emcal_r.push_back( sqrt(pow(_em_emcal_x,2) + pow(_em_emcal_y,2)) );

  TCanvas *can_xy = new TCanvas(Form("can_xy"), "can_xy", 800, 800);
  can_xy->SetLeftMargin(0.15);
  can_xy->SetRightMargin(0.05);
  can_xy->cd();

  TGraph *gr_clus_xy = new TGraph(vec_cluster_x.size(), vec_cluster_x.data(), vec_cluster_y.data());
  gr_clus_xy->SetTitle(Form("Run %d Event %d;X [cm];Y [cm]",runnumber,eventid));
  gr_clus_xy->GetXaxis()->SetLimits(-(emcal_radius+10.), (emcal_radius+10.));
  gr_clus_xy->SetMinimum(-(emcal_radius+10.));
  gr_clus_xy->SetMaximum((emcal_radius+10.));
  gr_clus_xy->SetMarkerSize(0.5);
  gr_clus_xy->SetMarkerStyle(20);
  gr_clus_xy->SetMarkerColor(kBlue);
  gr_clus_xy->Draw("AP");

  TGraph *gr_track_proj_xy = new TGraph(vec_track_proj_x.size(), vec_track_proj_x.data(), vec_track_proj_y.data());
  gr_track_proj_xy->SetMarkerColor(kRed);
  gr_track_proj_xy->Draw("P*,same");

  TGraph *gr_emcal_xy = new TGraph(vec_emcal_x.size(), vec_emcal_x.data(), vec_emcal_y.data());
  gr_emcal_xy->SetMarkerColor(kViolet);
  gr_emcal_xy->SetMarkerStyle(22);
  gr_emcal_xy->Draw("P,same");

  // Draw ellipses using TGraph
  DrawEllipseWithTGraph(0., 0., emcal_radius - 3., kBlack);
  DrawEllipseWithTGraph(0., 0., emcal_radius + 3., kBlack);

  can_xy->Update();
  can_xy->SaveAs(Form("figure/EvtDisplay_run%d_event%d_xy.pdf",runnumber,eventid));

  TCanvas *can_rz = new TCanvas(Form("can_rz"), "can_rz", 800, 800);
  can_rz->SetLeftMargin(0.15);
  can_rz->SetRightMargin(0.05);
  can_rz->cd();

  TGraph *gr_clus_rz = new TGraph(vec_cluster_z.size(), vec_cluster_z.data(), vec_cluster_r.data());
  gr_clus_rz->SetTitle(Form("Run %d Event %d;Z [cm];R [cm]",runnumber,eventid));
  gr_clus_rz->GetXaxis()->SetLimits(-(emcal_radius+10.), (emcal_radius+10.));
  gr_clus_rz->SetMinimum(20);
  gr_clus_rz->SetMaximum(emcal_radius+10.);
  gr_clus_rz->SetMarkerSize(0.5);
  gr_clus_rz->SetMarkerStyle(20);
  gr_clus_rz->SetMarkerColor(kBlue);
  gr_clus_rz->Draw("AP");
 
  TGraph *gr_track_proj_rz = new TGraph(vec_track_proj_z.size(), vec_track_proj_z.data(), vec_track_proj_r.data());
  gr_track_proj_rz->SetMarkerColor(kRed);
  gr_track_proj_rz->Draw("P*,same");

  TGraph *gr_emcal_rz = new TGraph(vec_emcal_z.size(), vec_emcal_z.data(), vec_emcal_r.data());
  gr_emcal_rz->SetMarkerColor(kViolet);
  gr_emcal_rz->SetMarkerStyle(22);
  gr_emcal_rz->Draw("P,same");

  // Draw ellipses using TGraph
  DrawHLineWithTGraph(emcal_radius - 3., -150, 150, kBlack);
  DrawHLineWithTGraph(emcal_radius + 3., -150, 150, kBlack);

  can_rz->Update();
  can_rz->SaveAs(Form("figure/EvtDisplay_run%d_event%d_rz.pdf",runnumber,eventid));


}
