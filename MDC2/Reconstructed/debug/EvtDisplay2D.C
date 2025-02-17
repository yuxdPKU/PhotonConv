#include <filesystem>
#include "../utilities.h"
#include <sPhenixStyle.C>

void EvtDisplay2D(){

  //SetsPhenixStyle();

  TChain* chain = new TChain("tree");
  //chain->Add(Form("../checkevent/eop_44_kfp_v2.root"));
  chain->Add(Form("../checkevent611/eop_15_kfp_v2.root"));

  int _runNumber, _eventNumber;
  std::vector<float> *_ep_clus_x=0, *_ep_clus_y=0, *_ep_clus_z=0;
  std::vector<float> *_em_clus_x=0, *_em_clus_y=0, *_em_clus_z=0;
  float _ep_x_projemc, _ep_y_projemc, _ep_z_projemc;
  float _em_x_projemc, _em_y_projemc, _em_z_projemc;
  float _ep_emcal_x, _ep_emcal_y, _ep_emcal_z;
  float _em_emcal_x, _em_emcal_y, _em_emcal_z;
  float _ep_px, _ep_py, _ep_pz;
  float _em_px, _em_py, _em_pz;
  float _ep_x, _ep_y, _ep_z;
  float _em_x, _em_y, _em_z;
  float _ep_x_corr2, _ep_y_corr2, _ep_z_corr2;
  float _em_x_corr2, _em_y_corr2, _em_z_corr2;
  std::vector<float> *_true_gamma_decay_x=0, *_true_gamma_decay_y=0, *_true_gamma_decay_z=0;

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
  chain->SetBranchAddress("_ep_px", &_ep_px);
  chain->SetBranchAddress("_em_px", &_em_px);
  chain->SetBranchAddress("_ep_py", &_ep_py);
  chain->SetBranchAddress("_em_py", &_em_py);
  chain->SetBranchAddress("_ep_pz", &_ep_pz);
  chain->SetBranchAddress("_em_pz", &_em_pz);
  chain->SetBranchAddress("_ep_x", &_ep_x);
  chain->SetBranchAddress("_em_x", &_em_x);
  chain->SetBranchAddress("_ep_y", &_ep_y);
  chain->SetBranchAddress("_em_y", &_em_y);
  chain->SetBranchAddress("_ep_z", &_ep_z);
  chain->SetBranchAddress("_em_z", &_em_z);
  chain->SetBranchAddress("_ep_x_corr2", &_ep_x_corr2);
  chain->SetBranchAddress("_em_x_corr2", &_em_x_corr2);
  chain->SetBranchAddress("_ep_y_corr2", &_ep_y_corr2);
  chain->SetBranchAddress("_em_y_corr2", &_em_y_corr2);
  chain->SetBranchAddress("_ep_z_corr2", &_ep_z_corr2);
  chain->SetBranchAddress("_em_z_corr2", &_em_z_corr2);
  chain->SetBranchAddress("_true_gamma_decay_x", &_true_gamma_decay_x);
  chain->SetBranchAddress("_true_gamma_decay_y", &_true_gamma_decay_y);
  chain->SetBranchAddress("_true_gamma_decay_z", &_true_gamma_decay_z);

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
  std::vector<float> vec_track_x;
  std::vector<float> vec_track_y;
  std::vector<float> vec_track_z;
  std::vector<float> vec_track_r;
  std::vector<float> vec_track_x_corr2;
  std::vector<float> vec_track_y_corr2;
  std::vector<float> vec_track_z_corr2;
  std::vector<float> vec_track_r_corr2;
  std::vector<float> vec_truth_x;
  std::vector<float> vec_truth_y;
  std::vector<float> vec_truth_z;
  std::vector<float> vec_truth_r;

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
  vec_track_x.clear();
  vec_track_y.clear();
  vec_track_z.clear();
  vec_track_r.clear();
  vec_track_x_corr2.clear();
  vec_track_y_corr2.clear();
  vec_track_z_corr2.clear();
  vec_track_r_corr2.clear();
  vec_truth_x.clear();
  vec_truth_y.clear();
  vec_truth_z.clear();
  vec_truth_r.clear();

  int runnumber = 15;
  //int eventid = 179;
  //int eventid = 951;
  //int eventid = 829;
  //int eventid = 86;
  int eventid = 430;

  for(int i = 0; i < nevent; i++)
  //for(int i = 0; i < 1; i++)
  //for(int i = 1; i < 2; i++)
  {
  //  if (i % (nevent / 10) == 0) cout << "Processing progress: " << i / (nevent / 10) << "0%" << endl;
    chain->GetEntry(i);

    if (_runNumber!=runnumber || _eventNumber!=eventid) continue;
cout<<"px = "<<_ep_px<<endl;
cout<<"py = "<<_ep_py<<endl;
cout<<"pz = "<<_ep_pz<<endl;

    for (int j = 0; j < (_ep_clus_x->size()); j++)
    {
      vec_cluster_x.push_back(_ep_clus_x->at(j));
      vec_cluster_y.push_back(_ep_clus_y->at(j));
      vec_cluster_r.push_back( sqrt(pow(_ep_clus_x->at(j),2) + pow(_ep_clus_y->at(j),2)) );
      vec_cluster_z.push_back(_ep_clus_z->at(j));
    }

    for (int j = 0; j < (_em_clus_x->size()); j++)
    {
      vec_cluster_x.push_back(_em_clus_x->at(j));
      vec_cluster_y.push_back(_em_clus_y->at(j));
      vec_cluster_r.push_back( sqrt(pow(_em_clus_x->at(j),2) + pow(_em_clus_y->at(j),2)) );
      vec_cluster_z.push_back(_em_clus_z->at(j));
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
cout<<"_ep_emcal_x = "<<_ep_emcal_x<<" _ep_emcal_y = "<<_ep_emcal_y<<" _ep_emcal_z = "<<_ep_emcal_z<<endl;
cout<<"_em_emcal_x = "<<_em_emcal_x<<" _em_emcal_y = "<<_em_emcal_y<<" _em_emcal_z = "<<_em_emcal_z<<endl;
    vec_emcal_r.push_back( sqrt(pow(_ep_emcal_x,2) + pow(_ep_emcal_y,2)) );
    vec_emcal_r.push_back( sqrt(pow(_em_emcal_x,2) + pow(_em_emcal_y,2)) );
    vec_track_x.push_back(_ep_x);
    vec_track_x.push_back(_em_x);
    vec_track_y.push_back(_ep_y);
    vec_track_y.push_back(_em_y);
    vec_track_z.push_back(_ep_z);
    vec_track_z.push_back(_em_z);
    vec_track_r.push_back( sqrt(_ep_x*_ep_x+_ep_y*_ep_y) );
    vec_track_r.push_back( sqrt(_em_x*_em_x+_em_y*_em_y) );
    vec_track_x_corr2.push_back(_ep_x_corr2);
    vec_track_x_corr2.push_back(_em_x_corr2);
    vec_track_y_corr2.push_back(_ep_y_corr2);
    vec_track_y_corr2.push_back(_em_y_corr2);
    vec_track_z_corr2.push_back(_ep_z_corr2);
    vec_track_z_corr2.push_back(_em_z_corr2);
    vec_track_r_corr2.push_back( sqrt(_ep_x_corr2*_ep_x_corr2+_ep_y_corr2*_ep_y_corr2) );
    vec_track_r_corr2.push_back( sqrt(_em_x_corr2*_em_x_corr2+_em_y_corr2*_em_y_corr2) );
    vec_truth_x.push_back(_true_gamma_decay_x->at(0));
    vec_truth_y.push_back(_true_gamma_decay_y->at(0));
    vec_truth_z.push_back(_true_gamma_decay_z->at(0));
    vec_truth_r.push_back( sqrt(_true_gamma_decay_x->at(0)*_true_gamma_decay_x->at(0)+_true_gamma_decay_y->at(0)*_true_gamma_decay_y->at(0)) );
  }

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

  TGraph *gr_track_xy = new TGraph(vec_track_x.size(), vec_track_x.data(), vec_track_y.data());
  gr_track_xy->SetMarkerColor(kRed);
  gr_track_xy->SetMarkerStyle(23);
  gr_track_xy->Draw("P,same");

  TGraph *gr_track_corr2_xy = new TGraph(vec_track_x_corr2.size(), vec_track_x_corr2.data(), vec_track_y_corr2.data());
  gr_track_corr2_xy->SetMarkerColor(kGreen);
  gr_track_corr2_xy->SetMarkerStyle(23);
  gr_track_corr2_xy->Draw("P,same");

  TGraph *gr_truth_xy = new TGraph(vec_truth_x.size(), vec_truth_x.data(), vec_truth_y.data());
  gr_truth_xy->SetMarkerColor(kOrange);
  gr_truth_xy->SetMarkerStyle(24);
  gr_truth_xy->Draw("P,same");

  // Draw ellipses using TGraph
  DrawEllipseWithTGraph(0., 0., emcal_radius - 3., kBlack);
  DrawEllipseWithTGraph(0., 0., emcal_radius + 3., kBlack);

  can_xy->Update();
  can_xy->SaveAs(Form("figure/EvtDisplay_xy_run%d_event%d.pdf",runnumber,eventid));

  TCanvas *can_rz = new TCanvas(Form("can_rz"), "can_rz", 800, 800);
  can_rz->SetLeftMargin(0.15);
  can_rz->SetRightMargin(0.05);
  can_rz->cd();

  TGraph *gr_clus_rz = new TGraph(vec_cluster_z.size(), vec_cluster_z.data(), vec_cluster_r.data());
  gr_clus_rz->SetTitle(Form("Run %d Event %d;Z [cm];R [cm]",runnumber,eventid));
  gr_clus_rz->GetXaxis()->SetLimits(-(emcal_radius+30.), (emcal_radius+30.));
  gr_clus_rz->SetMinimum(5);
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

  TGraph *gr_track_rz = new TGraph(vec_track_z.size(), vec_track_z.data(), vec_track_r.data());
  gr_track_rz->SetMarkerColor(kRed);
  gr_track_rz->SetMarkerStyle(23);
  gr_track_rz->Draw("P,same");

  TGraph *gr_track_corr2_rz = new TGraph(vec_track_z_corr2.size(), vec_track_z_corr2.data(), vec_track_r_corr2.data());
  gr_track_corr2_rz->SetMarkerColor(kGreen);
  gr_track_corr2_rz->SetMarkerStyle(23);
  gr_track_corr2_rz->Draw("P,same");

  TGraph *gr_truth_rz = new TGraph(vec_truth_z.size(), vec_truth_z.data(), vec_truth_r.data());
  gr_truth_rz->SetMarkerColor(kOrange);
  gr_truth_rz->SetMarkerStyle(24);
  gr_truth_rz->Draw("P,same");

  // Draw ellipses using TGraph
  DrawHLineWithTGraph(emcal_radius - 3., -150, 150, kBlack);
  DrawHLineWithTGraph(emcal_radius + 3., -150, 150, kBlack);

  can_rz->Update();
  can_rz->SaveAs(Form("figure/EvtDisplay_rz_run%d_event%d.pdf",runnumber,eventid));


}
