#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void EMCalResidual(int runnumber)
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TString inputfile = Form("%d/final_%d_ana.root",runnumber,runnumber);

  //TFile *file = new TFile(inputfile, "READ");
  //TTree *tree = (TTree*)file->Get("tree");

  TChain* chain = new TChain("tree");
  chain->Add(Form("./%d/TrackCalo_*_ana.root",runnumber));

  setBranch(chain);

  std::vector<int> matched_eventid;
  std::vector<std::vector<float>> matched_track_emcal_dphi;
  std::vector<std::vector<float>> matched_track_emcal_dz;
  std::vector<std::vector<float>> matched_track_emcal_deta;
  std::vector<std::vector<int>> matched_track_charge;
  std::vector<std::vector<float>> matched_emcal_y;
  std::vector<std::vector<float>> matched_emcal_z;
  std::vector<std::vector<float>> matched_emcal_phi;
  std::vector<std::vector<float>> matched_emcal_e;
  std::vector<std::vector<float>> matched_track_y;
  std::vector<std::vector<float>> matched_track_z;
  std::vector<std::vector<float>> matched_track_phi;
  std::vector<std::vector<float>> matched_track_p;
  std::vector<std::vector<float>> matched_track_eta;
  std::vector<std::vector<float>> matched_track_z_origin;
  std::vector<std::vector<float>> matched_track_phi_tilt;
  std::vector<float> vertex_z;
  std::vector<float> mbd_z;
  std::vector<std::vector<float>> matched_track_hcal_dphi;
  std::vector<std::vector<float>> matched_track_hcal_dz;
  std::vector<std::vector<float>> matched_track_hcal_deta;
  std::vector<std::vector<float>> matched_hcal_y;
  std::vector<std::vector<float>> matched_hcal_z;
  std::vector<std::vector<float>> matched_hcal_phi;
  std::vector<std::vector<float>> matched_hcal_e;

  matched_eventid.clear();
  matched_track_emcal_dphi.clear();
  matched_track_emcal_dz.clear();
  matched_track_emcal_deta.clear();
  matched_track_charge.clear();
  matched_emcal_y.clear();
  matched_emcal_z.clear();
  matched_emcal_e.clear();
  matched_emcal_phi.clear();
  matched_track_y.clear();
  matched_track_z.clear();
  matched_track_phi.clear();
  matched_track_p.clear();
  matched_track_eta.clear();
  matched_track_z_origin.clear();
  matched_track_phi_tilt.clear();
  vertex_z.clear();
  mbd_z.clear();
  matched_track_hcal_dphi.clear();
  matched_track_hcal_dz.clear();
  matched_track_hcal_deta.clear();
  matched_hcal_y.clear();
  matched_hcal_z.clear();
  matched_hcal_e.clear();
  matched_hcal_phi.clear();

  double ntrack=0;
  double nemcal=0;
  double nhcal=0;
  for(int i = 0; i < chain->GetEntries(); i++)
  {
    chain->GetEntry(i);

    if(_ntracks->at(0) > 800)
    {
      std::cout << "Event: " << i << " rejected by number of tracks: " << _ntracks->at(0) << std::endl;
      continue;
    }

    matched_eventid.push_back(i);
    for (int iv = 0; iv < _vertex_z->size(); iv++)
    {
      if (isnan(_vertex_z->at(iv))) continue;
      vertex_z.push_back(_vertex_z->at(iv));
    }
    for (int iv = 0; iv < _mbd_z->size(); iv++)
    {
      if (isnan(_mbd_z->at(iv))) continue;
      mbd_z.push_back(_mbd_z->at(iv));
    }

    std::vector<float> vec_track_emcal_residual_phi;
    std::vector<float> vec_track_emcal_residual_z;
    std::vector<float> vec_track_emcal_residual_eta;
    std::vector<int> vec_track_charge;
    std::vector<float> vec_emcal_y;
    std::vector<float> vec_emcal_z;
    std::vector<float> vec_emcal_phi;
    std::vector<float> vec_emcal_e;
    std::vector<float> vec_track_y;
    std::vector<float> vec_track_z;
    std::vector<float> vec_track_phi;
    std::vector<float> vec_track_p;
    std::vector<float> vec_track_eta;
    std::vector<float> vec_track_z_origin;
    std::vector<float> vec_track_phi_tilt;
    std::vector<float> vec_track_hcal_residual_phi;
    std::vector<float> vec_track_hcal_residual_z;
    std::vector<float> vec_track_hcal_residual_eta;
    std::vector<float> vec_hcal_y;
    std::vector<float> vec_hcal_z;
    std::vector<float> vec_hcal_phi;
    std::vector<float> vec_hcal_e;
    vec_track_emcal_residual_phi.clear();
    vec_track_emcal_residual_z.clear();
    vec_track_emcal_residual_eta.clear();
    vec_track_charge.clear();
    vec_emcal_y.clear();
    vec_emcal_z.clear();
    vec_emcal_phi.clear();
    vec_emcal_e.clear();
    vec_track_y.clear();
    vec_track_z.clear();
    vec_track_phi.clear();
    vec_track_p.clear();
    vec_track_eta.clear();
    vec_track_z_origin.clear();
    vec_track_phi_tilt.clear();
    vec_track_hcal_residual_phi.clear();
    vec_track_hcal_residual_z.clear();
    vec_track_hcal_residual_eta.clear();
    vec_hcal_y.clear();
    vec_hcal_z.clear();
    vec_hcal_phi.clear();
    vec_hcal_e.clear();

    // loop over all tracks
    for(unsigned int itrack = 0; itrack < _track_ptq->size(); itrack++)
    {
      if(isnan(_track_ptq->at(itrack))) continue;
      // cut good tpc tracks
      if(_track_nc_tpc->at(itrack) < 22) continue;
      if(_track_quality->at(itrack) > 100) continue;

      // we only need emcal matching, hcal is not necessary
      if(isnan(_track_phi_emc->at(itrack)))
      {
        continue;
      }

      // project tpc tracks to emcal and hcal
      std::pair<float, float> TrackProjsEMCal;
      std::pair<float, float> TrackProjsHCal;
      TrackProjsEMCal = std::make_pair(_track_phi_emc->at(itrack), _track_z_emc->at(itrack));
      TrackProjsHCal = std::make_pair(_track_phi_ohc->at(itrack), _track_z_ohc->at(itrack));
      TVector3 p3_track(_track_px->at(itrack),_track_py->at(itrack),_track_pz->at(itrack));
      TVector3 R3_track_emc(_track_x_emc->at(itrack),_track_y_emc->at(itrack),0);

      float track_p = sqrt(pow(_track_px->at(itrack),2) + pow(_track_py->at(itrack),2) + pow(_track_pz->at(itrack),2));
      float track_eta = _track_eta->at(itrack);

      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
        // cut good emcal clusters
        //if(_emcal_e->at(iem) < min_EMCal_E) continue;
        //if(fabs(_emcal_eta->at(iem)) > 1.1) continue;
        std::pair<float, float> EMCalPos;
        float emcal_phi = atan2(_emcal_y->at(iem), _emcal_x->at(iem));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        EMCalPos = std::make_pair(emcal_phi, radius_scale * _emcal_z->at(iem));

        float dphi = TrackProjsEMCal.first - EMCalPos.first;
        float dz = TrackProjsEMCal.second - EMCalPos.second;
        float deta = _track_eta_emc->at(itrack) - _emcal_eta->at(iem);

        vec_track_emcal_residual_phi.push_back(PiRange(dphi));
        vec_track_emcal_residual_z.push_back(dz);
        vec_track_emcal_residual_eta.push_back(deta);
        vec_emcal_phi.push_back(EMCalPos.first);
        vec_emcal_y.push_back(_emcal_y->at(iem) * radius_scale);
        vec_emcal_z.push_back(EMCalPos.second);
        vec_emcal_e.push_back(_emcal_e->at(iem));
        vec_track_phi.push_back(TrackProjsEMCal.first);
        vec_track_y.push_back(_track_y_emc->at(itrack));
        vec_track_z.push_back(TrackProjsEMCal.second);
        vec_track_p.push_back(track_p);
        vec_track_eta.push_back(track_eta);
        vec_track_z_origin.push_back(_track_z_origin->at(itrack));
        vec_track_phi_tilt.push_back(p3_track.Dot(R3_track_emc.Cross(z_direction)) / p3_track.Dot(R3_track_emc));
        if (_track_ptq->at(itrack)>0)
        {
          vec_track_charge.push_back(1);
        }
        else if (_track_ptq->at(itrack)<0)
        {
          vec_track_charge.push_back(-1);
        }
        else
        {
          cout<<"what??? ptq = "<<_track_ptq->at(itrack)<<endl;
        }
      }

      // loop all hcal clusters to match with tpc tracks
      for(unsigned int ihad = 0; ihad < _hcal_e->size(); ihad++)
      {
        // cut good hcal clusters
        //if(_hcal_e->at(ihad) < min_HCal_E) continue;
        //if(fabs(_hcal_eta->at(ihad)) > 1.1) continue;
        std::pair<float, float> HCalPos;
        HCalPos = std::make_pair(_hcal_phi->at(ihad), _hcal_z->at(ihad));

        float dphi = TrackProjsHCal.first - HCalPos.first;
        float dz = TrackProjsHCal.second - HCalPos.second;
        float deta = _track_eta_ohc->at(itrack) - _hcal_eta->at(ihad);

        vec_track_hcal_residual_phi.push_back(PiRange(dphi));
        vec_track_hcal_residual_z.push_back(dz);
        vec_track_hcal_residual_eta.push_back(deta);
        vec_hcal_phi.push_back(HCalPos.first);
        vec_hcal_y.push_back(_hcal_y->at(ihad));
        vec_hcal_z.push_back(HCalPos.second);
        vec_hcal_e.push_back(_hcal_e->at(ihad));
      }
    }

    matched_track_emcal_dphi.push_back(vec_track_emcal_residual_phi);
    matched_track_emcal_dz.push_back(vec_track_emcal_residual_z);
    matched_track_emcal_deta.push_back(vec_track_emcal_residual_eta);
    matched_track_charge.push_back(vec_track_charge);
    matched_emcal_y.push_back(vec_emcal_y);
    matched_emcal_z.push_back(vec_emcal_z);
    matched_emcal_phi.push_back(vec_emcal_phi);
    matched_emcal_e.push_back(vec_emcal_e);
    matched_track_y.push_back(vec_track_y);
    matched_track_z.push_back(vec_track_z);
    matched_track_phi.push_back(vec_track_phi);
    matched_track_p.push_back(vec_track_p);
    matched_track_eta.push_back(vec_track_eta);
    matched_track_z_origin.push_back(vec_track_z_origin);
    matched_track_phi_tilt.push_back(vec_track_phi_tilt);
    matched_track_hcal_dphi.push_back(vec_track_hcal_residual_phi);
    matched_track_hcal_dz.push_back(vec_track_hcal_residual_z);
    matched_track_hcal_deta.push_back(vec_track_hcal_residual_eta);
    matched_hcal_y.push_back(vec_hcal_y);
    matched_hcal_z.push_back(vec_hcal_z);
    matched_hcal_phi.push_back(vec_hcal_phi);
    matched_hcal_e.push_back(vec_hcal_e);

    ntrack += _track_ptq->size();
    nemcal += _emcal_e->size();
    nhcal += _hcal_e->size();

  }
  ntrack /= chain->GetEntries();
  nemcal /= chain->GetEntries();
  nhcal /= chain->GetEntries();
  cout<<"average track per event = "<<ntrack<<endl;
  cout<<"average emcal per event = "<<nemcal<<endl;
  cout<<"average hcal per event = "<<nhcal<<endl;

  cout<<"number of matched event = "<<matched_track_charge.size()<<endl;

  TH2* h2_dphi_dz_track_emcal = new TH2F("h2_dphi_dz_track_emcal", "h2_dphi_dz_track_emcal", 100, -0.1, 0.1, 100, -20, 20);
  h2_dphi_dz_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h2_dphi_dz_track_emcal->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_dz_track_emcal->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_dphi_dz_track_emcal->GetZaxis()->SetTitle("Entries");

  //TH2* h2_dphi_dz_track_emcal_pos = new TH2F("h2_dphi_dz_track_emcal_pos", "h2_dphi_dz_track_emcal_pos", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_dphi_dz_track_emcal_pos = new TH2F("h2_dphi_dz_track_emcal_pos", "h2_dphi_dz_track_emcal_pos", 100, -0.2, 0.2, 100, -20, 20);
  h2_dphi_dz_track_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h2_dphi_dz_track_emcal_pos->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_dz_track_emcal_pos->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_dphi_dz_track_emcal_pos->GetZaxis()->SetTitle("Entries");
  //h2_dphi_dz_track_emcal_pos->GetZaxis()->SetTitleOffset(1.5);

  //TH2* h2_dphi_dz_track_emcal_neg = new TH2F("h2_dphi_dz_track_emcal_neg", "h2_dphi_dz_track_emcal_neg", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_dphi_dz_track_emcal_neg = new TH2F("h2_dphi_dz_track_emcal_neg", "h2_dphi_dz_track_emcal_neg", 100, -0.2, 0.2, 100, -20, 20);
  h2_dphi_dz_track_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h2_dphi_dz_track_emcal_neg->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_dz_track_emcal_neg->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_dphi_dz_track_emcal_neg->GetZaxis()->SetTitle("Entries");
  //h2_dphi_dz_track_emcal_neg->GetZaxis()->SetTitleOffset(1.5);

  TH2* h2_dphi_deta_track_emcal_pos = new TH2F("h2_dphi_deta_track_emcal_pos", "h2_dphi_deta_track_emcal_pos", 100, -M_PI, M_PI, 100, -3, 3);
  h2_dphi_deta_track_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h2_dphi_deta_track_emcal_pos->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_deta_track_emcal_pos->GetYaxis()->SetTitle("#Delta#eta");
  h2_dphi_deta_track_emcal_pos->GetZaxis()->SetTitle("Entries");

  TH2* h2_dphi_deta_track_emcal_neg = new TH2F("h2_dphi_deta_track_emcal_neg", "h2_dphi_deta_track_emcal_neg", 100, -M_PI, M_PI, 100, -3, 3);
  h2_dphi_deta_track_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h2_dphi_deta_track_emcal_neg->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_deta_track_emcal_neg->GetYaxis()->SetTitle("#Delta#eta");
  h2_dphi_deta_track_emcal_neg->GetZaxis()->SetTitle("Entries");

  TH1* h1_dphi_track_emcal = new TH1F("h1_dphi_track_emcal", "h1_dphi_track_emcal", 100, -.1, .1);
  h1_dphi_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_dphi_track_emcal->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_emcal->GetYaxis()->SetTitle(Form("Events / %.3f rad",0.2/100.));
  h1_dphi_track_emcal->SetMinimum(0);

  TH1* h1_dphi_track_emcal_pos = new TH1F("h1_dphi_track_emcal_pos", "h1_dphi_track_emcal_pos", 100, -.1, .1);
  h1_dphi_track_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h1_dphi_track_emcal_pos->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_emcal_pos->GetYaxis()->SetTitle(Form("Events / %.3f rad",0.2/100.));
  h1_dphi_track_emcal_pos->SetMinimum(0);

  TH1* h1_dphi_track_emcal_neg = new TH1F("h1_dphi_track_emcal_neg", "h1_dphi_track_emcal_neg", 100, -.1, .1);
  h1_dphi_track_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h1_dphi_track_emcal_neg->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_emcal_neg->GetYaxis()->SetTitle(Form("Events / %.3f rad",0.2/100.));
  h1_dphi_track_emcal_neg->SetMinimum(0);

  TH1* h1_dz_track_emcal = new TH1F("h1_dz_track_emcal", "h1_dz_track_emcal", 100, -50, 50);
  h1_dz_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_dz_track_emcal->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_emcal->GetYaxis()->SetTitle(Form("Events / %1f cm",100./100.));
  h1_dz_track_emcal->SetMinimum(0);

  TH1* h1_dz_track_emcal_pos = new TH1F("h1_dz_track_emcal_pos", "h1_dz_track_emcal_pos", 100, -50, 50);
  h1_dz_track_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h1_dz_track_emcal_pos->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_emcal_pos->GetYaxis()->SetTitle(Form("Events / %1f cm",100./100.));
  h1_dz_track_emcal_pos->SetMinimum(0);

  TH1* h1_dz_track_emcal_neg = new TH1F("h1_dz_track_emcal_neg", "h1_dz_track_emcal_neg", 100, -50, 50);
  h1_dz_track_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h1_dz_track_emcal_neg->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_emcal_neg->GetYaxis()->SetTitle(Form("Events / %1f cm",100./100.));
  h1_dz_track_emcal_neg->SetMinimum(0);

  TH1* h1_deta_track_emcal_pos = new TH1F("h1_deta_track_emcal_pos", "h1_deta_track_emcal_pos", 100, -5, 5);
  h1_deta_track_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h1_deta_track_emcal_pos->GetXaxis()->SetTitle("#Delta#eta");
  h1_deta_track_emcal_pos->GetYaxis()->SetTitle(Form("Events / %.1f",10./100.));
  h1_deta_track_emcal_pos->SetMinimum(0);

  TH1* h1_deta_track_emcal_neg = new TH1F("h1_deta_track_emcal_neg", "h1_deta_track_emcal_neg", 100, -5, 5);
  h1_deta_track_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h1_deta_track_emcal_neg->GetXaxis()->SetTitle("#Delta#eta");
  h1_deta_track_emcal_neg->GetYaxis()->SetTitle(Form("Events / %.1f",10./100.));
  h1_deta_track_emcal_neg->SetMinimum(0);

  //TH2* h2_z_dz_emcal = new TH2F("h2_z_dz_emcal", "h2_z_dz_emcal", 100, -300, 300, 100, -300, 300);
  //TH2* h2_z_dz_emcal = new TH2F("h2_z_dz_emcal", "h2_z_dz_emcal", 50, -150, 150, 50, -100, 100);
  TH2* h2_z_dz_emcal = new TH2F("h2_z_dz_emcal", "h2_z_dz_emcal", 50, -150, 150, 40, -20, 20);
  h2_z_dz_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_dz_emcal->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_dz_emcal->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_z_dz_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dz_emcal->GetZaxis()->SetTitle("Entries");
  h2_z_dz_emcal->GetZaxis()->SetTitleOffset(1);

  //TH2* h2_z_dz_emcal_pos = new TH2F("h2_z_dz_emcal_pos", "h2_z_dz_emcal_pos", 100, -300, 300, 100, -300, 300);
  TH2* h2_z_dz_emcal_pos = new TH2F("h2_z_dz_emcal_pos", "h2_z_dz_emcal_pos", 50, -150, 150, 50, -100, 100);
  h2_z_dz_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive Only",runnumber));
  h2_z_dz_emcal_pos->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_dz_emcal_pos->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_z_dz_emcal_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dz_emcal_pos->GetZaxis()->SetTitle("Entries");
  h2_z_dz_emcal_pos->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dz_emcal_neg = new TH2F("h2_z_dz_emcal_neg", "h2_z_dz_emcal_neg", 100, -300, 300, 100, -300, 300);
  TH2* h2_z_dz_emcal_neg = new TH2F("h2_z_dz_emcal_neg", "h2_z_dz_emcal_neg", 50, -150, 150, 50, -100, 100);
  h2_z_dz_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative Only",runnumber));
  h2_z_dz_emcal_neg->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_dz_emcal_neg->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_z_dz_emcal_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dz_emcal_neg->GetZaxis()->SetTitle("Entries");
  h2_z_dz_emcal_neg->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dz_track = new TH2F("h2_z_dz_track", "h2_z_dz_track", 100, -300, 300, 100, -300, 300);
  TH2* h2_z_dz_track = new TH2F("h2_z_dz_track", "h2_z_dz_track", 50, -150, 150, 50, -100, 100);
  h2_z_dz_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_dz_track->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_dz_track->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_z_dz_track->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dz_track->GetZaxis()->SetTitle("Entries");
  h2_z_dz_track->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dz_track_pos = new TH2F("h2_z_dz_track_pos", "h2_z_dz_track_pos", 100, -300, 300, 100, -300, 300);
  TH2* h2_z_dz_track_pos = new TH2F("h2_z_dz_track_pos", "h2_z_dz_track_pos", 50, -150, 150, 50, -100, 100);
  h2_z_dz_track_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive Only",runnumber));
  h2_z_dz_track_pos->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_dz_track_pos->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_z_dz_track_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dz_track_pos->GetZaxis()->SetTitle("Entries");
  h2_z_dz_track_pos->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dz_track_neg = new TH2F("h2_z_dz_track_neg", "h2_z_dz_track_neg", 100, -300, 300, 100, -300, 300);
  TH2* h2_z_dz_track_neg = new TH2F("h2_z_dz_track_neg", "h2_z_dz_track_neg", 50, -150, 150, 100, -100, 100);
  h2_z_dz_track_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative Only",runnumber));
  h2_z_dz_track_neg->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_dz_track_neg->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_z_dz_track_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dz_track_neg->GetZaxis()->SetTitle("Entries");
  h2_z_dz_track_neg->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dphi_emcal = new TH2F("h2_z_dphi_emcal", "h2_z_dphi_emcal", 100, -300, 300, 100, -M_PI, M_PI);
  TH2* h2_z_dphi_emcal = new TH2F("h2_z_dphi_emcal", "h2_z_dphi_emcal", 50, -150, 150, 20, -0.2, 0.2);
  h2_z_dphi_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_dphi_emcal->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_dphi_emcal->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_z_dphi_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dphi_emcal->GetZaxis()->SetTitle("Entries");
  h2_z_dphi_emcal->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dphi_emcal_pos = new TH2F("h2_z_dphi_emcal_pos", "h2_z_dphi_emcal_pos", 100, -300, 300, 100, -M_PI, M_PI);
  TH2* h2_z_dphi_emcal_pos = new TH2F("h2_z_dphi_emcal_pos", "h2_z_dphi_emcal_pos", 50, -150, 150, 20, -0.2, 0.2);
  h2_z_dphi_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive Only",runnumber));
  h2_z_dphi_emcal_pos->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_dphi_emcal_pos->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_z_dphi_emcal_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dphi_emcal_pos->GetZaxis()->SetTitle("Entries");
  h2_z_dphi_emcal_pos->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dphi_emcal_neg = new TH2F("h2_z_dphi_emcal_neg", "h2_z_dphi_emcal_neg", 100, -300, 300, 100, -M_PI, M_PI);
  TH2* h2_z_dphi_emcal_neg = new TH2F("h2_z_dphi_emcal_neg", "h2_z_dphi_emcal_neg", 50, -150, 150, 20, -0.2, 0.2);
  h2_z_dphi_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative Only",runnumber));
  h2_z_dphi_emcal_neg->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_dphi_emcal_neg->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_z_dphi_emcal_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dphi_emcal_neg->GetZaxis()->SetTitle("Entries");
  h2_z_dphi_emcal_neg->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dphi_track = new TH2F("h2_z_dphi_track", "h2_z_dphi_track", 100, -300, 300, 100, -M_PI, M_PI);
  TH2* h2_z_dphi_track = new TH2F("h2_z_dphi_track", "h2_z_dphi_track", 50, -150, 150, 20, -0.2, 0.2);
  h2_z_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_dphi_track->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_z_dphi_track->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dphi_track->GetZaxis()->SetTitle("Entries");
  h2_z_dphi_track->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dphi_track_pos = new TH2F("h2_z_dphi_track_pos", "h2_z_dphi_track_pos", 100, -300, 300, 100, -M_PI, M_PI);
  TH2* h2_z_dphi_track_pos = new TH2F("h2_z_dphi_track_pos", "h2_z_dphi_track_pos", 50, -150, 150, 20, -0.2, 0.2);
  h2_z_dphi_track_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive Only",runnumber));
  h2_z_dphi_track_pos->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_dphi_track_pos->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_z_dphi_track_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dphi_track_pos->GetZaxis()->SetTitle("Entries");
  h2_z_dphi_track_pos->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_z_dphi_track_neg = new TH2F("h2_z_dphi_track_neg", "h2_z_dphi_track_neg", 100, -300, 300, 100, -M_PI, M_PI);
  TH2* h2_z_dphi_track_neg = new TH2F("h2_z_dphi_track_neg", "h2_z_dphi_track_neg", 50, -150, 150, 20, -0.2, 0.2);
  h2_z_dphi_track_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative Only",runnumber));
  h2_z_dphi_track_neg->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_dphi_track_neg->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_z_dphi_track_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_z_dphi_track_neg->GetZaxis()->SetTitle("Entries");
  h2_z_dphi_track_neg->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dphi_emcal = new TH2F("h2_phi_dphi_emcal", "h2_phi_dphi_emcal", 100, -M_PI, M_PI, 100, -M_PI, M_PI);
  TH2* h2_phi_dphi_emcal = new TH2F("h2_phi_dphi_emcal", "h2_phi_dphi_emcal", 50, -M_PI, M_PI, 50, -0.2, 0.2);
  h2_phi_dphi_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_phi_dphi_emcal->GetXaxis()->SetTitle("Calo #Phi [rad]");
  h2_phi_dphi_emcal->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_phi_dphi_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dphi_emcal->GetZaxis()->SetTitle("Entries");
  h2_phi_dphi_emcal->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dphi_track = new TH2F("h2_phi_dphi_track", "h2_phi_dphi_track", 100, -M_PI, M_PI, 100, -M_PI, M_PI);
  TH2* h2_phi_dphi_track = new TH2F("h2_phi_dphi_track", "h2_phi_dphi_track", 50, -M_PI, M_PI, 50, -0.2, 0.2);
  h2_phi_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_phi_dphi_track->GetXaxis()->SetTitle("Track #Phi [rad]");
  h2_phi_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_phi_dphi_track->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dphi_track->GetZaxis()->SetTitle("Entries");
  h2_phi_dphi_track->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_tilt_dphi_track = new TH2F("h2_phi_tilt_dphi_track", "h2_phi_tilt_dphi_track", 50, -1, 1, 50, -0.0, 0.05);
  TH2* h2_phi_tilt_dphi_track = new TH2F("h2_phi_tilt_dphi_track", "h2_phi_tilt_dphi_track", 50, -1, 1, 50, -0.2, 0.2);
  h2_phi_tilt_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_phi_tilt_dphi_track->GetXaxis()->SetTitle("\\vec{p}\\cdot (\\vec{R}\\times\\vec{z})/\\vec{p}\\cdot \\vec{R}");
  h2_phi_tilt_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_phi_tilt_dphi_track->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_tilt_dphi_track->GetZaxis()->SetTitle("Entries");
  h2_phi_tilt_dphi_track->GetZaxis()->SetTitleOffset(1);

  TH2* h2_phi_tilt_dphi_track_pos = new TH2F("h2_phi_tilt_dphi_track_pos", "h2_phi_tilt_dphi_track_pos", 50, -1, 1, 50, -0.2, 0.2);
  h2_phi_tilt_dphi_track_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive only",runnumber));
  h2_phi_tilt_dphi_track_pos->GetXaxis()->SetTitle("#vec{p}\\cdot (#vec{R}#times#vec{z})/#vec{p}\\cdot #vec{R}");
  h2_phi_tilt_dphi_track_pos->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_phi_tilt_dphi_track_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_tilt_dphi_track_pos->GetZaxis()->SetTitle("Entries");
  h2_phi_tilt_dphi_track_pos->GetZaxis()->SetTitleOffset(1.2);

  TH2* h2_phi_tilt_dphi_track_neg = new TH2F("h2_phi_tilt_dphi_track_neg", "h2_phi_tilt_dphi_track_neg", 50, -1, 1, 50, -0.2, 0.2);
  h2_phi_tilt_dphi_track_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative only",runnumber));
  h2_phi_tilt_dphi_track_neg->GetXaxis()->SetTitle("#vec{p}\\cdot (#vec{R}#times#vec{z})/#vec{p}\\cdot #vec{R}");
  h2_phi_tilt_dphi_track_neg->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_phi_tilt_dphi_track_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_tilt_dphi_track_neg->GetZaxis()->SetTitle("Entries");
  h2_phi_tilt_dphi_track_neg->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dz_emcal = new TH2F("h2_phi_dz_emcal", "h2_phi_dz_emcal", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_phi_dz_emcal = new TH2F("h2_phi_dz_emcal", "h2_phi_dz_emcal", 50, -M_PI, M_PI, 50, -100, 100);
  h2_phi_dz_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_phi_dz_emcal->GetXaxis()->SetTitle("Calo #Phi [rad]");
  h2_phi_dz_emcal->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_phi_dz_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dz_emcal->GetZaxis()->SetTitle("Entries");
  h2_phi_dz_emcal->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dz_emcal_pos = new TH2F("h2_phi_dz_emcal_pos", "h2_phi_dz_emcal_pos", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_phi_dz_emcal_pos = new TH2F("h2_phi_dz_emcal_pos", "h2_phi_dz_emcal_pos", 50, -M_PI, M_PI, 50, -100, 100);
  h2_phi_dz_emcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h2_phi_dz_emcal_pos->GetXaxis()->SetTitle("Calo #Phi [rad]");
  h2_phi_dz_emcal_pos->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_phi_dz_emcal_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dz_emcal_pos->GetZaxis()->SetTitle("Entries");
  h2_phi_dz_emcal_pos->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dz_emcal_neg = new TH2F("h2_phi_dz_emcal_neg", "h2_phi_dz_emcal_neg", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_phi_dz_emcal_neg = new TH2F("h2_phi_dz_emcal_neg", "h2_phi_dz_emcal_neg", 50, -M_PI, M_PI, 50, -100, 100);
  h2_phi_dz_emcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h2_phi_dz_emcal_neg->GetXaxis()->SetTitle("Calo #Phi [rad]");
  h2_phi_dz_emcal_neg->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_phi_dz_emcal_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dz_emcal_neg->GetZaxis()->SetTitle("Entries");
  h2_phi_dz_emcal_neg->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dz_track = new TH2F("h2_phi_dz_track", "h2_phi_dz_track", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_phi_dz_track = new TH2F("h2_phi_dz_track", "h2_phi_dz_track", 50, -M_PI, M_PI, 50, -100, 100);
  h2_phi_dz_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_phi_dz_track->GetXaxis()->SetTitle("Track #Phi [rad]");
  h2_phi_dz_track->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_phi_dz_track->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dz_track->GetZaxis()->SetTitle("Entries");
  h2_phi_dz_track->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dz_track_pos = new TH2F("h2_phi_dz_track_pos", "h2_phi_dz_track_pos", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_phi_dz_track_pos = new TH2F("h2_phi_dz_track_pos", "h2_phi_dz_track_pos", 50, -M_PI, M_PI, 50, -100, 100);
  h2_phi_dz_track_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h2_phi_dz_track_pos->GetXaxis()->SetTitle("Track #Phi [rad]");
  h2_phi_dz_track_pos->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_phi_dz_track_pos->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dz_track_pos->GetZaxis()->SetTitle("Entries");
  h2_phi_dz_track_pos->GetZaxis()->SetTitleOffset(1.2);

  //TH2* h2_phi_dz_track_neg = new TH2F("h2_phi_dz_track_neg", "h2_phi_dz_track_neg", 100, -M_PI, M_PI, 100, -300, 300);
  TH2* h2_phi_dz_track_neg = new TH2F("h2_phi_dz_track_neg", "h2_phi_dz_track_neg", 50, -M_PI, M_PI, 50, -100, 100);
  h2_phi_dz_track_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h2_phi_dz_track_neg->GetXaxis()->SetTitle("Track #Phi [rad]");
  h2_phi_dz_track_neg->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_phi_dz_track_neg->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_dz_track_neg->GetZaxis()->SetTitle("Entries");
  h2_phi_dz_track_neg->GetZaxis()->SetTitleOffset(1.2);

  TH2* h2_phi_track_emcal = new TH2F("h2_phi_track_emcal", "h2_phi_track_emcal", 100, -M_PI, M_PI, 100, -M_PI, M_PI);
  h2_phi_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_phi_track_emcal->GetXaxis()->SetTitle("Track #Phi [rad]");
  h2_phi_track_emcal->GetYaxis()->SetTitle("Calo #Phi [rad]");
  h2_phi_track_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_phi_track_emcal->GetZaxis()->SetTitle("Entries");
  h2_phi_track_emcal->GetZaxis()->SetTitleOffset(1.2);

  TH2* h2_z_track_emcal = new TH2F("h2_z_track_emcal", "h2_z_track_emcal", 100, -180, 180, 100, -140, 140);
  h2_z_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_track_emcal->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_track_emcal->GetYaxis()->SetTitle("Calo Z [cm]");
  h2_z_track_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_z_track_emcal->GetZaxis()->SetTitle("Entries");
  h2_z_track_emcal->GetZaxis()->SetTitleOffset(1.2);

  TH1* h1_phi_track = new TH1F("h1_phi_track", "h1_phi_track", 100, -M_PI, M_PI);
  h1_phi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_phi_track->GetXaxis()->SetTitle("#Phi [rad]");
  h1_phi_track->GetYaxis()->SetTitle(Form("Events / %.2f rad",(2*M_PI)/100.));
  h1_phi_track->SetMinimum(0);

  TH2* h2_eOp_p = new TH2F("h2_eOp_p", "h2_eOp_p", 100, 0, 15, 100, 0, 2);
  h2_eOp_p->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_eOp_p->GetXaxis()->SetTitle("Track p [GeV/c]");
  h2_eOp_p->GetYaxis()->SetTitle("E/p");
  h2_eOp_p->GetYaxis()->SetTitleOffset(1.2);
  h2_eOp_p->GetZaxis()->SetTitle("Entries");
  h2_eOp_p->GetZaxis()->SetTitleOffset(1.2);

  TH1* h1_eOp = new TH1F("h1_eOp", "h1_eOp", 100, 0, 2);
  h1_eOp->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_eOp->GetXaxis()->SetTitle("E/p");
  h1_eOp->GetYaxis()->SetTitle(Form("Events / %.2f",(2.-0.)/100.));
  h1_eOp->SetMinimum(0);

  TH2* h2_z_y_emcal = new TH2F("h2_z_y_emcal", "h2_z_y_emcal", 100, -300, 300, 100, -300, 300);
  h2_z_y_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_y_emcal->GetXaxis()->SetTitle("Calo Z [cm]");
  h2_z_y_emcal->GetYaxis()->SetTitle("Calo Y [cm]");
  h2_z_y_emcal->GetYaxis()->SetTitleOffset(1.2);
  h2_z_y_emcal->GetZaxis()->SetTitle("Entries");
  h2_z_y_emcal->GetZaxis()->SetTitleOffset(1.2);

  TH2* h2_z_y_track = new TH2F("h2_z_y_track", "h2_z_y_track", 100, -300, 300, 100, -300, 300);
  h2_z_y_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_z_y_track->GetXaxis()->SetTitle("Track Z [cm]");
  h2_z_y_track->GetYaxis()->SetTitle("Track Y [cm]");
  h2_z_y_track->GetYaxis()->SetTitleOffset(1.2);
  h2_z_y_track->GetZaxis()->SetTitle("Entries");
  h2_z_y_track->GetZaxis()->SetTitleOffset(1.2);

  TH1* h1_track_phi_tilt = new TH1F("h1_track_phi_tilt", "h1_track_phi_tilt", 100, -1, 1);
  h1_track_phi_tilt->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_track_phi_tilt->GetXaxis()->SetTitle("#vec{p}\\cdot (#vec{R}#times#vec{z})/#vec{p}\\cdot #vec{R}");
  h1_track_phi_tilt->GetYaxis()->SetTitle(Form("Events / %.2f",2./100.));
  h1_track_phi_tilt->SetMinimum(0);

  TH1* h1_z_origin_track = new TH1F("h1_z_origin_track", "h1_z_origin_track", 100, -150, 150);
  h1_z_origin_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_z_origin_track->GetXaxis()->SetTitle("Track Z @ R=0 [cm]");
  h1_z_origin_track->GetYaxis()->SetTitle(Form("Events / %1f cm",300./100.));
  h1_z_origin_track->SetMinimum(0);

  TH1* h1_z_vertex = new TH1F("h1_z_vertex", "h1_z_vertex", 100, -300, 300);
  h1_z_vertex->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_z_vertex->GetXaxis()->SetTitle("Vertex Z [cm]");
  h1_z_vertex->GetYaxis()->SetTitle(Form("Events / %.1f cm",400./100.));
  h1_z_vertex->SetMinimum(0);

  TH1* h1_z_mbd = new TH1F("h1_z_mbd", "h1_z_mbd", 100, -300, 300);
  h1_z_mbd->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_z_mbd->GetXaxis()->SetTitle("MBD z-vertex [cm]");
  h1_z_mbd->GetYaxis()->SetTitle(Form("Events / %.1f cm",400./100.));
  h1_z_mbd->SetMinimum(0);

  TH1* h1_e_emcal = new TH1F("h1_e_emcal", "h1_e_emcal", 100, 0, 5);
  h1_e_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_e_emcal->GetXaxis()->SetTitle("Calo E [GeV]");
  h1_e_emcal->GetYaxis()->SetTitle(Form("Events / %.2f GeV",5./100.));
  h1_e_emcal->SetMinimum(0);

  TH2* h2_dphi_dz_track_hcal_pos = new TH2F("h2_dphi_dz_track_hcal_pos", "h2_dphi_dz_track_hcal_pos", 100, -M_PI, M_PI, 100, -300, 300);
  h2_dphi_dz_track_hcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h2_dphi_dz_track_hcal_pos->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_dz_track_hcal_pos->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_dphi_dz_track_hcal_pos->GetZaxis()->SetTitle("Entries");

  TH2* h2_dphi_dz_track_hcal_neg = new TH2F("h2_dphi_dz_track_hcal_neg", "h2_dphi_dz_track_hcal_neg", 100, -M_PI, M_PI, 100, -300, 300);
  h2_dphi_dz_track_hcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h2_dphi_dz_track_hcal_neg->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_dz_track_hcal_neg->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_dphi_dz_track_hcal_neg->GetZaxis()->SetTitle("Entries");

  TH1* h1_dphi_track_hcal = new TH1F("h1_dphi_track_hcal", "h1_dphi_track_hcal", 100, -.5, .5);
  h1_dphi_track_hcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_dphi_track_hcal->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_hcal->GetYaxis()->SetTitle(Form("Events / %.1f rad",10./100.));
  h1_dphi_track_hcal->SetMinimum(0);

  TH1* h1_dphi_track_hcal_pos = new TH1F("h1_dphi_track_hcal_pos", "h1_dphi_track_hcal_pos", 100, -.5, .5);
  h1_dphi_track_hcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h1_dphi_track_hcal_pos->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_hcal_pos->GetYaxis()->SetTitle(Form("Events / %.1f rad",10./100.));
  h1_dphi_track_hcal_pos->SetMinimum(0);

  TH1* h1_dphi_track_hcal_neg = new TH1F("h1_dphi_track_hcal_neg", "h1_dphi_track_hcal_neg", 100, -.5, .5);
  h1_dphi_track_hcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h1_dphi_track_hcal_neg->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_hcal_neg->GetYaxis()->SetTitle(Form("Events / %.1f rad",10./100.));
  h1_dphi_track_hcal_neg->SetMinimum(0);

  TH1* h1_dz_track_hcal = new TH1F("h1_dz_track_hcal", "h1_dz_track_hcal", 100, -100, 100);
  h1_dz_track_hcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_dz_track_hcal->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_hcal->GetYaxis()->SetTitle(Form("Events / %1f cm",200./100.));
  h1_dz_track_hcal->SetMinimum(0);

  TH1* h1_dz_track_hcal_pos = new TH1F("h1_dz_track_hcal_pos", "h1_dz_track_hcal_pos", 100, -100, 100);
  h1_dz_track_hcal_pos->SetTitle(Form("sPHENIX Internal, Run %d, Positive charge",runnumber));
  h1_dz_track_hcal_pos->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_hcal_pos->GetYaxis()->SetTitle(Form("Events / %1f cm",200./100.));
  h1_dz_track_hcal_pos->SetMinimum(0);

  TH1* h1_dz_track_hcal_neg = new TH1F("h1_dz_track_hcal_neg", "h1_dz_track_hcal_neg", 100, -100, 100);
  h1_dz_track_hcal_neg->SetTitle(Form("sPHENIX Internal, Run %d, Negative charge",runnumber));
  h1_dz_track_hcal_neg->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_hcal_neg->GetYaxis()->SetTitle(Form("Events / %1f cm",200./100.));
  h1_dz_track_hcal_neg->SetMinimum(0);

  TH2* h2_track_p_dphi_track = new TH2F("h2_track_p_dphi_track", "h2_track_p_dphi_track", 50, 0, 10, 50, -0.2, 0.2);
  h2_track_p_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_track_p_dphi_track->GetXaxis()->SetTitle("p [GeV/#it{c}]");
  h2_track_p_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_track_p_dphi_track->GetYaxis()->SetTitleOffset(1.2);
  h2_track_p_dphi_track->GetZaxis()->SetTitle("Entries");
  h2_track_p_dphi_track->GetZaxis()->SetTitleOffset(1.2);

  TH2* h2_track_eta_dphi_track = new TH2F("h2_track_eta_dphi_track", "h2_track_eta_dphi_track", 50, -2, 2, 50, -0.2, 0.2);
  h2_track_eta_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_track_eta_dphi_track->GetXaxis()->SetTitle("#eta");
  h2_track_eta_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
  h2_track_eta_dphi_track->GetYaxis()->SetTitleOffset(1.2);
  h2_track_eta_dphi_track->GetZaxis()->SetTitle("Entries");
  h2_track_eta_dphi_track->GetZaxis()->SetTitleOffset(1.2);

  for (int index = 0; index<(matched_track_emcal_dphi.size()); index++)
  {
    for (int i = 0; i < (matched_track_emcal_dphi.at(index).size()); i++)
    {

      h2_phi_dphi_emcal->Fill(matched_emcal_phi.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
      h2_phi_dphi_track->Fill(matched_track_phi.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
      h2_phi_track_emcal->Fill(matched_track_phi.at(index).at(i),matched_emcal_phi.at(index).at(i));
      h2_z_track_emcal->Fill(matched_track_z.at(index).at(i),matched_emcal_z.at(index).at(i));

      h2_z_dphi_emcal->Fill(matched_emcal_z.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
      h2_z_dphi_track->Fill(matched_track_z.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
      h2_dphi_dz_track_emcal->Fill(matched_track_emcal_dphi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));

      if (matched_track_charge.at(index).at(i)>0)
      {
        h2_dphi_dz_track_emcal_pos->Fill(matched_track_emcal_dphi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        h2_dphi_deta_track_emcal_pos->Fill(matched_track_emcal_dphi.at(index).at(i),matched_track_emcal_deta.at(index).at(i));
        h1_dphi_track_emcal_pos->Fill(matched_track_emcal_dphi.at(index).at(i));
        h1_dz_track_emcal_pos->Fill(matched_track_emcal_dz.at(index).at(i));
        h1_deta_track_emcal_pos->Fill(matched_track_emcal_deta.at(index).at(i));
        h2_z_dphi_emcal_pos->Fill(matched_emcal_z.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        h2_z_dphi_track_pos->Fill(matched_track_z.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
      }
      else if (matched_track_charge.at(index).at(i)<0)
      {
        h2_dphi_dz_track_emcal_neg->Fill(matched_track_emcal_dphi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        h2_dphi_deta_track_emcal_neg->Fill(matched_track_emcal_dphi.at(index).at(i),matched_track_emcal_deta.at(index).at(i));
        h1_dphi_track_emcal_neg->Fill(matched_track_emcal_dphi.at(index).at(i));
        h1_dz_track_emcal_neg->Fill(matched_track_emcal_dz.at(index).at(i));
        h1_deta_track_emcal_neg->Fill(matched_track_emcal_deta.at(index).at(i));
        h1_deta_track_emcal_neg->Fill(matched_track_emcal_deta.at(index).at(i));
        h2_z_dphi_emcal_neg->Fill(matched_emcal_z.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        h2_z_dphi_track_neg->Fill(matched_track_z.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
      }

      // add dphi cut them draw dz plots
      if (fabs(matched_track_emcal_dphi.at(index).at(i))<0.1)
      {
        h1_dz_track_emcal->Fill(matched_track_emcal_dz.at(index).at(i));
        h2_z_dz_emcal->Fill(matched_emcal_z.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        h2_phi_dz_emcal->Fill(matched_emcal_phi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        h2_z_dz_track->Fill(matched_track_z.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        h2_phi_dz_track->Fill(matched_track_phi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        if (matched_track_charge.at(index).at(i)<0)
        {
          h2_z_dz_emcal_neg->Fill(matched_emcal_z.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
          h2_phi_dz_emcal_neg->Fill(matched_emcal_phi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
          h2_z_dz_track_neg->Fill(matched_track_z.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
          h2_phi_dz_track_neg->Fill(matched_track_phi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        }
        else if (matched_track_charge.at(index).at(i)>0)
        {
          h2_z_dz_emcal_pos->Fill(matched_emcal_z.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
          h2_phi_dz_emcal_pos->Fill(matched_emcal_phi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
          h2_z_dz_track_pos->Fill(matched_track_z.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
          h2_phi_dz_track_pos->Fill(matched_track_phi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
        }
      }

      if (fabs(matched_track_emcal_dz.at(index).at(i))<20)
      {
        h1_dphi_track_emcal->Fill(matched_track_emcal_dphi.at(index).at(i));
      }

      if (fabs(matched_track_emcal_dphi.at(index).at(i))<0.1 && fabs(matched_track_emcal_dz.at(index).at(i)<20))
      {
        h2_eOp_p->Fill(matched_track_p.at(index).at(i), matched_emcal_e.at(index).at(i) / matched_track_p.at(index).at(i));
        h1_eOp->Fill(matched_emcal_e.at(index).at(i) / matched_track_p.at(index).at(i));
        h1_e_emcal->Fill(matched_emcal_e.at(index).at(i));
      }

      if (fabs(matched_track_emcal_dz.at(index).at(i))<20)
      {
        h2_phi_tilt_dphi_track->Fill(matched_track_phi_tilt.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        h2_track_p_dphi_track->Fill(matched_track_p.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        h2_track_eta_dphi_track->Fill(matched_track_eta.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        if (matched_track_charge.at(index).at(i)>0)
        {
          h2_phi_tilt_dphi_track_pos->Fill(matched_track_phi_tilt.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        }
        else if (matched_track_charge.at(index).at(i)<0)
        {
          h2_phi_tilt_dphi_track_neg->Fill(matched_track_phi_tilt.at(index).at(i),matched_track_emcal_dphi.at(index).at(i));
        }
      }
    }

    if (matched_track_phi.at(index).size()!=0)
    {
      h1_phi_track->Fill(matched_track_phi.at(index).at(0));
    }
    if (matched_emcal_z.at(index).size()!=0)
    {
      h2_z_y_emcal->Fill(matched_emcal_z.at(index).at(0), matched_emcal_y.at(index).at(0));
    }
    if (matched_track_z.at(index).size()!=0)
    {
      h2_z_y_track->Fill(matched_track_z.at(index).at(0), matched_track_y.at(index).at(0));
    }
    if (matched_track_z_origin.at(index).size()!=0)
    {
      h1_z_origin_track->Fill(matched_track_z_origin.at(index).at(0));
    }
    if (matched_track_phi_tilt.at(index).size()!=0)
    {
      h1_track_phi_tilt->Fill(matched_track_phi_tilt.at(index).at(0));
    }
  }
  for (int index = 0; index<(vertex_z.size()); index++)
  {
    h1_z_vertex->Fill(vertex_z.at(index));
  }
  for (int index = 0; index<(mbd_z.size()); index++)
  {
    h1_z_mbd->Fill(mbd_z.at(index));
  }
  for (int index = 0; index<(matched_track_hcal_dphi.size()); index++)
  {
    for (int i = 0; i < (matched_track_hcal_dphi.at(index).size()); i++)
    {
      // must have EMcal cluster
      if (matched_track_charge.at(index).size()==0)
      {
        continue;
      }
      if (matched_track_charge.at(index).at(0)>0)
      {
        h2_dphi_dz_track_hcal_pos->Fill(matched_track_hcal_dphi.at(index).at(i),matched_track_hcal_dz.at(index).at(i));
        if (fabs(matched_track_hcal_dz.at(index).at(i))<100)
        {
          h1_dphi_track_hcal_pos->Fill(matched_track_hcal_dphi.at(index).at(i));
        }
        if (fabs(matched_track_hcal_dphi.at(index).at(i))<0.5)
        {
          h1_dz_track_hcal_pos->Fill(matched_track_hcal_dz.at(index).at(i));
        }
      }
      else if (matched_track_charge.at(index).at(0)<0)
      {
        h2_dphi_dz_track_hcal_neg->Fill(matched_track_hcal_dphi.at(index).at(i),matched_track_hcal_dz.at(index).at(i));
        if (fabs(matched_track_hcal_dz.at(index).at(i))<100)
        {
          h1_dphi_track_hcal_neg->Fill(matched_track_hcal_dphi.at(index).at(i));
        }
        if (fabs(matched_track_hcal_dphi.at(index).at(i))<0.5)
        {
          h1_dz_track_hcal_neg->Fill(matched_track_hcal_dz.at(index).at(i));
        }
      }
      if (fabs(matched_track_hcal_dz.at(index).at(i))<100)
      {
        h1_dphi_track_hcal->Fill(matched_track_hcal_dphi.at(index).at(i));
      }
      if (fabs(matched_track_hcal_dphi.at(index).at(i))<0.5)
      {
        h1_dz_track_hcal->Fill(matched_track_hcal_dz.at(index).at(i));
      }
    }
  }

  TPaveText *pt = new TPaveText(.55, .72, .85, .92, "NDC");
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt->AddText("p+p #sqrt{s}=200 GeV");
  pt->AddText(Form("Run %d",runnumber));

  TCanvas *can = new TCanvas(Form("can"), "can", 800, 800);
  //can->SetLeftMargin(0.12);
  //can->SetRightMargin(0.13);
  can->cd();
  h2_dphi_dz_track_emcal->Draw("COLZ");
  pt->Draw("same");

  TCanvas *can_pos = new TCanvas(Form("can_pos"), "can_pos", 800, 600);
  //can_pos->SetLeftMargin(0.12);
  can_pos->SetRightMargin(0.15);
  can_pos->cd();
  h2_dphi_dz_track_emcal_pos->Draw("COLZ");
  pt->Draw("same");

  TCanvas *can_neg = new TCanvas(Form("can_neg"), "can_neg", 800, 600);
  //can_neg->SetLeftMargin(0.12);
  can_neg->SetRightMargin(0.15);
  can_neg->cd();
  h2_dphi_dz_track_emcal_neg->Draw("COLZ");
  pt->Draw("same");

  TCanvas *can2 = new TCanvas(Form("can2"), "can2", 800, 800);
  //can2->SetLeftMargin(0.12);
  //can2->SetRightMargin(0.05);
  can2->cd();
  h1_dphi_track_emcal->Draw();

  TCanvas *can2_pos = new TCanvas(Form("can2_pos"), "can2_pos", 800, 800);
  //can2_pos->SetLeftMargin(0.12);
  //can2_pos->SetRightMargin(0.05);
  can2_pos->cd();
  h1_dphi_track_emcal_pos->Draw();

  TCanvas *can2_neg = new TCanvas(Form("can2_neg"), "can2_neg", 800, 800);
  //can2_neg->SetLeftMargin(0.12);
  //can2_neg->SetRightMargin(0.05);
  can2_neg->cd();
  h1_dphi_track_emcal_neg->Draw();

  TCanvas *can3 = new TCanvas(Form("can3"), "can3", 800, 800);
  //can3->SetLeftMargin(0.12);
  //can3->SetRightMargin(0.05);
  can3->cd();
  h1_dz_track_emcal->Draw();

  TCanvas *can3_pos = new TCanvas(Form("can3_pos"), "can3_pos", 800, 800);
  //can3_pos->SetLeftMargin(0.12);
  //can3_pos->SetRightMargin(0.05);
  can3_pos->cd();
  h1_dz_track_emcal_pos->Draw();

  TCanvas *can3_neg = new TCanvas(Form("can3_neg"), "can3_neg", 800, 800);
  //can3_neg->SetLeftMargin(0.12);
  //can3_neg->SetRightMargin(0.05);
  can3_neg->cd();
  h1_dz_track_emcal_neg->Draw();

  TPaveText *pt4 = new TPaveText(.25, .70, .55, .90, "NDC");
  pt4->SetFillColor(0);
  pt4->SetFillStyle(0);
  pt4->SetLineColor(0);
  pt4->SetBorderSize(0);
  pt4->SetTextColor(kWhite);
  pt4->AddText("#it{#bf{sPHENIX}} Internal");
  pt4->AddText("p+p #sqrt{s}=200 GeV");
  pt4->AddText(Form("Run %d",runnumber));

  TCanvas *can4 = new TCanvas(Form("can4"), "can4", 800, 600);
  //can4->SetLeftMargin(0.10);
  can4->SetRightMargin(0.15);
  can4->cd();
  h2_z_dz_emcal->Draw("COLZ");
  pt4->Draw("same");

  TCanvas *can4_pos = new TCanvas(Form("can4_pos"), "can4_pos", 800, 800);
  //can4_pos->SetLeftMargin(0.10);
  //can4_pos->SetRightMargin(0.16);
  can4_pos->cd();
  h2_z_dz_emcal_pos->Draw("COLZ");

  TCanvas *can4_neg = new TCanvas(Form("can4_neg"), "can4_neg", 800, 800);
  //can4_neg->SetLeftMargin(0.10);
  //can4_neg->SetRightMargin(0.16);
  can4_neg->cd();
  h2_z_dz_emcal_neg->Draw("COLZ");
  pt->Draw("same");

  TCanvas *can5 = new TCanvas(Form("can5"), "can5", 800, 800);
  //can5->SetLeftMargin(0.10);
  //can5->SetRightMargin(0.16);
  can5->cd();
  h2_phi_dphi_emcal->Draw("COLZ");

  TCanvas *can6 = new TCanvas(Form("can6"), "can6", 800, 800);
  //can6->SetLeftMargin(0.10);
  //can6->SetRightMargin(0.16);
  can6->cd();
  h2_phi_track_emcal->Draw("COLZ");

  TCanvas *can7 = new TCanvas(Form("can7"), "can7", 800, 800);
  //can7->SetLeftMargin(0.10);
  //can7->SetRightMargin(0.05);
  can7->cd();
  h1_phi_track->Draw();

  TCanvas *can8_pos = new TCanvas(Form("can8_pos"), "can8_pos", 800, 800);
  //can8_pos->SetLeftMargin(0.12);
  //can8_pos->SetRightMargin(0.13);
  can8_pos->cd();
  h2_dphi_deta_track_emcal_pos->Draw("COLZ");

  TCanvas *can8_neg = new TCanvas(Form("can8_neg"), "can8_neg", 800, 800);
  //can8_neg->SetLeftMargin(0.12);
  //can8_neg->SetRightMargin(0.13);
  can8_neg->cd();
  h2_dphi_deta_track_emcal_neg->Draw("COLZ");

  TCanvas *can9_pos = new TCanvas(Form("can9_pos"), "can9_pos", 800, 800);
  //can9_pos->SetLeftMargin(0.12);
  //can9_pos->SetRightMargin(0.05);
  can9_pos->cd();
  h1_deta_track_emcal_pos->Draw();

  TCanvas *can9_neg = new TCanvas(Form("can9_neg"), "can9_neg", 800, 800);
  //can9_neg->SetLeftMargin(0.12);
  //can9_neg->SetRightMargin(0.05);
  can9_neg->cd();
  h1_deta_track_emcal_neg->Draw();

  TCanvas *can10 = new TCanvas(Form("can10"), "can10", 800, 800);
  //can10->SetLeftMargin(0.10);
  //can10->SetRightMargin(0.16);
  can10->cd();
  h2_eOp_p->Draw("COLZ");

  TCanvas *can11 = new TCanvas(Form("can11"), "can11", 800, 800);
  //can11->SetLeftMargin(0.10);
  //can11->SetRightMargin(0.16);
  can11->cd();
  h2_z_track_emcal->Draw("COLZ");

  TCanvas *can12 = new TCanvas(Form("can12"), "can12", 800, 800);
  //can12->SetLeftMargin(0.12);
  //can12->SetRightMargin(0.05);
  can12->cd();
  h1_eOp->Draw("COLZ");

  TCanvas *can13 = new TCanvas(Form("can13"), "can13", 800, 800);
  //can13->SetLeftMargin(0.10);
  //can13->SetRightMargin(0.16);
  can13->cd();
  h2_z_dz_track->Draw("COLZ");

  TCanvas *can13_pos = new TCanvas(Form("can13_pos"), "can13_pos", 800, 800);
  //can13_pos->SetLeftMargin(0.10);
  //can13_pos->SetRightMargin(0.16);
  can13_pos->cd();
  h2_z_dz_track_pos->Draw("COLZ");

  TCanvas *can13_neg = new TCanvas(Form("can13_neg"), "can13_neg", 800, 800);
  //can13_neg->SetLeftMargin(0.10);
  //can13_neg->SetRightMargin(0.16);
  can13_neg->cd();
  h2_z_dz_track_neg->Draw("COLZ");

  TCanvas *can14 = new TCanvas(Form("can14"), "can14", 800, 800);
  //can14->SetLeftMargin(0.10);
  //can14->SetRightMargin(0.16);
  can14->cd();
  h2_phi_dz_emcal->Draw("COLZ");

  TCanvas *can14_pos = new TCanvas(Form("can14_pos"), "can14_pos", 800, 800);
  //can14_pos->SetLeftMargin(0.10);
  //can14_pos->SetRightMargin(0.16);
  can14_pos->cd();
  h2_phi_dz_emcal_pos->Draw("COLZ");

  TCanvas *can14_neg = new TCanvas(Form("can14_neg"), "can14_neg", 800, 800);
  //can14_neg->SetLeftMargin(0.10);
  //can14_neg->SetRightMargin(0.16);
  can14_neg->cd();
  h2_phi_dz_emcal_neg->Draw("COLZ");

  TCanvas *can15 = new TCanvas(Form("can15"), "can15", 800, 800);
  //can15->SetLeftMargin(0.10);
  //can15->SetRightMargin(0.16);
  can15->cd();
  h2_phi_dz_track->Draw("COLZ");

  TCanvas *can15_pos = new TCanvas(Form("can15_pos"), "can15_pos", 800, 800);
  //can15_pos->SetLeftMargin(0.10);
  //can15_pos->SetRightMargin(0.16);
  can15_pos->cd();
  h2_phi_dz_track_pos->Draw("COLZ");

  TCanvas *can15_neg = new TCanvas(Form("can15_neg"), "can15_neg", 800, 800);
  //can15_neg->SetLeftMargin(0.10);
  //can15_neg->SetRightMargin(0.16);
  can15_neg->cd();
  h2_phi_dz_track_neg->Draw("COLZ");

  TCanvas *can16 = new TCanvas(Form("can16"), "can16", 800, 800);
  //can16->SetLeftMargin(0.10);
  //can16->SetRightMargin(0.16);
  can16->cd();
  h2_z_dphi_emcal->Draw("COLZ");

  TCanvas *can16_pos = new TCanvas(Form("can16_pos"), "can16_pos", 800, 800);
  //can16_pos->SetLeftMargin(0.10);
  //can16_pos->SetRightMargin(0.16);
  can16_pos->cd();
  h2_z_dphi_emcal_pos->Draw("COLZ");

  TCanvas *can16_neg = new TCanvas(Form("can16_neg"), "can16_neg", 800, 800);
  //can16_neg->SetLeftMargin(0.10);
  //can16_neg->SetRightMargin(0.16);
  can16_neg->cd();
  h2_z_dphi_emcal_neg->Draw("COLZ");

  TCanvas *can17 = new TCanvas(Form("can17"), "can17", 800, 800);
  //can17->SetLeftMargin(0.10);
  //can17->SetRightMargin(0.16);
  can17->cd();
  h2_z_dphi_track->Draw("COLZ");

  TCanvas *can17_pos = new TCanvas(Form("can17_pos"), "can17_pos", 800, 800);
  //can17_pos->SetLeftMargin(0.10);
  //can17_pos->SetRightMargin(0.16);
  can17_pos->cd();
  h2_z_dphi_track_pos->Draw("COLZ");

  TCanvas *can17_neg = new TCanvas(Form("can17_neg"), "can17_neg", 800, 800);
  //can17_neg->SetLeftMargin(0.10);
  //can17_neg->SetRightMargin(0.16);
  can17_neg->cd();
  h2_z_dphi_track_neg->Draw("COLZ");

  TCanvas *can18 = new TCanvas(Form("can18"), "can18", 800, 800);
  //can18->SetLeftMargin(0.10);
  //can18->SetRightMargin(0.16);
  can18->cd();
  h2_phi_dphi_track->Draw("COLZ");

  TCanvas *can19 = new TCanvas(Form("can19"), "can19", 800, 800);
  //can19->SetLeftMargin(0.10);
  //can19->SetRightMargin(0.16);
  can19->cd();
  h2_z_y_emcal->Draw("COLZ");

  TCanvas *can20 = new TCanvas(Form("can20"), "can20", 800, 800);
  //can20->SetLeftMargin(0.10);
  //can20->SetRightMargin(0.16);
  can20->cd();
  h2_z_y_track->Draw("COLZ");

  TCanvas *can21 = new TCanvas(Form("can21"), "can21", 800, 800);
  //can21->SetLeftMargin(0.12);
  //can21->SetRightMargin(0.05);
  can21->cd();
  h1_z_origin_track->Draw();

  TCanvas *can22 = new TCanvas(Form("can22"), "can22", 800, 800);
  //can22->SetLeftMargin(0.12);
  //can22->SetRightMargin(0.05);
  can22->cd();
  h1_z_vertex->Draw();

  TCanvas *can22_a = new TCanvas(Form("can22_a"), "can22_a", 800, 800);
  //can22_a->SetLeftMargin(0.12);
  //can22_a->SetRightMargin(0.05);
  can22_a->cd();
  h1_z_mbd->Draw();

  TCanvas *can23 = new TCanvas(Form("can23"), "can23", 800, 800);
  //can23->SetLeftMargin(0.12);
  //can23->SetRightMargin(0.05);
  can23->cd();
  h1_e_emcal->Draw();

  TCanvas *can24_pos = new TCanvas(Form("can24_pos"), "can24_pos", 800, 800);
  //can24_pos->SetLeftMargin(0.12);
  //can24_pos->SetRightMargin(0.13);
  can24_pos->cd();
  h2_dphi_dz_track_hcal_pos->Draw("COLZ");

  TCanvas *can24_neg = new TCanvas(Form("can24_neg"), "can24_neg", 800, 800);
  //can24_neg->SetLeftMargin(0.12);
  //can24_neg->SetRightMargin(0.13);
  can24_neg->cd();
  h2_dphi_dz_track_hcal_neg->Draw("COLZ");

  TCanvas *can25 = new TCanvas(Form("can25"), "can25", 800, 800);
  //can25->SetLeftMargin(0.12);
  //can25->SetRightMargin(0.05);
  can25->cd();
  h1_dphi_track_hcal->Draw();

  TCanvas *can25_pos = new TCanvas(Form("can25_pos"), "can25_pos", 800, 800);
  //can25_pos->SetLeftMargin(0.12);
  //can25_pos->SetRightMargin(0.05);
  can25_pos->cd();
  h1_dphi_track_hcal_pos->Draw();

  TCanvas *can25_neg = new TCanvas(Form("can25_neg"), "can25_neg", 800, 800);
  //can25_neg->SetLeftMargin(0.12);
  //can25_neg->SetRightMargin(0.05);
  can25_neg->cd();
  h1_dphi_track_hcal_neg->Draw();

  TCanvas *can26 = new TCanvas(Form("can26"), "can26", 800, 800);
  //can26->SetLeftMargin(0.12);
  //can26->SetRightMargin(0.05);
  can26->cd();
  h1_dz_track_hcal->Draw();

  TCanvas *can26_pos = new TCanvas(Form("can26_pos"), "can26_pos", 800, 800);
  //can26_pos->SetLeftMargin(0.12);
  //can26_pos->SetRightMargin(0.05);
  can26_pos->cd();
  h1_dz_track_hcal_pos->Draw();

  TCanvas *can26_neg = new TCanvas(Form("can26_neg"), "can26_neg", 800, 800);
  //can26_neg->SetLeftMargin(0.12);
  //can26_neg->SetRightMargin(0.05);
  can26_neg->cd();
  h1_dz_track_hcal_neg->Draw();

  TCanvas *can27 = new TCanvas(Form("can27"), "can27", 800, 800);
  //can27->SetLeftMargin(0.12);
  //can27->SetRightMargin(0.05);
  can27->cd();
  h1_track_phi_tilt->Draw();

  TPaveText *pt28 = new TPaveText(.30, .70, .60, .90, "NDC");
  pt28->SetFillColor(0);
  pt28->SetFillStyle(0);
  pt28->SetLineColor(0);
  pt28->SetBorderSize(0);
  pt28->SetTextColor(kWhite);
  pt28->AddText("#it{#bf{sPHENIX}} Internal");
  pt28->AddText("p+p #sqrt{s}=200 GeV");
  pt28->AddText(Form("Run %d",runnumber));

  TCanvas *can28 = new TCanvas(Form("can28"), "can28", 800, 600);
  //can28->SetLeftMargin(0.10);
  can28->SetRightMargin(0.15);
  can28->cd();
  h2_phi_tilt_dphi_track->Draw("COLZ");
  pt28->Draw("same");

  TCanvas *can28_pos = new TCanvas(Form("can28_pos"), "can28_pos", 800, 800);
  //can28_pos->SetLeftMargin(0.10);
  //can28_pos->SetRightMargin(0.16);
  can28_pos->cd();
  h2_phi_tilt_dphi_track_pos->Draw("COLZ");

  TCanvas *can28_neg = new TCanvas(Form("can28_neg"), "can28_neg", 800, 800);
  //can28_neg->SetLeftMargin(0.10);
  //can28_neg->SetRightMargin(0.16);
  can28_neg->cd();
  h2_phi_tilt_dphi_track_neg->Draw("COLZ");

  TCanvas *can29 = new TCanvas(Form("can29"), "can29", 800, 800);
  //can29->SetLeftMargin(0.10);
  //can29->SetRightMargin(0.16);
  can29->cd();
  can29->SetLogz(1);
  h2_track_p_dphi_track->Draw("COLZ");

  TCanvas *can30 = new TCanvas(Form("can30"), "can30", 800, 800);
  //can30->SetLeftMargin(0.10);
  //can30->SetRightMargin(0.16);
  can30->cd();
  can30->SetLogz(1);
  h2_track_eta_dphi_track->Draw("COLZ");

  fs::path dir = Form("figure/%d",runnumber);

  if (!fs::exists(dir)) {
      if (fs::create_directory(dir)) {
          std::cout << Form("Directory 'figure/%d' created successfully.\n",runnumber);
      } else {
          std::cerr << "Failed to create directory 'aaa'.\n";
      }
  } else {
      std::cout << Form("Directory 'figure/%d' already exists.\n",runnumber);
  }

  can->Update();
  can->SaveAs(Form("figure/%d/TrackEMcal_dphi_dz_run%d.pdf",runnumber,runnumber));

  can_pos->Update();
  can_pos->SaveAs(Form("figure/%d/TrackEMcal_dphi_dz_run%d_pos.pdf",runnumber,runnumber));

  can_neg->Update();
  can_neg->SaveAs(Form("figure/%d/TrackEMcal_dphi_dz_run%d_neg.pdf",runnumber,runnumber));

  can2->Update();
  can2->SaveAs(Form("figure/%d/TrackEMcal_dphi_run%d.pdf",runnumber,runnumber));

  can2_pos->Update();
  can2_pos->SaveAs(Form("figure/%d/TrackEMcal_dphi_run%d_pos.pdf",runnumber,runnumber));

  can2_neg->Update();
  can2_neg->SaveAs(Form("figure/%d/TrackEMcal_dphi_run%d_neg.pdf",runnumber,runnumber));

  can3->Update();
  can3->SaveAs(Form("figure/%d/TrackEMcal_dz_run%d.pdf",runnumber,runnumber));

  can3_pos->Update();
  can3_pos->SaveAs(Form("figure/%d/TrackEMcal_dz_run%d_pos.pdf",runnumber,runnumber));

  can3_neg->Update();
  can3_neg->SaveAs(Form("figure/%d/TrackEMcal_dz_run%d_neg.pdf",runnumber,runnumber));

  can4->Update();
  can4->SaveAs(Form("figure/%d/TrackEMcal_z_dz_emcal_run%d.pdf",runnumber,runnumber));

  can4_pos->Update();
  can4_pos->SaveAs(Form("figure/%d/TrackEMcal_z_dz_emcal_run%d_pos.pdf",runnumber,runnumber));

  can4_neg->Update();
  can4_neg->SaveAs(Form("figure/%d/TrackEMcal_z_dz_emcal_run%d_neg.pdf",runnumber,runnumber));

  can5->Update();
  can5->SaveAs(Form("figure/%d/TrackEMcal_phi_dphi_emcal_run%d.pdf",runnumber,runnumber));

  can6->Update();
  can6->SaveAs(Form("figure/%d/TrackEMcal_phi_trackemcal_run%d.pdf",runnumber,runnumber));

  can7->Update();
  can7->SaveAs(Form("figure/%d/TrackEMcal_phi_track_run%d.pdf",runnumber,runnumber));

  can8_pos->Update();
  can8_pos->SaveAs(Form("figure/%d/TrackEMcal_dphi_deta_run%d_pos.pdf",runnumber,runnumber));

  can8_neg->Update();
  can8_neg->SaveAs(Form("figure/%d/TrackEMcal_dphi_deta_run%d_neg.pdf",runnumber,runnumber));

  can9_pos->Update();
  can9_pos->SaveAs(Form("figure/%d/TrackEMcal_deta_run%d_pos.pdf",runnumber,runnumber));

  can9_neg->Update();
  can9_neg->SaveAs(Form("figure/%d/TrackEMcal_deta_run%d_neg.pdf",runnumber,runnumber));

  can10->Update();
  can10->SaveAs(Form("figure/%d/TrackEMcal_eOp_p_run%d.pdf",runnumber,runnumber));

  can11->Update();
  can11->SaveAs(Form("figure/%d/TrackEMcal_z_trackemcal_run%d.pdf",runnumber,runnumber));

  can12->Update();
  can12->SaveAs(Form("figure/%d/TrackEMcal_eOp_run%d.pdf",runnumber,runnumber));

  can13->Update();
  can13->SaveAs(Form("figure/%d/TrackEMcal_z_dz_track_run%d.pdf",runnumber,runnumber));

  can13_pos->Update();
  can13_pos->SaveAs(Form("figure/%d/TrackEMcal_z_dz_track_run%d_pos.pdf",runnumber,runnumber));

  can13_neg->Update();
  can13_neg->SaveAs(Form("figure/%d/TrackEMcal_z_dz_track_run%d_neg.pdf",runnumber,runnumber));

  can14->Update();
  can14->SaveAs(Form("figure/%d/TrackEMcal_phi_dz_emcal_run%d.pdf",runnumber,runnumber));

  can14_pos->Update();
  can14_pos->SaveAs(Form("figure/%d/TrackEMcal_phi_dz_emcal_run%d_pos.pdf",runnumber,runnumber));

  can14_neg->Update();
  can14_neg->SaveAs(Form("figure/%d/TrackEMcal_phi_dz_emcal_run%d_neg.pdf",runnumber,runnumber));

  can15->Update();
  can15->SaveAs(Form("figure/%d/TrackEMcal_phi_dz_track_run%d.pdf",runnumber,runnumber));

  can15_pos->Update();
  can15_pos->SaveAs(Form("figure/%d/TrackEMcal_phi_dz_track_run%d_pos.pdf",runnumber,runnumber));

  can15_neg->Update();
  can15_neg->SaveAs(Form("figure/%d/TrackEMcal_phi_dz_track_run%d_neg.pdf",runnumber,runnumber));

  can16->Update();
  can16->SaveAs(Form("figure/%d/TrackEMcal_z_dphi_emcal_run%d.pdf",runnumber,runnumber));

  can16_pos->Update();
  can16_pos->SaveAs(Form("figure/%d/TrackEMcal_z_dphi_emcal_run%d_pos.pdf",runnumber,runnumber));

  can16_neg->Update();
  can16_neg->SaveAs(Form("figure/%d/TrackEMcal_z_dphi_emcal_run%d_neg.pdf",runnumber,runnumber));

  can17->Update();
  can17->SaveAs(Form("figure/%d/TrackEMcal_z_dphi_track_run%d.pdf",runnumber,runnumber));

  can17_pos->Update();
  can17_pos->SaveAs(Form("figure/%d/TrackEMcal_z_dphi_track_run%d_pos.pdf",runnumber,runnumber));

  can17_neg->Update();
  can17_neg->SaveAs(Form("figure/%d/TrackEMcal_z_dphi_track_run%d_neg.pdf",runnumber,runnumber));

  can18->Update();
  can18->SaveAs(Form("figure/%d/TrackEMcal_phi_dphi_track_run%d.pdf",runnumber,runnumber));

  can19->Update();
  can19->SaveAs(Form("figure/%d/TrackEMcal_z_y_emcal_run%d.pdf",runnumber,runnumber));

  can20->Update();
  can20->SaveAs(Form("figure/%d/TrackEMcal_z_y_track_run%d.pdf",runnumber,runnumber));

  can21->Update();
  can21->SaveAs(Form("figure/%d/TrackEMcal_z_origin_track_run%d.pdf",runnumber,runnumber));

  can22->Update();
  can22->SaveAs(Form("figure/%d/TrackEMcal_z_vertex_run%d.pdf",runnumber,runnumber));

  can22_a->Update();
  can22_a->SaveAs(Form("figure/%d/TrackEMcal_z_mbd_run%d.pdf",runnumber,runnumber));

  can23->Update();
  can23->SaveAs(Form("figure/%d/TrackEMcal_e_emcal_run%d.pdf",runnumber,runnumber));

  can24_pos->Update();
  can24_pos->SaveAs(Form("figure/%d/TrackHcal_dphi_dz_run%d_pos.pdf",runnumber,runnumber));

  can24_neg->Update();
  can24_neg->SaveAs(Form("figure/%d/TrackHcal_dphi_dz_run%d_neg.pdf",runnumber,runnumber));

  can25->Update();
  can25->SaveAs(Form("figure/%d/TrackHcal_dphi_run%d.pdf",runnumber,runnumber));

  can25_pos->Update();
  can25_pos->SaveAs(Form("figure/%d/TrackHcal_dphi_run%d_pos.pdf",runnumber,runnumber));

  can25_neg->Update();
  can25_neg->SaveAs(Form("figure/%d/TrackHcal_dphi_run%d_neg.pdf",runnumber,runnumber));

  can26->Update();
  can26->SaveAs(Form("figure/%d/TrackHcal_dz_run%d.pdf",runnumber,runnumber));

  can26_pos->Update();
  can26_pos->SaveAs(Form("figure/%d/TrackHcal_dz_run%d_pos.pdf",runnumber,runnumber));

  can26_neg->Update();
  can26_neg->SaveAs(Form("figure/%d/TrackHcal_dz_run%d_neg.pdf",runnumber,runnumber));

  can27->Update();
  can27->SaveAs(Form("figure/%d/TrackEMcal_phi_tilt_track_run%d.pdf",runnumber,runnumber));

  can28->Update();
  can28->SaveAs(Form("figure/%d/TrackEMcal_phi_tilt_dphi_track_run%d.png",runnumber,runnumber));

  can28_pos->Update();
  can28_pos->SaveAs(Form("figure/%d/TrackEMcal_phi_tilt_dphi_track_run%d_pos.pdf",runnumber,runnumber));

  can28_neg->Update();
  can28_neg->SaveAs(Form("figure/%d/TrackEMcal_phi_tilt_dphi_track_run%d_neg.pdf",runnumber,runnumber));

  can29->Update();
  can29->SaveAs(Form("figure/%d/TrackEMcal_track_p_dphi_run%d.pdf",runnumber,runnumber));

  can30->Update();
  can30->SaveAs(Form("figure/%d/TrackEMcal_track_eta_dphi_run%d.pdf",runnumber,runnumber));

  TCanvas *can233 = new TCanvas(Form("can233"), "can233", 800, 800);
  //can233->SetLeftMargin(0.12);
  //can233->SetRightMargin(0.05);
  can233->cd();
  h1_z_vertex->SetLineColor(kBlue);
  h1_z_vertex->GetXaxis()->SetTitle("Vertex Z [cm]");
  h1_z_vertex->Draw();
//  h1_z_origin_track->SetLineColor(kRed);
//  h1_z_origin_track->Scale(h1_z_vertex->Integral() / h1_z_origin_track->Integral());
//  h1_z_origin_track->Draw("same");
  h1_z_mbd->SetMarkerSize(0);
  h1_z_mbd->SetLineColor(kRed);
  //h1_z_mbd->SetFillColor(kRed);
  h1_z_mbd->Scale(h1_z_vertex->Integral() / h1_z_mbd->Integral());
  h1_z_mbd->Draw("hist,same");
  TLegend *legend = new TLegend(0.62, 0.75, 0.92, 0.9);
  legend->AddEntry(h1_z_vertex, "Svtx Vertex", "l");
//  legend->AddEntry(h1_z_origin_track, "Track project to R=0", "lep");
  legend->AddEntry(h1_z_mbd, "MBD Vertex", "l");
  legend->Draw();

  can233->Update();
  can233->SaveAs(Form("figure/%d/TrackEMcal_z_vertex_svtxVSmbd_run%d.pdf",runnumber,runnumber));

}
