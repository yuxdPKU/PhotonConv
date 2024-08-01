#include <filesystem>
#include "utilities.h"

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void saveroot(int runnumber)
{
  gStyle->SetOptStat(0);

  TString inputfile = Form("%d/final_%d_ana.root",runnumber,runnumber);

  TFile *file = new TFile(inputfile, "READ");
  TTree *tree = (TTree*)file->Get("tree");

  std::vector<float> *_vertex_x = 0;
  std::vector<float> *_vertex_y = 0;
  std::vector<float> *_vertex_z = 0;
  std::vector<int> *_track_id = 0;
  std::vector<float> *_track_quality = 0;
  std::vector<float> *_track_dcaxy = 0;
  std::vector<float> *_track_dcaz = 0;
  std::vector<int> *_track_nc_tpc = 0;
  std::vector<int> *_track_nc_mvtx = 0;
  std::vector<int> *_track_nc_intt = 0;
  std::vector<int> *_track_bc = 0;
  std::vector<float> *_track_phi = 0;
  std::vector<float> *_track_eta = 0;
  std::vector<float> *_track_x = 0;
  std::vector<float> *_track_y = 0;
  std::vector<float> *_track_z = 0;
  std::vector<float> *_track_ptq = 0;
  std::vector<float> *_track_px = 0;
  std::vector<float> *_track_py = 0;
  std::vector<float> *_track_pz = 0;
  std::vector<float> *_track_phi_origin = 0;
  std::vector<float> *_track_eta_origin = 0;
  std::vector<float> *_track_x_origin = 0;
  std::vector<float> *_track_y_origin = 0;
  std::vector<float> *_track_z_origin = 0;
  std::vector<float> *_track_phi_emc = 0;
  std::vector<float> *_track_eta_emc = 0;
  std::vector<float> *_track_x_emc = 0;
  std::vector<float> *_track_y_emc = 0;
  std::vector<float> *_track_z_emc = 0;
  std::vector<float> *_track_phi_ihc = 0;
  std::vector<float> *_track_eta_ihc = 0;
  std::vector<float> *_track_x_ihc = 0;
  std::vector<float> *_track_y_ihc = 0;
  std::vector<float> *_track_z_ihc = 0;
  std::vector<float> *_track_phi_ohc = 0;
  std::vector<float> *_track_eta_ohc = 0;
  std::vector<float> *_track_x_ohc = 0;
  std::vector<float> *_track_y_ohc = 0;
  std::vector<float> *_track_z_ohc = 0;
  std::vector<int> *_trClus_track_id = 0;
  std::vector<int> *_trClus_type = 0;
  std::vector<float> *_trClus_x = 0;
  std::vector<float> *_trClus_y = 0;
  std::vector<float> *_trClus_z = 0;
  std::vector<int> *_emcal_id = 0;
  std::vector<float> *_emcal_phi = 0;
  std::vector<float> *_emcal_eta = 0;
  std::vector<float> *_emcal_x = 0;
  std::vector<float> *_emcal_y = 0;
  std::vector<float> *_emcal_z = 0;
  std::vector<float> *_emcal_e = 0;
  std::vector<int> *_emcal_tower_cluster_id = 0;
  std::vector<float> *_emcal_tower_e = 0;
  std::vector<float> *_emcal_tower_phi = 0;
  std::vector<float> *_emcal_tower_eta = 0;
  std::vector<int> *_emcal_tower_status = 0;
  std::vector<int> *_hcal_id = 0;
  std::vector<float> *_hcal_phi = 0;
  std::vector<float> *_hcal_eta = 0;
  std::vector<float> *_hcal_x = 0;
  std::vector<float> *_hcal_y = 0;
  std::vector<float> *_hcal_z = 0;
  std::vector<float> *_hcal_e = 0;
  std::vector<int> *_hcal_tower_cluster_id = 0;
  std::vector<float> *_hcal_tower_e = 0;
  std::vector<float> *_hcal_tower_phi = 0;
  std::vector<float> *_hcal_tower_eta = 0;
  std::vector<int> *_hcal_tower_status = 0;
  std::vector<int> *_ntracks = 0;
  std::vector<float> *_mbd_z = 0;
  std::vector<int> *_triggers = 0;

  tree->SetBranchAddress("_vertex_x", &_vertex_x);
  tree->SetBranchAddress("_vertex_y", &_vertex_y);
  tree->SetBranchAddress("_vertex_z", &_vertex_z);
  tree->SetBranchAddress("_track_id", &_track_id);
  tree->SetBranchAddress("_track_quality", &_track_quality);
  tree->SetBranchAddress("_track_dcaxy", &_track_dcaxy);
  tree->SetBranchAddress("_track_dcaz", &_track_dcaz);
  tree->SetBranchAddress("_track_nc_tpc", &_track_nc_tpc);
  tree->SetBranchAddress("_track_nc_mvtx", &_track_nc_mvtx);
  tree->SetBranchAddress("_track_nc_intt", &_track_nc_intt);
  tree->SetBranchAddress("_track_bc", &_track_bc);
  tree->SetBranchAddress("_track_phi", &_track_phi);
  tree->SetBranchAddress("_track_eta", &_track_eta);
  tree->SetBranchAddress("_track_x", &_track_x);
  tree->SetBranchAddress("_track_y", &_track_y);
  tree->SetBranchAddress("_track_z", &_track_z);
  tree->SetBranchAddress("_track_ptq", &_track_ptq);
  tree->SetBranchAddress("_track_px", &_track_px);
  tree->SetBranchAddress("_track_py", &_track_py);
  tree->SetBranchAddress("_track_pz", &_track_pz);
  tree->SetBranchAddress("_track_phi_origin", &_track_phi_origin);
  tree->SetBranchAddress("_track_eta_origin", &_track_eta_origin);
  tree->SetBranchAddress("_track_x_origin", &_track_x_origin);
  tree->SetBranchAddress("_track_y_origin", &_track_y_origin);
  tree->SetBranchAddress("_track_z_origin", &_track_z_origin);
  tree->SetBranchAddress("_track_phi_emc", &_track_phi_emc);
  tree->SetBranchAddress("_track_eta_emc", &_track_eta_emc);
  tree->SetBranchAddress("_track_x_emc", &_track_x_emc);
  tree->SetBranchAddress("_track_y_emc", &_track_y_emc);
  tree->SetBranchAddress("_track_z_emc", &_track_z_emc);
  tree->SetBranchAddress("_track_phi_ihc", &_track_phi_ihc);
  tree->SetBranchAddress("_track_eta_ihc", &_track_eta_ihc);
  tree->SetBranchAddress("_track_x_ihc", &_track_x_ihc);
  tree->SetBranchAddress("_track_y_ihc", &_track_y_ihc);
  tree->SetBranchAddress("_track_z_ihc", &_track_z_ihc);
  tree->SetBranchAddress("_track_phi_ohc", &_track_phi_ohc);
  tree->SetBranchAddress("_track_eta_ohc", &_track_eta_ohc);
  tree->SetBranchAddress("_track_x_ohc", &_track_x_ohc);
  tree->SetBranchAddress("_track_y_ohc", &_track_y_ohc);
  tree->SetBranchAddress("_track_z_ohc", &_track_z_ohc);
  tree->SetBranchAddress("_trClus_track_id", &_trClus_track_id);
  tree->SetBranchAddress("_trClus_type", &_trClus_type);
  tree->SetBranchAddress("_trClus_x", &_trClus_x);
  tree->SetBranchAddress("_trClus_y", &_trClus_y);
  tree->SetBranchAddress("_trClus_z", &_trClus_z);
  tree->SetBranchAddress("_emcal_id", &_emcal_id);
  tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
  tree->SetBranchAddress("_emcal_eta", &_emcal_eta);
  tree->SetBranchAddress("_emcal_x", &_emcal_x);
  tree->SetBranchAddress("_emcal_y", &_emcal_y);
  tree->SetBranchAddress("_emcal_z", &_emcal_z);
  tree->SetBranchAddress("_emcal_e", &_emcal_e);
  tree->SetBranchAddress("_emcal_tower_cluster_id", &_emcal_tower_cluster_id);
  tree->SetBranchAddress("_emcal_tower_e", &_emcal_tower_e);
  tree->SetBranchAddress("_emcal_tower_phi", &_emcal_tower_phi);
  tree->SetBranchAddress("_emcal_tower_eta", &_emcal_tower_eta);
  tree->SetBranchAddress("_emcal_tower_status", &_emcal_tower_status);
  tree->SetBranchAddress("_hcal_id", &_hcal_id);
  tree->SetBranchAddress("_hcal_phi", &_hcal_phi);
  tree->SetBranchAddress("_hcal_eta", &_hcal_eta);
  tree->SetBranchAddress("_hcal_x", &_hcal_x);
  tree->SetBranchAddress("_hcal_y", &_hcal_y);
  tree->SetBranchAddress("_hcal_z", &_hcal_z);
  tree->SetBranchAddress("_hcal_e", &_hcal_e);
  tree->SetBranchAddress("_hcal_tower_cluster_id", &_hcal_tower_cluster_id);
  tree->SetBranchAddress("_hcal_tower_e", &_hcal_tower_e);
  tree->SetBranchAddress("_hcal_tower_phi", &_hcal_tower_phi);
  tree->SetBranchAddress("_hcal_tower_eta", &_hcal_tower_eta);
  tree->SetBranchAddress("_hcal_tower_status", &_hcal_tower_status);
  tree->SetBranchAddress("_mbd_z", &_mbd_z);
  tree->SetBranchAddress("_triggers", &_triggers);
  tree->SetBranchAddress("_ntracks", &_ntracks);

  TFile* outputfile = new TFile("var.root","recreate");
  TTree* outtree = new TTree("tree","");

  float m_dphi, m_dz;
  int m_charge;
  float m_phi_tilt;

  outtree->Branch("dphi",&m_dphi,"dphi/F");
  outtree->Branch("dz",&m_dz,"dz/F");
  outtree->Branch("charge",&m_charge,"charge/I");
  outtree->Branch("phi_tilt",&m_phi_tilt,"phi_tilt/F");

  for(int i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);

    if(_ntracks->at(0) > 800)
    {
      std::cout << "Event: " << i << " rejected by number of tracks: " << _ntracks->at(0) << std::endl;
      continue;
    }

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
      TVector3 p3_track(_track_px->at(itrack),_track_py->at(itrack),_track_pz->at(itrack));
      TVector3 R3_track_emc(_track_x_emc->at(itrack),_track_y_emc->at(itrack),0);

      float track_p = sqrt(pow(_track_px->at(itrack),2) + pow(_track_py->at(itrack),2) + pow(_track_pz->at(itrack),2));

      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
        // cut good emcal clusters
        //if(_emcal_e->at(iem) < min_EMCal_E) continue;
        //if(fabs(_emcal_eta->at(iem)) > 1.1) continue;
        std::pair<float, float> EMCalPos;
        float emcal_phi = atan2(_emcal_y->at(iem), _emcal_x->at(iem));
        EMCalPos = std::make_pair(emcal_phi, _emcal_z->at(iem));

        float dphi = PiRange(TrackProjsEMCal.first - EMCalPos.first);
        float dz = TrackProjsEMCal.second - EMCalPos.second;
        float deta = _track_eta_emc->at(itrack) - _emcal_eta->at(iem);

        m_dphi = dphi;
        m_dz = dz;
        m_charge = _track_ptq->at(itrack) > 0? 1 : -1;
        m_phi_tilt = p3_track.Dot(R3_track_emc.Cross(z_direction)) / p3_track.Dot(R3_track_emc);
        outtree->Fill();
      }
    }
  }

  outputfile->Write();
  outputfile->Close();

}
