#include <filesystem>
#include "utilities.h"

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void saveroot(int runnumber)
{
  gStyle->SetOptStat(0);

  //TString inputfile = Form("/sphenix/u/xyu3/hftg01/PhotonConv/macro/DVC_0_14/Reconstructed/%d/final_%d_ana.root",runnumber,runnumber);
  //TString inputfile = Form("%d_woRcorr/final_%d_ana.root",runnumber,runnumber);
  //TString inputfile = Form("%d_0p077m/final_%d_ana.root",runnumber,runnumber);
  TString inputfile = Form("%d/final_%d_ana.root",runnumber,runnumber);
  //TString inputfile = Form("%d_old/final_%d_ana.root",runnumber,runnumber);
  //TString inputfile = Form("%d/TrackCalo_0_ana.root",runnumber);

  TFile *file = new TFile(inputfile, "READ");
  TTree *tree = (TTree*)file->Get("tree");

  setBranch(tree);

  TFile* outputfile = new TFile(Form("var_%d.root",runnumber),"recreate");
  TTree* outtree = new TTree("tree","");

  float m_dphi, m_dz;
  float m_calo_z;
  float m_track_z;
  int m_charge;
  float m_phi_tilt;
  float m_p, m_eta, m_e;

  outtree->Branch("dphi",&m_dphi,"dphi/F");
  outtree->Branch("dz",&m_dz,"dz/F");
  outtree->Branch("calo_z",&m_calo_z,"calo_z/F");
  outtree->Branch("track_z",&m_track_z,"track_z/F");
  outtree->Branch("charge",&m_charge,"charge/I");
  outtree->Branch("phi_tilt",&m_phi_tilt,"phi_tilt/F");
  outtree->Branch("p",&m_p,"p/F");
  outtree->Branch("e",&m_e,"e/F");
  outtree->Branch("eta",&m_eta,"eta/F");

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
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        EMCalPos = std::make_pair(emcal_phi, radius_scale*_emcal_z->at(iem));

        float dphi = PiRange(TrackProjsEMCal.first - EMCalPos.first);
        float dz = TrackProjsEMCal.second - EMCalPos.second;
        float deta = _track_eta_emc->at(itrack) - _emcal_eta->at(iem);

        m_dphi = dphi;
        m_dz = dz;
        m_charge = _track_ptq->at(itrack) > 0? 1 : -1;
        m_phi_tilt = p3_track.Dot(R3_track_emc.Cross(z_direction)) / p3_track.Dot(R3_track_emc);
        m_p = sqrt( pow(_track_px->at(itrack),2) + pow(_track_py->at(itrack),2) + pow(_track_pz->at(itrack),2) );
        m_e = _emcal_e->at(iem);
        m_calo_z = EMCalPos.second;
        m_track_z = TrackProjsEMCal.second;
        m_eta = _track_eta->at(itrack);
        outtree->Fill();
      }
    }
  }

  outputfile->Write();
  outputfile->Close();

}
