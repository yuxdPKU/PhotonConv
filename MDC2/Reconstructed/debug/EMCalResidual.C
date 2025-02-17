#include <filesystem>
#include "../utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void EMCalResidual(int runnumber=15)
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree");
  chain->Add("../15/TrackCalo_160_ana.root");

  setBranch(chain);

  std::vector<int> matched_eventid;
  std::vector<std::vector<float>> matched_track_emcal_dphi;
  std::vector<std::vector<float>> matched_track_emcal_dz;
  std::vector<std::vector<float>> matched_track_emcal_deta;
  std::vector<std::vector<int>> matched_track_charge;
  std::vector<std::vector<float>> matched_emcal_x;
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
  matched_emcal_x.clear();
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

  //int ievent = 722;
  //for(int i = ievent; i < (ievent+1); i++)
  for(int i = 0; i < chain->GetEntries(); i++)
  {
    chain->GetEntry(i);
if (_eventNumber!=974) continue;
cout<<"Begin entry "<<i<<endl;

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
    std::vector<float> vec_emcal_x;
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
    vec_emcal_x.clear();
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
      //if(_track_quality->at(itrack) > 100) continue;

      // we only need emcal matching, hcal is not necessary
      if(isnan(_track_phi_emc->at(itrack)))
      {
        continue;
      }

      // project tpc tracks to emcal and hcal
      std::pair<float, float> TrackProjsEMCal;
      std::pair<float, float> TrackProjsHCal;
      TrackProjsEMCal = std::make_pair( cal_phi(_track_x_emc->at(itrack), _track_y_emc->at(itrack)), _track_z_emc->at(itrack));
      TrackProjsHCal = std::make_pair( cal_phi(_track_x_ohc->at(itrack), _track_y_ohc->at(itrack)), _track_z_ohc->at(itrack));
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
        float emcal_phi = cal_phi(_emcal_x->at(iem), _emcal_y->at(iem));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        EMCalPos = std::make_pair(emcal_phi, radius_scale * _emcal_z->at(iem));

        float dphi = TrackProjsEMCal.first - EMCalPos.first;
        float dz = TrackProjsEMCal.second - EMCalPos.second;
        float deta = _track_eta_emc->at(itrack) - _emcal_eta->at(iem);

        vec_track_emcal_residual_phi.push_back(PiRange(dphi));
        vec_track_emcal_residual_z.push_back(dz);
        vec_track_emcal_residual_eta.push_back(deta);
        vec_emcal_phi.push_back(EMCalPos.first);
        vec_emcal_x.push_back(_emcal_x->at(iem) * radius_scale);
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
    matched_emcal_x.push_back(vec_emcal_x);
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
  }

for (int i = 0; i < (matched_track_emcal_dphi.size()); i++)
{
  cout<<"event "<<matched_eventid.at(i)<<endl;
for (int j = 0; j < (matched_track_emcal_dphi.at(i).size()); j++)
{
  cout<<"matched pair "<<j<<", ";
  cout<<"dphi = "<<matched_track_emcal_dphi.at(i).at(j)<<", ";
  cout<<"dz = "<<matched_track_emcal_dz.at(i).at(j)<<", ";
  cout<<"deta = "<<matched_track_emcal_deta.at(i).at(j)<<", ";
  cout<<"emcal_x = "<<matched_emcal_x.at(i).at(j)<<", ";
  cout<<"emcal_y = "<<matched_emcal_y.at(i).at(j)<<", ";
  cout<<"emcal_z = "<<matched_emcal_z.at(i).at(j)<<", ";
  cout<<"emcal_phi = "<<matched_emcal_phi.at(i).at(j)<<", ";
  cout<<"emcal_e = "<<matched_emcal_e.at(i).at(j)<<", ";
  cout<<"track_phi = "<<matched_track_phi.at(i).at(j)<<", ";
  cout<<"track_p = "<<matched_track_p.at(i).at(j)<<", ";
  cout<<"track_eta = "<<matched_track_eta.at(i).at(j)<<endl;
}
}

}
