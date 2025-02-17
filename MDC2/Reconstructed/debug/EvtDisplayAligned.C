#include <filesystem>
#include "../utilities.h"

namespace fs = std::filesystem;

void EvtDisplayAligned(int runnumber=15)
{
  gStyle->SetOptStat(0);

  TString inputfile = "../15/TrackCalo_13_ana.root";

  TFile *file = new TFile(inputfile, "READ");
  TTree *tree = (TTree*)file->Get("tree");

  setBranch(tree);

  std::vector<std::vector<float>> all_tpc_cluster_x;
  std::vector<std::vector<float>> all_tpc_cluster_y;
  std::vector<std::vector<float>> all_tpc_cluster_z;
  std::vector<std::vector<float>> all_tpc_cluster_R;
  std::vector<std::vector<float>> unmatched_tpc_cluster_x;
  std::vector<std::vector<float>> unmatched_tpc_cluster_y;
  std::vector<std::vector<float>> unmatched_tpc_cluster_z;
  std::vector<std::vector<float>> unmatched_tpc_cluster_R;
  std::vector<std::vector<float>> unmatched_track_project_emcal_x;
  std::vector<std::vector<float>> unmatched_track_project_emcal_y;
  std::vector<std::vector<float>> unmatched_track_project_emcal_z;
  std::vector<std::vector<float>> unmatched_track_project_emcal_R;
  std::vector<std::vector<float>> unmatched_emcal_cluster_x;
  std::vector<std::vector<float>> unmatched_emcal_cluster_y;
  std::vector<std::vector<float>> unmatched_emcal_cluster_z;
  std::vector<std::vector<float>> unmatched_emcal_cluster_R;
  std::vector<std::vector<float>> matched_tpc_cluster_x;
  std::vector<std::vector<float>> matched_tpc_cluster_y;
  std::vector<std::vector<float>> matched_tpc_cluster_z;
  std::vector<std::vector<float>> matched_tpc_cluster_R;
  std::vector<std::vector<float>> matched_track_project_emcal_x;
  std::vector<std::vector<float>> matched_track_project_emcal_y;
  std::vector<std::vector<float>> matched_track_project_emcal_z;
  std::vector<std::vector<float>> matched_track_project_emcal_R;
  std::vector<std::vector<float>> matched_emcal_cluster_x;
  std::vector<std::vector<float>> matched_emcal_cluster_y;
  std::vector<std::vector<float>> matched_emcal_cluster_z;
  std::vector<std::vector<float>> matched_emcal_cluster_R;
  std::vector<int> eventid;
  std::vector<float> diff_dptdpz;
  std::vector<float> avg_dptdpz;
  std::vector<float> diff_phi;
  std::vector<float> matched_track_emcal_pair_E;
  std::vector<float> matched_track_emcal_pair_p;

  matched_tpc_cluster_x.clear();
  matched_tpc_cluster_y.clear();
  matched_tpc_cluster_z.clear();
  matched_tpc_cluster_R.clear();
  matched_track_project_emcal_x.clear();
  matched_track_project_emcal_y.clear();
  matched_track_project_emcal_z.clear();
  matched_track_project_emcal_R.clear();
  matched_emcal_cluster_x.clear();
  matched_emcal_cluster_y.clear();
  matched_emcal_cluster_z.clear();
  matched_emcal_cluster_R.clear();
  all_tpc_cluster_x.clear();
  all_tpc_cluster_y.clear();
  all_tpc_cluster_z.clear();
  all_tpc_cluster_R.clear();
  unmatched_tpc_cluster_x.clear();
  unmatched_tpc_cluster_y.clear();
  unmatched_tpc_cluster_z.clear();
  unmatched_tpc_cluster_R.clear();
  unmatched_track_project_emcal_x.clear();
  unmatched_track_project_emcal_y.clear();
  unmatched_track_project_emcal_z.clear();
  unmatched_track_project_emcal_R.clear();
  unmatched_emcal_cluster_x.clear();
  unmatched_emcal_cluster_y.clear();
  unmatched_emcal_cluster_z.clear();
  unmatched_emcal_cluster_R.clear();
  eventid.clear();
  diff_dptdpz.clear();
  avg_dptdpz.clear();
  diff_phi.clear();
  matched_track_emcal_pair_E.clear();
  matched_track_emcal_pair_p.clear();

  float minpt = 100;
  float maxpt = 0;

  int ncharge_pp = 0;
  int ncharge_pn = 0;
  int ncharge_nn = 0;

  int ievent = 427;
  for(int i = ievent; i < (ievent+1); i++)
  //for(int i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
cout<<"_eventNumber = "<<_eventNumber<<endl;

    if(_ntracks->at(0) > 800)
    {
      std::cout << "Event: " << i << " rejected by number of tracks: " << _ntracks->at(0) << std::endl;
      continue;
    }

    std::vector<float> vec_track_cluster_x;
    std::vector<float> vec_track_cluster_y;
    std::vector<float> vec_track_cluster_z;
    std::vector<float> vec_track_cluster_R;
    std::vector<float> vec_track_pt;
    std::vector<float> vec_track_pz;
    std::vector<float> vec_track_phi;
    std::vector<float> vec_track_z;
    std::vector<int> vec_track_charge;
    std::vector<float> vec_matched_track_cluster_x;
    std::vector<float> vec_matched_track_cluster_y;
    std::vector<float> vec_matched_track_cluster_z;
    std::vector<float> vec_matched_track_cluster_R;
    vec_track_cluster_x.clear();
    vec_track_cluster_y.clear();
    vec_track_cluster_z.clear();
    vec_track_cluster_R.clear();
    vec_track_pt.clear();
    vec_track_pz.clear();
    vec_track_phi.clear();
    vec_track_z.clear();
    vec_track_charge.clear();
    vec_matched_track_cluster_x.clear();
    vec_matched_track_cluster_y.clear();
    vec_matched_track_cluster_z.clear();
    vec_matched_track_cluster_R.clear();

    std::vector<float> vec_cluster_x;
    std::vector<float> vec_cluster_y;
    std::vector<float> vec_cluster_z;
    std::vector<float> vec_cluster_R;
    vec_cluster_x.clear();
    vec_cluster_y.clear();
    vec_cluster_z.clear();
    vec_cluster_R.clear();

    std::vector<float> vec_emcal_x;
    std::vector<float> vec_emcal_y;
    std::vector<float> vec_emcal_z;
    std::vector<float> vec_emcal_R;
    std::vector<float> vec_matched_emcal_x;
    std::vector<float> vec_matched_emcal_y;
    std::vector<float> vec_matched_emcal_z;
    std::vector<float> vec_matched_emcal_R;
    vec_emcal_x.clear();
    vec_emcal_y.clear();
    vec_emcal_z.clear();
    vec_emcal_R.clear();
    vec_matched_emcal_x.clear();
    vec_matched_emcal_y.clear();
    vec_matched_emcal_z.clear();
    vec_matched_emcal_R.clear();

    std::vector<float> vec_hcal_x;
    std::vector<float> vec_hcal_y;
    std::vector<float> vec_hcal_z;
    std::vector<float> vec_hcal_R;
    std::vector<float> vec_matched_hcal_x;
    std::vector<float> vec_matched_hcal_y;
    std::vector<float> vec_matched_hcal_z;
    std::vector<float> vec_matched_hcal_R;
    vec_hcal_x.clear();
    vec_hcal_y.clear();
    vec_hcal_z.clear();
    vec_hcal_R.clear();
    vec_matched_hcal_x.clear();
    vec_matched_hcal_y.clear();
    vec_matched_hcal_z.clear();
    vec_matched_hcal_R.clear();

    std::vector<float> vec_track_project_emcal_x;
    std::vector<float> vec_track_project_emcal_y;
    std::vector<float> vec_track_project_emcal_z;
    std::vector<float> vec_track_project_emcal_R;
    std::vector<float> vec_matched_track_project_emcal_x;
    std::vector<float> vec_matched_track_project_emcal_y;
    std::vector<float> vec_matched_track_project_emcal_z;
    std::vector<float> vec_matched_track_project_emcal_R;
    vec_track_project_emcal_x.clear();
    vec_track_project_emcal_y.clear();
    vec_track_project_emcal_z.clear();
    vec_track_project_emcal_R.clear();
    vec_matched_track_project_emcal_x.clear();
    vec_matched_track_project_emcal_y.clear();
    vec_matched_track_project_emcal_z.clear();
    vec_matched_track_project_emcal_R.clear();

    std::vector<float> vec_track_project_hcal_x;
    std::vector<float> vec_track_project_hcal_y;
    std::vector<float> vec_track_project_hcal_z;
    std::vector<float> vec_track_project_hcal_R;
    std::vector<float> vec_matched_track_project_hcal_x;
    std::vector<float> vec_matched_track_project_hcal_y;
    std::vector<float> vec_matched_track_project_hcal_z;
    std::vector<float> vec_matched_track_project_hcal_R;
    vec_track_project_hcal_x.clear();
    vec_track_project_hcal_y.clear();
    vec_track_project_hcal_z.clear();
    vec_track_project_hcal_R.clear();
    vec_matched_track_project_hcal_x.clear();
    vec_matched_track_project_hcal_y.clear();
    vec_matched_track_project_hcal_z.clear();
    vec_matched_track_project_hcal_R.clear();

    // track + emcal matching
    if(_track_ptq->size() < 1)
    {
      continue;
    }
    std::vector<int> vtrack_index;
    std::vector<int> vcalo_index;
    vtrack_index.clear();
    vcalo_index.clear();
    for(unsigned int itrack = 0; itrack < _track_ptq->size(); itrack++)
    {
      // cut good tpc tracks
      if(_track_nc_tpc->at(itrack) < 30) continue;
      float track_p = sqrt(pow(_track_px->at(itrack),2) + pow(_track_py->at(itrack),2) + pow(_track_pz->at(itrack),2));
      //if(track_p < 1.5) continue;

      if(isnan(_track_phi_emc->at(itrack)))
      {
        continue;
      }

      // project tpc tracks to emcal and hcal
      std::pair<float, float> TrackProjsEMCal;
      TrackProjsEMCal = std::make_pair( cal_phi(_track_x_emc->at(itrack), _track_y_emc->at(itrack)), _track_z_emc->at(itrack));

      bool is_match = false;
      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
        if (_emcal_e->at(iem)<0.5)
        {
          continue;
        }
        std::pair<float, float> EMCalPos;
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        EMCalPos = std::make_pair(_emcal_phi->at(iem), radius_scale * _emcal_z->at(iem));

        float dphi = PiRange(TrackProjsEMCal.first - EMCalPos.first);
        float dz = TrackProjsEMCal.second - EMCalPos.second;
//cout<<"dphi = "<<dphi<<" , dz = "<<dz<<endl;
        if(fabs(dphi)<0.2 && fabs(dz)<5)
        {
          is_match = true;
          vcalo_index.push_back(iem);
        }
      }
      if (is_match)
      {
        vtrack_index.push_back(itrack);
      }
    }

    if (vtrack_index.size()<2)
    {
      continue;
    }
    cout<<"Has matched pairs in event "<<i<<endl;

    // Draw
    // loop over all tracks
    for(unsigned int itrack = 0; itrack < _track_ptq->size(); itrack++)
    {
      float trk_project_emcal_x = _track_x_emc->at(itrack);
      float trk_project_emcal_y = _track_y_emc->at(itrack);
      float trk_project_emcal_z = _track_z_emc->at(itrack);
      float trk_project_emcal_R = sqrt(pow(trk_project_emcal_x,2)+pow(trk_project_emcal_y,2));
      float trk_project_hcal_x = _track_x_ohc->at(itrack);
      float trk_project_hcal_y = _track_y_ohc->at(itrack);
      float trk_project_hcal_z = _track_z_ohc->at(itrack);
      float trk_project_hcal_R = sqrt(pow(trk_project_hcal_x,2)+pow(trk_project_hcal_y,2));
      float trk_pt = sqrt(pow(_track_px->at(itrack),2)+pow(_track_py->at(itrack),2));
      float trk_pz = _track_pz->at(itrack);
      float trk_phi = _track_phi->at(itrack);
      float trk_x = _track_x_emc->at(itrack);
      float trk_y = _track_y_emc->at(itrack);
      float trk_z = _track_z_emc->at(itrack);
      int trk_charge = _track_ptq->at(itrack) > 0 ? 1 : -1;
cout<<"trk_project_emcal x = "<<trk_project_emcal_x<<" , y = "<<trk_project_emcal_y<<" , z = "<<trk_project_emcal_z<<endl;

      auto it = std::find(vtrack_index.begin(), vtrack_index.end(), itrack);

      if (it != vtrack_index.end()) {

        vec_track_pt.push_back(trk_pt);
        vec_track_pz.push_back(trk_pz);
        vec_track_phi.push_back(trk_phi);
        vec_track_z.push_back(trk_z);
        vec_track_charge.push_back(trk_charge);

        for(unsigned int ic = 0; ic < _trClus_track_id->size(); ic++)
        {
          if(_track_id->at(itrack) != _trClus_track_id->at(ic)) continue;
          vec_matched_track_cluster_x.push_back(_trClus_x->at(ic));
          vec_matched_track_cluster_y.push_back(_trClus_y->at(ic));
          vec_matched_track_cluster_z.push_back(_trClus_z->at(ic));
          vec_matched_track_cluster_R.push_back(sqrt(pow(_trClus_x->at(ic),2)+pow(_trClus_y->at(ic),2)));
        }
        vec_matched_track_project_emcal_x.push_back(trk_project_emcal_x);
        vec_matched_track_project_emcal_y.push_back(trk_project_emcal_y);
        vec_matched_track_project_emcal_z.push_back(trk_project_emcal_z);
        vec_matched_track_project_emcal_R.push_back(trk_project_emcal_R);
        vec_matched_track_project_hcal_x.push_back(trk_project_hcal_x);
        vec_matched_track_project_hcal_y.push_back(trk_project_hcal_y);
        vec_matched_track_project_hcal_z.push_back(trk_project_hcal_z);
        vec_matched_track_project_hcal_R.push_back(trk_project_hcal_R);

      } else {

        for(unsigned int ic = 0; ic < _trClus_track_id->size(); ic++)
        {
          if(_track_id->at(itrack) != _trClus_track_id->at(ic)) continue;
          vec_track_cluster_x.push_back(_trClus_x->at(ic));
          vec_track_cluster_y.push_back(_trClus_y->at(ic));
          vec_track_cluster_z.push_back(_trClus_z->at(ic));
          vec_track_cluster_R.push_back(sqrt(pow(_trClus_x->at(ic),2)+pow(_trClus_y->at(ic),2)));
        }
        vec_track_project_emcal_x.push_back(trk_project_emcal_x);
        vec_track_project_emcal_y.push_back(trk_project_emcal_y);
        vec_track_project_emcal_z.push_back(trk_project_emcal_z);
        vec_track_project_emcal_R.push_back(trk_project_emcal_R);
        vec_track_project_hcal_x.push_back(trk_project_hcal_x);
        vec_track_project_hcal_y.push_back(trk_project_hcal_y);
        vec_track_project_hcal_z.push_back(trk_project_hcal_z);
        vec_track_project_hcal_R.push_back(trk_project_hcal_R);

      }

    }

    if (*(std::min_element(vec_track_pt.begin(), vec_track_pt.end())) < minpt) minpt = *(std::min_element(vec_track_pt.begin(), vec_track_pt.end()));
    if (*(std::max_element(vec_track_pt.begin(), vec_track_pt.end())) > maxpt) maxpt = *(std::max_element(vec_track_pt.begin(), vec_track_pt.end()));

    float min_diff_dptdpz = 1000;
    float avg_diff_dptdpz = -1;
    for (int i1 = 0; i1 < (vec_track_pt.size()-1); i1++)
    {
      for (int i2 = i1+1; i2 < (vec_track_pt.size()); i2++)
      {
        if (vec_track_charge.at(i1) > 0 && vec_track_charge.at(i2) > 0) ncharge_pp++;
        if (vec_track_charge.at(i1) < 0 && vec_track_charge.at(i2) < 0) ncharge_nn++;
        if (vec_track_charge.at(i1) * vec_track_charge.at(i2) == -1) ncharge_pn++;
        if (vec_track_charge.at(i1) * vec_track_charge.at(i2) != -1) continue;

        float dptdpz1 = vec_track_pt.at(i1) / vec_track_pz.at(i1);
        float dptdpz2 = vec_track_pt.at(i2) / vec_track_pz.at(i2);
        float diff_tmp = fabs(dptdpz1 - dptdpz2);
        float avg_tmp = (dptdpz1+dptdpz2) / 2.;
        float ratio_pt = vec_track_pt.at(i1) / vec_track_pt.at(i2);
        float diff_phi_tmp = fabs(vec_track_phi.at(i1) - vec_track_phi.at(i2));
        float p1 = sqrt(vec_track_pt.at(i1)*vec_track_pt.at(i1) + vec_track_pz.at(i1)*vec_track_pz.at(i1));
        float p2 = sqrt(vec_track_pt.at(i2)*vec_track_pt.at(i2) + vec_track_pz.at(i2)*vec_track_pz.at(i2));
        if (diff_tmp < min_diff_dptdpz)
        {
          min_diff_dptdpz = diff_tmp;
          avg_diff_dptdpz = avg_tmp;
        }

        //if (diff_tmp >= 0.015) continue;
        //if (ratio_pt > 1.2 || ratio_pt < 0.8) continue;
        diff_phi.push_back(diff_phi_tmp);

        for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
        {
          if (_emcal_e->at(iem)<0.5)
          {
            continue;
          }
          std::pair<float, float> EMCalPos;
          float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
          EMCalPos = std::make_pair(_emcal_phi->at(iem), radius_scale * _emcal_z->at(iem));
cout<<"emcal x = "<<radius_scale * _emcal_x->at(iem)<<" , y = "<<radius_scale * _emcal_y->at(iem)<<" , z = "<<radius_scale * _emcal_z->at(iem)<<endl;

          float dphi1 = PiRange(vec_track_phi.at(i1) - EMCalPos.first);
          float dz1 = vec_track_z.at(i1) - EMCalPos.second;
          if(fabs(dphi1)<0.2 && fabs(dz1)<5)
          {
            matched_track_emcal_pair_E.push_back(_emcal_e->at(iem));
            matched_track_emcal_pair_p.push_back(p1);
          }

          float dphi2 = PiRange(vec_track_phi.at(i2) - EMCalPos.first);
          float dz2 = vec_track_z.at(i2) - EMCalPos.second;
          if(fabs(dphi2)<0.2 && fabs(dz2)<5)
          {
            matched_track_emcal_pair_E.push_back(_emcal_e->at(iem));
            matched_track_emcal_pair_p.push_back(p2);
          }

        }

      }
    }

    // loop all emcal clusters to match with tpc tracks
    for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
    {
      float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
      float x = _emcal_x->at(iem) * radius_scale;
      float y = _emcal_y->at(iem) * radius_scale;
      float z = _emcal_z->at(iem) * radius_scale;
      float R = sqrt(x*x+y*y);

      auto it = std::find(vcalo_index.begin(), vcalo_index.end(), iem);

      if (it != vcalo_index.end()) {
        vec_matched_emcal_x.push_back(x);
        vec_matched_emcal_y.push_back(y);
        vec_matched_emcal_z.push_back(z);
        vec_matched_emcal_R.push_back(R);
      } else {
        vec_emcal_x.push_back(x);
        vec_emcal_y.push_back(y);
        vec_emcal_z.push_back(z);
        vec_emcal_R.push_back(R);
      }

    }

    // loop all tpc clusters no mater it belongs to tracks or not
    for(unsigned int iclu = 0; iclu < _cluster_x->size(); iclu++)
    {
      float x = _cluster_x->at(iclu);
      float y = _cluster_y->at(iclu);
      float z = _cluster_z->at(iclu);
      float R = sqrt(x*x+y*y);

      vec_cluster_x.push_back(x);
      vec_cluster_y.push_back(y);
      vec_cluster_z.push_back(z);
      vec_cluster_R.push_back(R);
    }

    all_tpc_cluster_x.push_back(vec_cluster_x);
    all_tpc_cluster_y.push_back(vec_cluster_y);
    all_tpc_cluster_z.push_back(vec_cluster_z);
    all_tpc_cluster_R.push_back(vec_cluster_R);
    unmatched_tpc_cluster_x.push_back(vec_track_cluster_x);
    unmatched_tpc_cluster_y.push_back(vec_track_cluster_y);
    unmatched_tpc_cluster_z.push_back(vec_track_cluster_z);
    unmatched_tpc_cluster_R.push_back(vec_track_cluster_R);
for (int ii=0; ii<(vec_emcal_x.size()); ii++)
{
  cout<<"emcal x = "<<vec_emcal_x.at(ii)<<" , y = "<<vec_emcal_y.at(ii)<<" , z = "<<vec_emcal_z.at(ii)<<endl;
}
    unmatched_emcal_cluster_x.push_back(vec_emcal_x);
    unmatched_emcal_cluster_y.push_back(vec_emcal_y);
    unmatched_emcal_cluster_z.push_back(vec_emcal_z);
    unmatched_emcal_cluster_R.push_back(vec_emcal_R);
    unmatched_track_project_emcal_x.push_back(vec_track_project_emcal_x);
    unmatched_track_project_emcal_y.push_back(vec_track_project_emcal_y);
    unmatched_track_project_emcal_z.push_back(vec_track_project_emcal_z);
    unmatched_track_project_emcal_R.push_back(vec_track_project_emcal_R);
    matched_tpc_cluster_x.push_back(vec_matched_track_cluster_x);
    matched_tpc_cluster_y.push_back(vec_matched_track_cluster_y);
    matched_tpc_cluster_z.push_back(vec_matched_track_cluster_z);
    matched_tpc_cluster_R.push_back(vec_matched_track_cluster_R);
    matched_emcal_cluster_x.push_back(vec_matched_emcal_x);
    matched_emcal_cluster_y.push_back(vec_matched_emcal_y);
    matched_emcal_cluster_z.push_back(vec_matched_emcal_z);
    matched_emcal_cluster_R.push_back(vec_matched_emcal_R);
    matched_track_project_emcal_x.push_back(vec_matched_track_project_emcal_x);
    matched_track_project_emcal_y.push_back(vec_matched_track_project_emcal_y);
    matched_track_project_emcal_z.push_back(vec_matched_track_project_emcal_z);
    matched_track_project_emcal_R.push_back(vec_matched_track_project_emcal_R);
    eventid.push_back(_eventNumber);
    diff_dptdpz.push_back(min_diff_dptdpz);
    avg_dptdpz.push_back(avg_diff_dptdpz);

  }
  cout<<"ncharge_pp = "<<ncharge_pp<<" , ncharge_pn = "<<ncharge_pn<<" , ncharge_nn = "<<ncharge_nn<<endl;
  std::cout << "min pt = " << minpt << std::endl;
  std::cout << "max pt = " << maxpt << std::endl;

  cout<<"number of matched event = "<<matched_tpc_cluster_x.size()<<endl;

  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.05);
  canvas->cd();

  //TH1F* h_diff_dptdpz = new TH1F("h_diff_dptdpz","h_diff_dptdpz",100, 0, 6);
  TH1F* h_diff_dptdpz = new TH1F("h_diff_dptdpz","h_diff_dptdpz",100, 0, 0.1);
  for (int i = 0; i < (diff_dptdpz.size()); i++)
  {
    h_diff_dptdpz->Fill(diff_dptdpz.at(i));
  }

  h_diff_dptdpz->SetTitle(Form("sPHENIX Internal, Run %d;|#Delta(pT/pz)|;Entries",runnumber));
  h_diff_dptdpz->SetMinimum(0.1);
  h_diff_dptdpz->Draw("hist");

  canvas->SetLogy(1);
  canvas->Update();
  canvas->SaveAs(Form("figure/Diff_dptdpz_run%d.pdf",runnumber));

  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 800, 800);
  canvas1->SetLeftMargin(0.15);
  canvas1->SetRightMargin(0.05);
  canvas1->cd();

  TH1F* h_diff_phi = new TH1F("h_diff_phi","h_diff_phi",100, 0, M_PI);
  for (int i = 0; i < (diff_phi.size()); i++)
  {
    h_diff_phi->Fill(diff_phi.at(i));
  }

  h_diff_phi->SetTitle(Form("sPHENIX Internal, Run %d;|#Phi_{1}-#Phi_{2}|;Entries",runnumber));
  h_diff_phi->SetMinimum(0.1);
  h_diff_phi->Draw("hist");

  canvas1->Update();
  canvas1->SaveAs(Form("figure/Diff_phi_run%d.pdf",runnumber));

  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 800, 800);
  canvas2->SetLeftMargin(0.15);
  canvas2->SetRightMargin(0.05);
  canvas2->cd();

  TH2F* h2_eop_p = new TH2F("h2_eop_p","h2_eop_p", 20, 0, 10, 20, 0, 2);
  for (int i = 0; i < (matched_track_emcal_pair_E.size()); i++)
  {
    float eop = matched_track_emcal_pair_E.at(i) / matched_track_emcal_pair_p.at(i);
    float p = matched_track_emcal_pair_p.at(i);
    h2_eop_p->Fill(p, eop);
  }

  h2_eop_p->SetTitle(Form("sPHENIX Internal, Run %d;P [GeV/c]; E/p;Entries",runnumber));
  h2_eop_p->Draw("colz");

  canvas2->Update();
  canvas2->SaveAs(Form("figure/eOp_p_run%d.pdf",runnumber));

  // draw x y plots
  int nevent = matched_tpc_cluster_x.size();
  for (int index = 0; index<nevent; index++)
  {
    TCanvas *can = new TCanvas(Form("can_%d",index), "can", 800, 800);
    can->SetLeftMargin(0.15);
    can->SetRightMargin(0.05);
    can->cd();

    TGraph *gr_tpc_all = new TGraph(all_tpc_cluster_x.at(index).size(), all_tpc_cluster_x.at(index).data(), all_tpc_cluster_y.at(index).data());
    gr_tpc_all->SetTitle(Form("sPHENIX Internal, Run %d Event %d;X [cm];Y [cm]",runnumber,eventid.at(index)));
    gr_tpc_all->GetXaxis()->SetLimits(-(emcal_radius+10.), (emcal_radius+10.));
    gr_tpc_all->SetMinimum(-(emcal_radius+10.));
    gr_tpc_all->SetMaximum((emcal_radius+10.));
    gr_tpc_all->SetMarkerSize(0.2);
    gr_tpc_all->SetMarkerStyle(20);
    gr_tpc_all->SetMarkerColor(kGray);
    gr_tpc_all->Draw("AP");

    TGraph *gr_tpc_unmatched = new TGraph(unmatched_tpc_cluster_x.at(index).size(), unmatched_tpc_cluster_x.at(index).data(), unmatched_tpc_cluster_y.at(index).data());
    gr_tpc_unmatched->SetTitle(Form("sPHENIX Internal, Run %d Event %d;X [cm];Y [cm]",runnumber,eventid.at(index)));
    gr_tpc_unmatched->GetXaxis()->SetLimits(-(emcal_radius+10.), (emcal_radius+10.));
    gr_tpc_unmatched->SetMinimum(-(emcal_radius+10.));
    gr_tpc_unmatched->SetMaximum((emcal_radius+10.));
    gr_tpc_unmatched->SetMarkerSize(0.2);
    gr_tpc_unmatched->SetMarkerStyle(20);
    gr_tpc_unmatched->SetMarkerColor(kBlue);
    gr_tpc_unmatched->Draw("P,same");

    TGraph *gr_tpc_matched = new TGraph(matched_tpc_cluster_x.at(index).size(), matched_tpc_cluster_x.at(index).data(), matched_tpc_cluster_y.at(index).data());
    gr_tpc_matched->SetMarkerStyle(20);
    gr_tpc_matched->SetMarkerSize(0.2);
    gr_tpc_matched->SetMarkerColor(kRed);
    gr_tpc_matched->Draw("P,same");

    TGraph *gr_tpc_proj_emcal_unmatched = new TGraph(unmatched_track_project_emcal_x.at(index).size(), unmatched_track_project_emcal_x.at(index).data(), unmatched_track_project_emcal_y.at(index).data());
    gr_tpc_proj_emcal_unmatched->SetMarkerColor(kBlue);
    gr_tpc_proj_emcal_unmatched->Draw("P*,same");

    TGraph *gr_tpc_proj_emcal_matched = new TGraph(matched_track_project_emcal_x.at(index).size(), matched_track_project_emcal_x.at(index).data(), matched_track_project_emcal_y.at(index).data());
    gr_tpc_proj_emcal_matched->SetMarkerColor(kRed);
    gr_tpc_proj_emcal_matched->Draw("P*,same");

    TGraph *gr_emcal_unmatched = new TGraph(unmatched_emcal_cluster_x.at(index).size(), unmatched_emcal_cluster_x.at(index).data(), unmatched_emcal_cluster_y.at(index).data());
    gr_emcal_unmatched->SetMarkerColor(kBlue);
    gr_emcal_unmatched->SetMarkerStyle(22);
    gr_emcal_unmatched->Draw("P,same");
  
    TGraph *gr_emcal_matched = new TGraph(matched_emcal_cluster_x.at(index).size(), matched_emcal_cluster_x.at(index).data(), matched_emcal_cluster_y.at(index).data());
    gr_emcal_matched->SetMarkerColor(kRed);
    gr_emcal_matched->SetMarkerStyle(22);
    gr_emcal_matched->Draw("P,same");
 
    // Draw ellipses using TGraph
    DrawEllipseWithTGraph(0., 0., emcal_radius - 3., kBlack);
    DrawEllipseWithTGraph(0., 0., emcal_radius + 3., kBlack);

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

    can->SaveAs(Form("figure/EvtDisplay_xy_run%d_event%d.pdf",runnumber,index));

  }

  for (int index = 0; index<nevent; index++)
  {
    TCanvas *can = new TCanvas(Form("can_%d",index), "can", 800, 800);
    can->SetLeftMargin(0.15);
    can->SetRightMargin(0.05);
    can->cd();

    TGraph *gr_tpc_all = new TGraph(all_tpc_cluster_R.at(index).size(), all_tpc_cluster_z.at(index).data(), all_tpc_cluster_R.at(index).data());
    gr_tpc_all->SetTitle(Form("sPHENIX Internal, Run %d Event %d;Z [cm];R [cm]",runnumber,eventid.at(index)));
    gr_tpc_all->GetXaxis()->SetLimits(-150, 150);
    gr_tpc_all->SetMinimum(20);
    gr_tpc_all->SetMaximum(emcal_radius+10.);
    gr_tpc_all->SetMarkerSize(0.2);
    gr_tpc_all->SetMarkerStyle(20);
    gr_tpc_all->SetMarkerColor(kGray);
    gr_tpc_all->Draw("AP");

    TGraph *gr_tpc_unmatched = new TGraph(unmatched_tpc_cluster_R.at(index).size(), unmatched_tpc_cluster_z.at(index).data(), unmatched_tpc_cluster_R.at(index).data());
    gr_tpc_unmatched->SetTitle(Form("sPHENIX Internal, Run %d Event %d;Z [cm];R [cm]",runnumber,eventid.at(index)));
    gr_tpc_unmatched->GetXaxis()->SetLimits(-150, 150);
    gr_tpc_unmatched->SetMinimum(20);
    gr_tpc_unmatched->SetMaximum(emcal_radius+10.);
    gr_tpc_unmatched->SetMarkerStyle(20);
    gr_tpc_unmatched->SetMarkerSize(0.2);
    gr_tpc_unmatched->SetMarkerColor(kBlue);
    gr_tpc_unmatched->Draw("P,same");

    TGraph *gr_tpc_matched = new TGraph(matched_tpc_cluster_R.at(index).size(), matched_tpc_cluster_z.at(index).data(), matched_tpc_cluster_R.at(index).data());
    gr_tpc_matched->SetMarkerStyle(20);
    gr_tpc_matched->SetMarkerSize(0.2);
    gr_tpc_matched->SetMarkerColor(kRed);
    gr_tpc_matched->Draw("P,same");

    TGraph *gr_tpc_proj_emcal_unmatched = new TGraph(unmatched_track_project_emcal_R.at(index).size(), unmatched_track_project_emcal_z.at(index).data(), unmatched_track_project_emcal_R.at(index).data());
    gr_tpc_proj_emcal_unmatched->SetMarkerColor(kBlue);
    gr_tpc_proj_emcal_unmatched->Draw("P*,same");

    TGraph *gr_tpc_proj_emcal_matched = new TGraph(matched_track_project_emcal_R.at(index).size(), matched_track_project_emcal_z.at(index).data(), matched_track_project_emcal_R.at(index).data());
    gr_tpc_proj_emcal_matched->SetMarkerColor(kRed);
    gr_tpc_proj_emcal_matched->Draw("P*,same");

    TGraph *gr_emcal_unmatched = new TGraph(unmatched_emcal_cluster_R.at(index).size(), unmatched_emcal_cluster_z.at(index).data(), unmatched_emcal_cluster_R.at(index).data());
    gr_emcal_unmatched->SetMarkerColor(kBlue);
    gr_emcal_unmatched->SetMarkerStyle(22);
    gr_emcal_unmatched->Draw("P,same");

    TGraph *gr_emcal_matched = new TGraph(matched_emcal_cluster_R.at(index).size(), matched_emcal_cluster_z.at(index).data(), matched_emcal_cluster_R.at(index).data());
    gr_emcal_matched->SetMarkerColor(kRed);
    gr_emcal_matched->SetMarkerStyle(22);
    gr_emcal_matched->Draw("P,same");

    DrawHLineWithTGraph(emcal_radius - 3., -150, 150, kBlack);
    DrawHLineWithTGraph(emcal_radius + 3., -150, 150, kBlack);

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

    can->SaveAs(Form("figure/EvtDisplay_rz_run%d_event%d.pdf",runnumber,index));

  }

}

