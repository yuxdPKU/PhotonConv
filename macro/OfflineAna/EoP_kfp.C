#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void EoP_kfp(int runnumber)
{
  int verbosity = 0;

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree_KFP");
  //chain->Add(Form("./%d/TrackCalo_*_ana.root",runnumber));
  chain->Add(Form("eop_all_kfp.root"));

  setBranch_kfp(chain);

  TFile* outputfile = new TFile(Form("./eop_%d_kfp.root",runnumber),"recreate");
  TTree* outputtree = new TTree("tree","tree with eop info");

  float teop_gamma_mass, teop_gamma_radius;
  float teop_gamma_vertex_volume, teop_gamma_SV_chi2_per_nDoF;
  float teop_gamma_quality, teop_track_quality;
  float teop_emcal_e, teop_emcal_phi, teop_emcal_z, teop_emcal_eta;
  float teop_track_p, teop_track_phi_projemc, teop_track_z_projemc;
  float teop_track_pt, teop_track_eta, teop_track_e;
  float teop_track_pt_raw, teop_track_p_raw;
  float teop_track_pt_unmoved, teop_track_p_unmoved, teop_track_e_unmoved;
  int teop_track_charge, teop_track_crossing;
  float teop_track_mass_1, teop_track_mass_2, teop_track_mass_3;
  float teop_track12_deta, teop_track12_dphi, teop_track12_dz_emc;
  float teop_trkemc_dphi, teop_trkemc_dz;
  float teop_track12_dca_2d, teop_track12_dca_3d;

  outputtree->Branch("_runNumber",&_runNumber,"_runNumber/I");
  outputtree->Branch("_eventNumber",&_eventNumber,"_eventNumber/I");
  outputtree->Branch("_gamma_mass",&teop_gamma_mass,"_gamma_mass/F");
  outputtree->Branch("_gamma_radius",&teop_gamma_radius,"_gamma_radius/F");
  outputtree->Branch("_gamma_quality",&teop_gamma_quality,"_gamma_quality/F");
  outputtree->Branch("_gamma_vertex_volume",&teop_gamma_vertex_volume,"_gamma_vertex_volume/F");
  outputtree->Branch("_gamma_SV_chi2_per_nDoF",&teop_gamma_SV_chi2_per_nDoF,"_gamma_SV_chi2_per_nDoF/F");
  outputtree->Branch("_emcal_e",&teop_emcal_e,"_emcal_e/F");
  outputtree->Branch("_emcal_phi",&teop_emcal_phi,"_emcal_phi/F");
  outputtree->Branch("_emcal_z",&teop_emcal_z,"_emcal_z/F");
  outputtree->Branch("_emcal_eta",&teop_emcal_eta,"_emcal_eta/F");
  outputtree->Branch("_track_p",&teop_track_p,"_track_p/F");
  outputtree->Branch("_track_p_raw",&teop_track_p_raw,"_track_p_raw/F");
  outputtree->Branch("_track_p_unmoved",&teop_track_p_unmoved,"_track_p_unmoved/F");
  outputtree->Branch("_track_pt",&teop_track_pt,"_track_pt/F");
  outputtree->Branch("_track_pt_raw",&teop_track_pt_raw,"_track_pt_raw/F");
  outputtree->Branch("_track_pt_unmoved",&teop_track_pt_unmoved,"_track_pt_unmoved/F");
  outputtree->Branch("_track_eta",&teop_track_eta,"_track_eta/F");
  outputtree->Branch("_track_e",&teop_track_e,"_track_e/F");
  outputtree->Branch("_track_e_unmoved",&teop_track_e_unmoved,"_track_e_unmoved/F");
  outputtree->Branch("_track_mass_1",&teop_track_mass_1,"_track_mass_1/F");
  outputtree->Branch("_track_mass_2",&teop_track_mass_2,"_track_mass_2/F");
  outputtree->Branch("_track_mass_3",&teop_track_mass_3,"_track_mass_3/F");
  outputtree->Branch("_track_phi_projemc",&teop_track_phi_projemc,"_track_phi_projemc/F");
  outputtree->Branch("_track_z_projemc",&teop_track_z_projemc,"_track_z_projemc/F");
  outputtree->Branch("_track_quality",&teop_track_quality,"_track_quality/F");
  outputtree->Branch("_track_charge",&teop_track_charge,"_track_charge/I");
  outputtree->Branch("_track_crossing",&teop_track_crossing,"_track_crossing/I");
  outputtree->Branch("_track12_deta",&teop_track12_deta,"_track12_deta/F");
  outputtree->Branch("_track12_dphi",&teop_track12_dphi,"_track12_dphi/F");
  outputtree->Branch("_track12_dz_emc",&teop_track12_dz_emc,"_track12_dz_emc/F");
  outputtree->Branch("_trkemc_dphi",&teop_trkemc_dphi,"_trkemc_dphi/F");
  outputtree->Branch("_trkemc_dz",&teop_trkemc_dz,"_trkemc_dz/F");
  outputtree->Branch("_track12_dca_2d",&teop_track12_dca_2d,"_track12_dca_2d/F");
  outputtree->Branch("_track12_dca_3d",&teop_track12_dca_3d,"_track12_dca_3d/F");

  std::vector<float> vec_emcal_matched_e;
  std::vector<float> vec_emcal_matched_phi;
  std::vector<float> vec_emcal_matched_z;
  std::vector<float> vec_emcal_matched_eta;
  std::vector<float> vec_track_emcal_residual_phi;
  std::vector<float> vec_track_emcal_residual_z;

  float max_e = 0;
  float min_dphi = 100;
  int index = -1;

  int nevent  = chain->GetEntries();
  cout<<"total nevent = "<<nevent<<endl;
  for(int i = 0; i < nevent; i++)
  {
    if (i % (nevent / 10) == 0) cout << "Processing progress: " << i / (nevent / 10) << "0%" << endl;
    chain->GetEntry(i);

    for (int ican = 0; ican < _numCan; ican++)
    {
      if (_em_pT->at(ican) < min_Track_Pt || _ep_pT->at(ican) < min_Track_Pt) continue;

      float quality_cut = 1000;
      if ( (_gamma_chi2->at(ican) / _gamma_nDoF->at(ican))>quality_cut || (_ep_chi2->at(ican) / _ep_nDoF->at(ican))>quality_cut || (_em_chi2->at(ican) / _em_nDoF->at(ican))>quality_cut)
      {
        if (verbosity>0) std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail gamma/e+/e- quality cut"<<std::endl;
        if (verbosity>0) std::cout<<"gamma quality = "<<_gamma_chi2->at(ican) / _gamma_nDoF->at(ican)<<" e+ quality = "<<_ep_chi2->at(ican) / _ep_nDoF->at(ican)<<" e- quality = "<<_em_chi2->at(ican) / _em_nDoF->at(ican)<<std::endl;
        //continue;
      }

      float gamma_r_cut = 0;
      float gamma_r = customsqrt(_gamma_x->at(ican)*_gamma_x->at(ican)+_gamma_y->at(ican)*_gamma_y->at(ican));
      if (gamma_r < gamma_r_cut)
      {  
        if (verbosity>0) std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail secondary vertex cut"<<std::endl;
        //continue;
      }

      float gamma_mass_cut = 1;
      if (_gamma_mass->at(ican) > gamma_mass_cut)
      {  
        if (verbosity>0) std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail mass cut"<<std::endl;
        //continue;
      }

      teop_gamma_mass = _gamma_mass->at(ican);
      teop_gamma_radius = gamma_r;
      teop_gamma_quality = _gamma_chi2->at(ican) / _gamma_nDoF->at(ican);
      teop_gamma_vertex_volume = _gamma_vertex_volume->at(ican);
      teop_gamma_SV_chi2_per_nDoF = _gamma_SV_chi2_per_nDoF->at(ican);

      // match e+ with calo
      vec_emcal_matched_e.clear();
      vec_emcal_matched_phi.clear();
      vec_emcal_matched_z.clear();
      vec_emcal_matched_eta.clear();
      vec_track_emcal_residual_phi.clear();
      vec_track_emcal_residual_z.clear();
      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
        if (verbosity>0)
        {
          cout<<"emcal "<<iem<<" e = "<<_emcal_e->at(iem)<<" x = "<<_emcal_x->at(iem)<<" y = "<<_emcal_y->at(iem)<<" z = "<<_emcal_z->at(iem)<<endl;
        }
        // cut good emcal clusters
        if(_emcal_e->at(iem) < min_EMCal_E) continue;
        //if(fabs(_emcal_eta->at(iem)) > 1.1) continue;
        float emcal_phi = cal_phi(_emcal_x->at(iem), _emcal_y->at(iem));
        float track_phi = cal_phi(_ep_x_emc->at(ican), _ep_y_emc->at(ican));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_z = radius_scale * _emcal_z->at(iem);
        float emcal_eta = cal_eta(_emcal_x->at(iem), _emcal_y->at(iem), _emcal_z->at(iem));

        float dphi = PiRange(track_phi - emcal_phi);
        float dz = _ep_z_emc->at(ican) - emcal_z;

        //matching
        //if (dphi>-0.02 && dphi<0.1 && fabs(dz)<20) // loose cut in dz
        //if (dphi>-0.1 && dphi<0.1 && fabs(dz)<20) // loose cut in all
        //if (dphi>0 && dphi<0.1 && fabs(dz)<5) // tight cut in all
        if (dphi>-0.2 && dphi<0.2) // no cut in all
        {
          vec_emcal_matched_e.push_back(_emcal_e->at(iem));
          vec_emcal_matched_phi.push_back(emcal_phi);
          vec_emcal_matched_z.push_back(emcal_z);
          vec_emcal_matched_eta.push_back(emcal_eta);
          vec_track_emcal_residual_phi.push_back(dphi);
          vec_track_emcal_residual_z.push_back(dz);
        }

      }

      if (vec_emcal_matched_e.size() == 0)
      {
        if (verbosity>0)cout<<"no matched???"<<endl;
        continue;
      }

      max_e = 0;
      min_dphi = 100;
      index = -1;
      for(unsigned int iem = 0; iem < vec_track_emcal_residual_phi.size(); iem++)
      {
        if (fabs(vec_track_emcal_residual_phi.at(iem))<min_dphi)
        {
          min_dphi = fabs(vec_track_emcal_residual_phi.at(iem));
          index = iem;
        }
      }

      teop_emcal_e = vec_emcal_matched_e.at(index);
      teop_emcal_phi = vec_emcal_matched_phi.at(index);
      teop_emcal_z = vec_emcal_matched_z.at(index);
      teop_emcal_eta = vec_emcal_matched_eta.at(index);
      teop_track_p = _ep_p->at(ican);
      teop_track_p_raw = _ep_p_raw->at(ican);
      teop_track_p_unmoved = _ep_p_unmoved->at(ican);
      teop_track_pt = _ep_pT->at(ican);
      teop_track_pt_raw = _ep_pT_raw->at(ican);
      teop_track_pt_unmoved = _ep_pT_unmoved->at(ican);
      teop_track_e_unmoved = _ep_pE_unmoved->at(ican);
      teop_track_eta = _ep_pseudorapidity->at(ican);
      teop_track_mass_1 = _ep_mass->at(ican);
      teop_track_mass_2 = customsqrt( pow(_ep_pE->at(ican),2) - pow(_ep_p->at(ican),2) );
      teop_track_mass_3 = customsqrt( pow(_ep_pE_unmoved->at(ican),2) - pow(_ep_p_unmoved->at(ican),2) );
      teop_track_phi_projemc = cal_phi(_ep_x_emc->at(ican), _ep_y_emc->at(ican));
      teop_track_z_projemc = _ep_z_emc->at(ican);
      teop_track_quality = _ep_chi2->at(ican) / _ep_nDoF->at(ican);
      teop_track_charge = 1;
      teop_track_crossing = _ep_crossing->at(ican);
      teop_track12_deta = _ep_pseudorapidity->at(ican) - _em_pseudorapidity->at(ican);
      teop_track12_dphi = _ep_phi->at(ican) - _em_phi->at(ican);
      teop_track12_dz_emc = _ep_z_emc->at(ican) - _em_z_emc->at(ican);
      teop_trkemc_dphi = vec_track_emcal_residual_phi.at(index);
      teop_trkemc_dz = vec_track_emcal_residual_z.at(index);
      teop_track12_dca_2d = _epem_DCA_2d->at(ican);
      teop_track12_dca_3d = _epem_DCA_3d->at(ican);

      outputtree->Fill();
      if (verbosity>0) cout<<"fill e+"<<endl;

      // match e- with calo
      vec_emcal_matched_e.clear();
      vec_emcal_matched_phi.clear();
      vec_emcal_matched_z.clear();
      vec_emcal_matched_eta.clear();
      vec_track_emcal_residual_phi.clear();
      vec_track_emcal_residual_z.clear();
      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
        if (verbosity>0) cout<<"emcal "<<iem<<" e = "<<_emcal_e->at(iem)<<" x = "<<_emcal_x->at(iem)<<" y = "<<_emcal_y->at(iem)<<" z = "<<_emcal_z->at(iem)<<endl;
        // cut good emcal clusters
        if(_emcal_e->at(iem) < min_EMCal_E) continue;
        //if(fabs(_emcal_eta->at(iem)) > 1.1) continue;
        float emcal_phi = cal_phi(_emcal_x->at(iem), _emcal_y->at(iem));
        float track_phi = cal_phi(_em_x_emc->at(ican), _em_y_emc->at(ican));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_z = radius_scale * _emcal_z->at(iem);
        float emcal_eta = cal_eta(_emcal_x->at(iem), _emcal_y->at(iem), _emcal_z->at(iem));

        float dphi = PiRange(track_phi - emcal_phi);
        float dz = _em_z_emc->at(ican) - emcal_z;

        //matching
        //if (dphi>-0.02 && dphi<0.1 && fabs(dz)<20) // loose cut in dz
        //if (dphi>-0.1 && dphi<0.1 && fabs(dz)<20) // loose cut in all
        //if (dphi>0 && dphi<0.1 && fabs(dz)<5) // tight cut in all
        if (dphi>-0.2 && dphi<0.2) // no cut in all
        {
          vec_emcal_matched_e.push_back(_emcal_e->at(iem));
          vec_emcal_matched_phi.push_back(emcal_phi);
          vec_emcal_matched_z.push_back(emcal_z);
          vec_emcal_matched_eta.push_back(emcal_eta);
          vec_track_emcal_residual_phi.push_back(dphi);
          vec_track_emcal_residual_z.push_back(dz);
        }

      }

      if (vec_emcal_matched_e.size() == 0)
      {
        if (verbosity>0)cout<<"no matched???"<<endl;
        continue;
      }

      max_e = 0;
      min_dphi = 100;
      index = -1;
      for(unsigned int iem = 0; iem < vec_track_emcal_residual_phi.size(); iem++)
      {
        if (fabs(vec_track_emcal_residual_phi.at(iem))<min_dphi)
        {
          min_dphi = fabs(vec_track_emcal_residual_phi.at(iem));
          index = iem;
        }
      }

      teop_emcal_e = vec_emcal_matched_e.at(index);
      teop_emcal_phi = vec_emcal_matched_phi.at(index);
      teop_emcal_z = vec_emcal_matched_z.at(index);
      teop_emcal_eta = vec_emcal_matched_eta.at(index);
      teop_track_p = _em_p->at(ican);
      teop_track_p_raw = _em_p_raw->at(ican);
      teop_track_p_unmoved = _em_p_unmoved->at(ican);
      teop_track_pt = _em_pT->at(ican);
      teop_track_pt_raw = _em_pT_raw->at(ican);
      teop_track_pt_unmoved = _em_pT_unmoved->at(ican);
      teop_track_e_unmoved = _em_pE_unmoved->at(ican);
      teop_track_eta = _em_pseudorapidity->at(ican);
      teop_track_mass_1 = _em_mass->at(ican);
      teop_track_mass_2 = customsqrt( pow(_em_pE->at(ican),2) - pow(_em_p->at(ican),2) );
      teop_track_mass_3 = customsqrt( pow(_em_pE_unmoved->at(ican),2) - pow(_em_p_unmoved->at(ican),2) );
      teop_track_phi_projemc = cal_phi(_em_x_emc->at(ican), _em_y_emc->at(ican));
      teop_track_z_projemc = _em_z_emc->at(ican);
      teop_track_quality = _em_chi2->at(ican) / _em_nDoF->at(ican);
      teop_track_charge = -1;
      teop_track_crossing = _em_crossing->at(ican);
      teop_track12_deta = _ep_pseudorapidity->at(ican) - _em_pseudorapidity->at(ican);
      teop_track12_dphi = _ep_phi->at(ican) - _em_phi->at(ican);
      teop_track12_dz_emc = _ep_z_emc->at(ican) - _em_z_emc->at(ican);
      teop_trkemc_dphi = vec_track_emcal_residual_phi.at(index);
      teop_trkemc_dz = vec_track_emcal_residual_z.at(index);
      teop_track12_dca_2d = _epem_DCA_2d->at(ican);
      teop_track12_dca_3d = _epem_DCA_3d->at(ican);

      outputtree->Fill();
      if (verbosity>0) cout<<"fill e-"<<endl;
    }
  }

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();
}
