#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void EoP_kfp_v2(int runnumber=0)
{
  int verbosity = 0;

  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree_KFP");
  //chain->Add(Form("./%d/TrackCalo_*_ana.root",runnumber));
  chain->Add(Form("eop_all_kfp.root"));

  setBranch_kfp(chain);

  TFile* outputfile = new TFile(Form("./eop_%d_kfp_v2.root",runnumber),"recreate");
  TTree* outputtree = new TTree("tree","tree with eop info");

  float teop_gamma_mass, teop_gamma_radius;
  float teop_gamma_x, teop_gamma_y, teop_gamma_z;
  float teop_gamma_pE;
  float teop_gamma_vertex_volume, teop_gamma_SV_chi2_per_nDoF;
  float teop_gamma_quality, teop_ep_quality, teop_em_quality;
  float teop_ep_emcal_e, teop_ep_emcal_phi, teop_ep_emcal_x, teop_ep_emcal_y, teop_ep_emcal_z, teop_ep_emcal_eta;
  float teop_em_emcal_e, teop_em_emcal_phi, teop_em_emcal_x, teop_em_emcal_y, teop_em_emcal_z, teop_em_emcal_eta;
  float teop_ep_p, teop_ep_phi_projemc, teop_ep_x_projemc, teop_ep_y_projemc, teop_ep_z_projemc;
  float teop_em_p, teop_em_phi_projemc, teop_em_x_projemc, teop_em_y_projemc, teop_em_z_projemc;
  float teop_ep_pt, teop_ep_eta, teop_ep_phi, teop_ep_pE;
  float teop_em_pt, teop_em_eta, teop_em_phi, teop_em_pE;
  float teop_ep_pt_raw, teop_ep_p_raw;
  float teop_em_pt_raw, teop_em_p_raw;
  float teop_ep_pt_unmoved, teop_ep_p_unmoved, teop_ep_pE_unmoved;
  float teop_em_pt_unmoved, teop_em_p_unmoved, teop_em_pE_unmoved;
  int teop_ep_charge, teop_ep_crossing;
  int teop_em_charge, teop_em_crossing;
  float teop_ep_mass_1, teop_ep_mass_2, teop_ep_mass_3;
  float teop_em_mass_1, teop_em_mass_2, teop_em_mass_3;
  float teop_track12_deta, teop_track12_dphi, teop_track12_dz_emc, teop_track12_dp;
  float teop_epemc_dphi, teop_epemc_dz;
  float teop_ememc_dphi, teop_ememc_dz;
  float teop_track12_dca_2d, teop_track12_dca_3d;
  float teop_eop_ep, teop_eop_em;
  std::vector<float> teop_ep_clus_x, teop_ep_clus_y, teop_ep_clus_z;
  std::vector<float> teop_em_clus_x, teop_em_clus_y, teop_em_clus_z;
  std::vector<float> teop_true_gamma_mass, teop_true_gamma_pE, teop_true_gamma_eta;
  std::vector<float> teop_true_gamma_prod_radius, teop_true_gamma_prod_x, teop_true_gamma_prod_y, teop_true_gamma_prod_z;
  std::vector<float> teop_true_gamma_decay_radius, teop_true_gamma_decay_x, teop_true_gamma_decay_y, teop_true_gamma_decay_z;
  std::vector<int> teop_true_gamma_mother_id;
  std::vector<float> teop_true_ep_pt, teop_true_ep_p, teop_true_ep_eta, teop_true_ep_phi, teop_true_ep_pE;
  std::vector<float> teop_true_em_pt, teop_true_em_p, teop_true_em_eta, teop_true_em_phi, teop_true_em_pE;
  //for truth-reco matching test
  std::vector<float> teop_ep_phi_recotruthDiff, teop_em_phi_recotruthDiff;
  std::vector<float> teop_ep_eta_recotruthDiff, teop_em_eta_recotruthDiff;
  std::vector<float> teop_SV_phi_recotruthDiff, teop_SV_z_recotruthDiff;

  outputtree->Branch("_runNumber",&_runNumber,"_runNumber/I");
  outputtree->Branch("_eventNumber",&_eventNumber,"_eventNumber/I");
  outputtree->Branch("_gamma_mass",&teop_gamma_mass,"_gamma_mass/F");
  outputtree->Branch("_gamma_radius",&teop_gamma_radius,"_gamma_radius/F");
  outputtree->Branch("_gamma_x",&teop_gamma_x,"_gamma_x/F");
  outputtree->Branch("_gamma_y",&teop_gamma_y,"_gamma_y/F");
  outputtree->Branch("_gamma_z",&teop_gamma_z,"_gamma_z/F");
  outputtree->Branch("_gamma_pE",&teop_gamma_pE,"_gamma_pE/F");
  outputtree->Branch("_gamma_quality",&teop_gamma_quality,"_gamma_quality/F");
  outputtree->Branch("_gamma_vertex_volume",&teop_gamma_vertex_volume,"_gamma_vertex_volume/F");
  outputtree->Branch("_gamma_SV_chi2_per_nDoF",&teop_gamma_SV_chi2_per_nDoF,"_gamma_SV_chi2_per_nDoF/F");
  outputtree->Branch("_ep_emcal_e",&teop_ep_emcal_e,"_ep_emcal_e/F");
  outputtree->Branch("_em_emcal_e",&teop_em_emcal_e,"_em_emcal_e/F");
  outputtree->Branch("_ep_emcal_phi",&teop_ep_emcal_phi,"_ep_emcal_phi/F");
  outputtree->Branch("_em_emcal_phi",&teop_em_emcal_phi,"_em_emcal_phi/F");
  outputtree->Branch("_ep_emcal_x",&teop_ep_emcal_x,"_ep_emcal_x/F");
  outputtree->Branch("_em_emcal_x",&teop_em_emcal_x,"_em_emcal_x/F");
  outputtree->Branch("_ep_emcal_y",&teop_ep_emcal_y,"_ep_emcal_y/F");
  outputtree->Branch("_em_emcal_y",&teop_em_emcal_y,"_em_emcal_y/F");
  outputtree->Branch("_ep_emcal_z",&teop_ep_emcal_z,"_ep_emcal_z/F");
  outputtree->Branch("_em_emcal_z",&teop_em_emcal_z,"_em_emcal_z/F");
  outputtree->Branch("_ep_emcal_eta",&teop_ep_emcal_eta,"_ep_emcal_eta/F");
  outputtree->Branch("_em_emcal_eta",&teop_em_emcal_eta,"_em_emcal_eta/F");
  outputtree->Branch("_ep_p",&teop_ep_p,"_ep_p/F");
  outputtree->Branch("_em_p",&teop_em_p,"_em_p/F");
  outputtree->Branch("_ep_p_raw",&teop_ep_p_raw,"_ep_p_raw/F");
  outputtree->Branch("_em_p_raw",&teop_em_p_raw,"_em_p_raw/F");
  outputtree->Branch("_ep_p_unmoved",&teop_ep_p_unmoved,"_ep_p_unmoved/F");
  outputtree->Branch("_em_p_unmoved",&teop_em_p_unmoved,"_em_p_unmoved/F");
  outputtree->Branch("_ep_pt",&teop_ep_pt,"_ep_pt/F");
  outputtree->Branch("_em_pt",&teop_em_pt,"_em_pt/F");
  outputtree->Branch("_ep_pt_raw",&teop_ep_pt_raw,"_ep_pt_raw/F");
  outputtree->Branch("_em_pt_raw",&teop_em_pt_raw,"_em_pt_raw/F");
  outputtree->Branch("_ep_pt_unmoved",&teop_ep_pt_unmoved,"_ep_pt_unmoved/F");
  outputtree->Branch("_em_pt_unmoved",&teop_em_pt_unmoved,"_em_pt_unmoved/F");
  outputtree->Branch("_ep_eta",&teop_ep_eta,"_ep_eta/F");
  outputtree->Branch("_em_eta",&teop_em_eta,"_em_eta/F");
  outputtree->Branch("_ep_phi",&teop_ep_phi,"_ep_phi/F");
  outputtree->Branch("_em_phi",&teop_em_phi,"_em_phi/F");
  outputtree->Branch("_ep_pE",&teop_ep_pE,"_ep_pE/F");
  outputtree->Branch("_em_pE",&teop_em_pE,"_em_pE/F");
  outputtree->Branch("_ep_pE_unmoved",&teop_ep_pE_unmoved,"_ep_pE_unmoved/F");
  outputtree->Branch("_em_pE_unmoved",&teop_em_pE_unmoved,"_em_pE_unmoved/F");
  outputtree->Branch("_ep_mass_1",&teop_ep_mass_1,"_ep_mass_1/F");
  outputtree->Branch("_em_mass_1",&teop_em_mass_1,"_em_mass_1/F");
  outputtree->Branch("_ep_mass_2",&teop_ep_mass_2,"_ep_mass_2/F");
  outputtree->Branch("_em_mass_2",&teop_em_mass_2,"_em_mass_2/F");
  outputtree->Branch("_ep_mass_3",&teop_ep_mass_3,"_ep_mass_3/F");
  outputtree->Branch("_em_mass_3",&teop_em_mass_3,"_em_mass_3/F");
  outputtree->Branch("_ep_phi_projemc",&teop_ep_phi_projemc,"_ep_phi_projemc/F");
  outputtree->Branch("_em_phi_projemc",&teop_em_phi_projemc,"_em_phi_projemc/F");
  outputtree->Branch("_ep_x_projemc",&teop_ep_x_projemc,"_ep_x_projemc/F");
  outputtree->Branch("_em_x_projemc",&teop_em_x_projemc,"_em_x_projemc/F");
  outputtree->Branch("_ep_y_projemc",&teop_ep_y_projemc,"_ep_y_projemc/F");
  outputtree->Branch("_em_y_projemc",&teop_em_y_projemc,"_em_y_projemc/F");
  outputtree->Branch("_ep_z_projemc",&teop_ep_z_projemc,"_ep_z_projemc/F");
  outputtree->Branch("_em_z_projemc",&teop_em_z_projemc,"_em_z_projemc/F");
  outputtree->Branch("_ep_quality",&teop_ep_quality,"_ep_quality/F");
  outputtree->Branch("_em_quality",&teop_em_quality,"_em_quality/F");
  outputtree->Branch("_ep_charge",&teop_ep_charge,"_ep_charge/I");
  outputtree->Branch("_em_charge",&teop_em_charge,"_em_charge/I");
  outputtree->Branch("_ep_crossing",&teop_ep_crossing,"_ep_crossing/I");
  outputtree->Branch("_em_crossing",&teop_em_crossing,"_em_crossing/I");
  outputtree->Branch("_track12_deta",&teop_track12_deta,"_track12_deta/F");
  outputtree->Branch("_track12_dphi",&teop_track12_dphi,"_track12_dphi/F");
  outputtree->Branch("_track12_dz_emc",&teop_track12_dz_emc,"_track12_dz_emc/F");
  outputtree->Branch("_track12_dp",&teop_track12_dp,"_track12_dp/F");
  outputtree->Branch("_epemc_dphi",&teop_epemc_dphi,"_epemc_dphi/F");
  outputtree->Branch("_ememc_dphi",&teop_ememc_dphi,"_ememc_dphi/F");
  outputtree->Branch("_epemc_dz",&teop_epemc_dz,"_epemc_dz/F");
  outputtree->Branch("_ememc_dz",&teop_ememc_dz,"_ememc_dz/F");
  outputtree->Branch("_track12_dca_2d",&teop_track12_dca_2d,"_track12_dca_2d/F");
  outputtree->Branch("_track12_dca_3d",&teop_track12_dca_3d,"_track12_dca_3d/F");
  outputtree->Branch("_ep_clus_x",&teop_ep_clus_x);
  outputtree->Branch("_em_clus_x",&teop_em_clus_x);
  outputtree->Branch("_ep_clus_y",&teop_ep_clus_y);
  outputtree->Branch("_em_clus_y",&teop_em_clus_y);
  outputtree->Branch("_ep_clus_z",&teop_ep_clus_z);
  outputtree->Branch("_em_clus_z",&teop_em_clus_z);
  outputtree->Branch("_eop_ep",&teop_eop_ep,"_eop_ep/F");
  outputtree->Branch("_eop_em",&teop_eop_em,"_eop_em/F");

  outputtree->Branch("_true_gamma_mass",&teop_true_gamma_mass);
  outputtree->Branch("_true_gamma_pE",&teop_true_gamma_pE);
  outputtree->Branch("_true_gamma_eta",&teop_true_gamma_eta);
  outputtree->Branch("_true_gamma_prod_radius",&teop_true_gamma_prod_radius);
  outputtree->Branch("_true_gamma_prod_x",&teop_true_gamma_prod_x);
  outputtree->Branch("_true_gamma_prod_y",&teop_true_gamma_prod_y);
  outputtree->Branch("_true_gamma_prod_z",&teop_true_gamma_prod_z);
  outputtree->Branch("_true_gamma_decay_radius",&teop_true_gamma_decay_radius);
  outputtree->Branch("_true_gamma_decay_x",&teop_true_gamma_decay_x);
  outputtree->Branch("_true_gamma_decay_y",&teop_true_gamma_decay_y);
  outputtree->Branch("_true_gamma_decay_z",&teop_true_gamma_decay_z);
  outputtree->Branch("_true_gamma_mother_id",&teop_true_gamma_mother_id);
  outputtree->Branch("_true_ep_pt",&teop_true_ep_pt);
  outputtree->Branch("_true_em_pt",&teop_true_em_pt);
  outputtree->Branch("_true_ep_p",&teop_true_ep_p);
  outputtree->Branch("_true_em_p",&teop_true_em_p);
  outputtree->Branch("_true_ep_eta",&teop_true_ep_eta);
  outputtree->Branch("_true_em_eta",&teop_true_em_eta);
  outputtree->Branch("_true_ep_phi",&teop_true_ep_phi);
  outputtree->Branch("_true_em_phi",&teop_true_em_phi);
  outputtree->Branch("_true_ep_pE",&teop_true_ep_pE);
  outputtree->Branch("_true_em_pE",&teop_true_em_pE);
  outputtree->Branch("_ep_phi_recotruthDiff",&teop_ep_phi_recotruthDiff);
  outputtree->Branch("_em_phi_recotruthDiff",&teop_em_phi_recotruthDiff);
  outputtree->Branch("_ep_eta_recotruthDiff",&teop_ep_eta_recotruthDiff);
  outputtree->Branch("_em_eta_recotruthDiff",&teop_em_eta_recotruthDiff);
  outputtree->Branch("_SV_phi_recotruthDiff",&teop_SV_phi_recotruthDiff);
  outputtree->Branch("_SV_z_recotruthDiff",&teop_SV_z_recotruthDiff);

  std::vector<int> vec_ep_emcal_matched_index;
  std::vector<float> vec_ep_emcal_matched_e;
  std::vector<float> vec_ep_emcal_matched_phi;
  std::vector<float> vec_ep_emcal_matched_x;
  std::vector<float> vec_ep_emcal_matched_y;
  std::vector<float> vec_ep_emcal_matched_z;
  std::vector<float> vec_ep_emcal_matched_eta;
  std::vector<float> vec_ep_emcal_residual_phi;
  std::vector<float> vec_ep_emcal_residual_z;

  std::vector<int> vec_em_emcal_matched_index;
  std::vector<float> vec_em_emcal_matched_e;
  std::vector<float> vec_em_emcal_matched_phi;
  std::vector<float> vec_em_emcal_matched_x;
  std::vector<float> vec_em_emcal_matched_y;
  std::vector<float> vec_em_emcal_matched_z;
  std::vector<float> vec_em_emcal_matched_eta;
  std::vector<float> vec_em_emcal_residual_phi;
  std::vector<float> vec_em_emcal_residual_z;

  float ep_max_e = 0;
  float ep_min_dphi = 100;
  float ep_min_distance = 100;
  int ep_index = -1;

  float em_max_e = 0;
  float em_min_dphi = 100;
  float em_min_distance = 100;
  int em_index = -1;

  int nevent  = chain->GetEntries();
  cout<<"total nevent = "<<nevent<<endl;
  for(int i = 0; i < nevent; i++)
  {
    if (i % (nevent / 10) == 0)
    {
      cout << "Processing progress: " << i / (nevent / 10) << "0%" << endl;
    }
    chain->GetEntry(i);

    // get truth info
    teop_true_gamma_mass.clear();
    teop_true_gamma_pE.clear();
    teop_true_gamma_eta.clear();
    teop_true_gamma_prod_radius.clear();
    teop_true_gamma_prod_x.clear();
    teop_true_gamma_prod_y.clear();
    teop_true_gamma_prod_z.clear();
    teop_true_gamma_decay_radius.clear();
    teop_true_gamma_decay_x.clear();
    teop_true_gamma_decay_y.clear();
    teop_true_gamma_decay_z.clear();
    teop_true_gamma_mother_id.clear();
    teop_true_ep_pt.clear();
    teop_true_em_pt.clear();
    teop_true_ep_p.clear();
    teop_true_em_p.clear();
    teop_true_ep_eta.clear();
    teop_true_em_eta.clear();
    teop_true_ep_phi.clear();
    teop_true_em_phi.clear();
    teop_true_ep_pE.clear();
    teop_true_em_pE.clear();
    for (int icant = 0; icant < _true_numCan; icant++)
    {
      float m2 = pow(_true_gamma_pE->at(icant),2) - pow(_true_gamma_px->at(icant),2) - pow(_true_gamma_py->at(icant),2) - pow(_true_gamma_pz->at(icant),2);
      float mass = m2 > 0 ? sqrt(m2) : sqrt(-m2);
      float gamma_radius = sqrt( pow(_true_gamma_x->at(icant),2) + pow(_true_gamma_y->at(icant),2) );
      float gamma_eta = cal_eta(_true_gamma_px->at(icant),_true_gamma_py->at(icant),_true_gamma_pz->at(icant));
      float epem_radius = sqrt( pow(_true_ep_x->at(icant),2) + pow(_true_ep_y->at(icant),2) );
      float ep_pt = sqrt( pow(_true_ep_px->at(icant),2) + pow(_true_ep_py->at(icant),2) );
      float em_pt = sqrt( pow(_true_em_px->at(icant),2) + pow(_true_em_py->at(icant),2) );
      float ep_p = sqrt( pow(_true_ep_px->at(icant),2) + pow(_true_ep_py->at(icant),2) + pow(_true_ep_pz->at(icant),2) );
      float em_p = sqrt( pow(_true_em_px->at(icant),2) + pow(_true_em_py->at(icant),2) + pow(_true_em_pz->at(icant),2) );

      if (_true_gamma_embedding_id->at(icant)!=1) {continue;}
      if (!isInRange(-1.1,gamma_eta,1.1)) { continue;}
      if (!isInRange(0,gamma_radius,1)) {continue;}
      if (!isInRange(-20,_true_gamma_z->at(icant),20)) {continue;}
      if (!isInRange(0,epem_radius,25)) {continue;}

      teop_true_gamma_mass.push_back(mass);
      teop_true_gamma_pE.push_back(_true_gamma_pE->at(icant));
      teop_true_gamma_eta.push_back( gamma_eta );
      teop_true_gamma_prod_radius.push_back(gamma_radius);
      teop_true_gamma_prod_x.push_back(_true_gamma_x->at(icant));
      teop_true_gamma_prod_y.push_back(_true_gamma_y->at(icant));
      teop_true_gamma_prod_z.push_back(_true_gamma_z->at(icant));
      teop_true_gamma_mother_id.push_back(_true_gamma_mother_id->at(icant));

      teop_true_gamma_decay_radius.push_back(epem_radius);
      teop_true_gamma_decay_x.push_back(_true_ep_x->at(icant));
      teop_true_gamma_decay_y.push_back(_true_ep_y->at(icant));
      teop_true_gamma_decay_z.push_back(_true_ep_z->at(icant));

      teop_true_ep_phi.push_back(_true_ep_phi->at(icant));
      teop_true_em_phi.push_back(_true_em_phi->at(icant));
      teop_true_ep_eta.push_back(_true_ep_eta->at(icant));
      teop_true_em_eta.push_back(_true_em_eta->at(icant));
      teop_true_ep_pE.push_back(_true_ep_pE->at(icant));
      teop_true_em_pE.push_back(_true_em_pE->at(icant));

      teop_true_ep_pt.push_back( ep_pt );
      teop_true_em_pt.push_back( em_pt );
      teop_true_ep_p.push_back( ep_p );
      teop_true_em_p.push_back( em_p );
    }
    if (teop_true_em_p.size()==0) continue;

    for (int ican = 0; ican < _numCan; ican++)
    {
      if (_em_pT->at(ican) < min_Track_Pt || _ep_pT->at(ican) < min_Track_Pt) continue;

      float quality_cut = 1000;
      if ( (_gamma_chi2->at(ican) / _gamma_nDoF->at(ican))>quality_cut || (_ep_chi2->at(ican) / _ep_nDoF->at(ican))>quality_cut || (_em_chi2->at(ican) / _em_nDoF->at(ican))>quality_cut)
      {
        if (verbosity>0)
        {
          std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail gamma/e+/e- quality cut"<<std::endl;
          std::cout<<"gamma quality = "<<_gamma_chi2->at(ican) / _gamma_nDoF->at(ican)<<" e+ quality = "<<_ep_chi2->at(ican) / _ep_nDoF->at(ican)<<" e- quality = "<<_em_chi2->at(ican) / _em_nDoF->at(ican)<<std::endl;
        }
        //continue;
      }

      float gamma_r_cut = 0;
      float gamma_r = customsqrt(_gamma_x->at(ican)*_gamma_x->at(ican)+_gamma_y->at(ican)*_gamma_y->at(ican));
      if (gamma_r < gamma_r_cut)
      {  
        if (verbosity>0)
        {
          std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail secondary vertex cut"<<std::endl;
        }
        //continue;
      }

      float gamma_mass_cut = 1;
      if (_gamma_mass->at(ican) > gamma_mass_cut)
      {  
        if (verbosity>0)
        {
          std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail mass cut"<<std::endl;
        }
        //continue;
      }

      teop_gamma_mass = _gamma_mass->at(ican);
      teop_gamma_radius = gamma_r;
      teop_gamma_x = _gamma_x->at(ican);
      teop_gamma_y = _gamma_y->at(ican);
      teop_gamma_z = _gamma_z->at(ican);
      teop_gamma_pE = _gamma_pE->at(ican);
      teop_gamma_quality = _gamma_chi2->at(ican) / _gamma_nDoF->at(ican);
      teop_gamma_vertex_volume = _gamma_vertex_volume->at(ican);
      teop_gamma_SV_chi2_per_nDoF = _gamma_SV_chi2_per_nDoF->at(ican);

      vec_ep_emcal_matched_index.clear();
      vec_ep_emcal_matched_e.clear();
      vec_ep_emcal_matched_phi.clear();
      vec_ep_emcal_matched_x.clear();
      vec_ep_emcal_matched_y.clear();
      vec_ep_emcal_matched_z.clear();
      vec_ep_emcal_matched_eta.clear();
      vec_ep_emcal_residual_phi.clear();
      vec_ep_emcal_residual_z.clear();
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
        float emcal_x = radius_scale * _emcal_x->at(iem);
        float emcal_y = radius_scale * _emcal_y->at(iem);
        float emcal_z = radius_scale * _emcal_z->at(iem);
        float emcal_eta = cal_eta(_emcal_x->at(iem), _emcal_y->at(iem), _emcal_z->at(iem));

        float dphi = PiRange(track_phi - emcal_phi);
        float dz = _ep_z_emc->at(ican) - emcal_z;

        //matching
        //if (dphi>-0.2 && dphi<0.2) // no cut in all
        if (isInRange(-0.03,dphi,0.03) && isInRange(-5,dz,5))
        {
          vec_ep_emcal_matched_index.push_back(iem);
          vec_ep_emcal_matched_e.push_back(_emcal_e->at(iem));
          vec_ep_emcal_matched_phi.push_back(emcal_phi);
          vec_ep_emcal_matched_x.push_back(emcal_x);
          vec_ep_emcal_matched_y.push_back(emcal_y);
          vec_ep_emcal_matched_z.push_back(emcal_z);
          vec_ep_emcal_matched_eta.push_back(emcal_eta);
          vec_ep_emcal_residual_phi.push_back(dphi);
          vec_ep_emcal_residual_z.push_back(dz);
        }
      }

      vec_em_emcal_matched_index.clear();
      vec_em_emcal_matched_e.clear();
      vec_em_emcal_matched_phi.clear();
      vec_em_emcal_matched_x.clear();
      vec_em_emcal_matched_y.clear();
      vec_em_emcal_matched_z.clear();
      vec_em_emcal_matched_eta.clear();
      vec_em_emcal_residual_phi.clear();
      vec_em_emcal_residual_z.clear();
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
        float track_phi = cal_phi(_em_x_emc->at(ican), _em_y_emc->at(ican));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_x = radius_scale * _emcal_x->at(iem);
        float emcal_y = radius_scale * _emcal_y->at(iem);
        float emcal_z = radius_scale * _emcal_z->at(iem);
        float emcal_eta = cal_eta(_emcal_x->at(iem), _emcal_y->at(iem), _emcal_z->at(iem));

        float dphi = PiRange(track_phi - emcal_phi);
        float dz = _em_z_emc->at(ican) - emcal_z;

        //matching
        //if (dphi>-0.2 && dphi<0.2) // no cut in all
        if (isInRange(-0.03,dphi,0.03) && isInRange(-5,dz,5))
        {
          vec_em_emcal_matched_index.push_back(iem);
          vec_em_emcal_matched_e.push_back(_emcal_e->at(iem));
          vec_em_emcal_matched_phi.push_back(emcal_phi);
          vec_em_emcal_matched_x.push_back(emcal_x);
          vec_em_emcal_matched_y.push_back(emcal_y);
          vec_em_emcal_matched_z.push_back(emcal_z);
          vec_em_emcal_matched_eta.push_back(emcal_eta);
          vec_em_emcal_residual_phi.push_back(dphi);
          vec_em_emcal_residual_z.push_back(dz);
        }
      }

      if (vec_ep_emcal_matched_e.size() == 0 || vec_em_emcal_matched_e.size() == 0)
      {
        if (verbosity>0)cout<<"no matched???"<<endl;
        continue;
      }

      ep_max_e = 0;
      ep_min_dphi = 100;
      ep_min_distance = 100;
      ep_index = -1;
      for(unsigned int iem = 0; iem < vec_ep_emcal_residual_phi.size(); iem++)
      {
        //if (fabs(vec_ep_emcal_residual_phi.at(iem))<ep_min_dphi)
        //{
        //  ep_min_dphi = fabs(vec_ep_emcal_residual_phi.at(iem));
        //  ep_index = iem;
        //}
        //if (fabs(vec_ep_emcal_matched_e.at(iem))>ep_max_e)
        //{
        //  ep_max_e = vec_ep_emcal_matched_e.at(iem);
        //  ep_index = iem;
        //}
        float R = sqrt(pow(emcal_radius*vec_ep_emcal_residual_phi.at(iem), 2) + pow(vec_ep_emcal_residual_z.at(iem), 2));
        if (R<ep_min_distance)
        {
          ep_min_distance = R;
          ep_index = iem;
        }
      }

      em_max_e = 0;
      em_min_dphi = 100;
      em_min_distance = 100;
      em_index = -1;
      for(unsigned int iem = 0; iem < vec_em_emcal_residual_phi.size(); iem++)
      {
        //if (fabs(vec_em_emcal_residual_phi.at(iem))<em_min_dphi)
        //{
        //  em_min_dphi = fabs(vec_em_emcal_residual_phi.at(iem));
        //  em_index = iem;
        //}
        //if (fabs(vec_em_emcal_matched_e.at(iem))>em_max_e)
        //{
        //  em_max_e = vec_em_emcal_matched_e.at(iem);
        //  em_index = iem;
        //}
        float R = sqrt(pow(emcal_radius*vec_em_emcal_residual_phi.at(iem), 2) + pow(vec_em_emcal_residual_z.at(iem), 2));
        if (R<em_min_distance)
        {
          em_min_distance = R;
          em_index = iem;
        }
      }

      teop_ep_clus_x.clear();
      teop_ep_clus_y.clear();
      teop_ep_clus_z.clear();
      for (int j = 0; j < (_ep_clus_x->size()); j++)
      {
        if (_ep_clus_ican->at(j) != ican) continue;
        teop_ep_clus_x.push_back(_ep_clus_x->at(j));
        teop_ep_clus_y.push_back(_ep_clus_y->at(j));
        teop_ep_clus_z.push_back(_ep_clus_z->at(j));
      }

      teop_em_clus_x.clear();
      teop_em_clus_y.clear();
      teop_em_clus_z.clear();
      for (int j = 0; j < (_em_clus_x->size()); j++)
      {
        if (_em_clus_ican->at(j) != ican) continue;
        teop_em_clus_x.push_back(_em_clus_x->at(j));
        teop_em_clus_y.push_back(_em_clus_y->at(j));
        teop_em_clus_z.push_back(_em_clus_z->at(j));
      }

      teop_ep_emcal_e = vec_ep_emcal_matched_e.at(ep_index);
      teop_ep_emcal_phi = vec_ep_emcal_matched_phi.at(ep_index);
      teop_ep_emcal_x = vec_ep_emcal_matched_x.at(ep_index);
      teop_ep_emcal_y = vec_ep_emcal_matched_y.at(ep_index);
      teop_ep_emcal_z = vec_ep_emcal_matched_z.at(ep_index);
      teop_ep_emcal_eta = vec_ep_emcal_matched_eta.at(ep_index);

      teop_em_emcal_e = vec_em_emcal_matched_e.at(em_index);
      teop_em_emcal_phi = vec_em_emcal_matched_phi.at(em_index);
      teop_em_emcal_x = vec_em_emcal_matched_x.at(em_index);
      teop_em_emcal_y = vec_em_emcal_matched_y.at(em_index);
      teop_em_emcal_z = vec_em_emcal_matched_z.at(em_index);
      teop_em_emcal_eta = vec_em_emcal_matched_eta.at(em_index);

      //if (vec_ep_emcal_matched_index.at(ep_index)==vec_em_emcal_matched_index.at(em_index)) continue;

      teop_ep_p = _ep_p->at(ican);
      teop_ep_p_raw = _ep_p_raw->at(ican);
      teop_ep_p_unmoved = _ep_p_unmoved->at(ican);
      teop_ep_pt = _ep_pT->at(ican);
      teop_ep_pt_raw = _ep_pT_raw->at(ican);
      teop_ep_pt_unmoved = _ep_pT_unmoved->at(ican);
      teop_ep_pE_unmoved = _ep_pE_unmoved->at(ican);
      teop_ep_eta = _ep_pseudorapidity->at(ican);
      teop_ep_phi = _ep_phi->at(ican);
      teop_ep_pE = _ep_pE->at(ican);
      teop_ep_mass_1 = _ep_mass->at(ican);
      teop_ep_mass_2 = customsqrt( pow(_ep_pE->at(ican),2) - pow(_ep_p->at(ican),2) );
      teop_ep_mass_3 = customsqrt( pow(_ep_pE_unmoved->at(ican),2) - pow(_ep_p_unmoved->at(ican),2) );
      teop_ep_phi_projemc = cal_phi(_ep_x_emc->at(ican), _ep_y_emc->at(ican));
      teop_ep_x_projemc = _ep_x_emc->at(ican);
      teop_ep_y_projemc = _ep_y_emc->at(ican);
      teop_ep_z_projemc = _ep_z_emc->at(ican);
      teop_ep_quality = _ep_chi2->at(ican) / _ep_nDoF->at(ican);
      teop_ep_charge = 1;
      teop_ep_crossing = _ep_crossing->at(ican);

      teop_em_p = _em_p->at(ican);
      teop_em_p_raw = _em_p_raw->at(ican);
      teop_em_p_unmoved = _em_p_unmoved->at(ican);
      teop_em_pt = _em_pT->at(ican);
      teop_em_pt_raw = _em_pT_raw->at(ican);
      teop_em_pt_unmoved = _em_pT_unmoved->at(ican);
      teop_em_pE_unmoved = _em_pE_unmoved->at(ican);
      teop_em_eta = _em_pseudorapidity->at(ican);
      teop_em_phi = _em_phi->at(ican);
      teop_em_pE = _em_pE->at(ican);
      teop_em_mass_1 = _em_mass->at(ican);
      teop_em_mass_2 = customsqrt( pow(_em_pE->at(ican),2) - pow(_em_p->at(ican),2) );
      teop_em_mass_3 = customsqrt( pow(_em_pE_unmoved->at(ican),2) - pow(_em_p_unmoved->at(ican),2) );
      teop_em_phi_projemc = cal_phi(_em_x_emc->at(ican), _em_y_emc->at(ican));
      teop_em_x_projemc = _em_x_emc->at(ican);
      teop_em_y_projemc = _em_y_emc->at(ican);
      teop_em_z_projemc = _em_z_emc->at(ican);
      teop_em_quality = _em_chi2->at(ican) / _em_nDoF->at(ican);
      teop_em_charge = -1;
      teop_em_crossing = _em_crossing->at(ican);

      teop_track12_deta = _ep_pseudorapidity->at(ican) - _em_pseudorapidity->at(ican);
      teop_track12_dphi = _ep_phi->at(ican) - _em_phi->at(ican);
      teop_track12_dz_emc = _ep_z_emc->at(ican) - _em_z_emc->at(ican);
      teop_track12_dp = _ep_p->at(ican) - _em_p->at(ican);
      teop_epemc_dphi = vec_ep_emcal_residual_phi.at(ep_index);
      teop_ememc_dphi = vec_em_emcal_residual_phi.at(em_index);
      teop_epemc_dz = vec_ep_emcal_residual_z.at(ep_index);
      teop_ememc_dz = vec_em_emcal_residual_z.at(em_index);
      teop_track12_dca_2d = _epem_DCA_2d->at(ican);
      teop_track12_dca_3d = _epem_DCA_3d->at(ican);

      float eop_ep = teop_ep_emcal_e / teop_ep_p;
      float eop_em = teop_em_emcal_e / teop_em_p;

      if (eop_ep > 0.8 && eop_ep < 1.3)
      {
        teop_eop_em = eop_em;
      }
      else
      {
        teop_eop_em = -1;
      }

      if (eop_em > 0.8 && eop_em < 1.3)
      {
        teop_eop_ep = eop_ep;
      }
      else
      {
        teop_eop_ep = -1;
      }

      teop_ep_phi_recotruthDiff.clear();
      teop_em_phi_recotruthDiff.clear();
      teop_ep_eta_recotruthDiff.clear();
      teop_em_eta_recotruthDiff.clear();
      teop_SV_phi_recotruthDiff.clear();
      teop_SV_z_recotruthDiff.clear();
      for (int icant = 0; icant < (teop_true_gamma_decay_radius.size()); icant++)
      {
        teop_ep_phi_recotruthDiff.push_back(PiRange(teop_ep_phi - teop_true_ep_phi.at(icant)));
        teop_em_phi_recotruthDiff.push_back(PiRange(teop_em_phi - teop_true_em_phi.at(icant)));
        teop_ep_eta_recotruthDiff.push_back(teop_ep_eta - teop_true_ep_eta.at(icant));
        teop_em_eta_recotruthDiff.push_back(teop_em_eta - teop_true_em_eta.at(icant));
        teop_SV_phi_recotruthDiff.push_back(PiRange(atan2(teop_gamma_y,teop_gamma_x) - atan2(teop_true_gamma_decay_y.at(icant),teop_true_gamma_decay_x.at(icant))));
        teop_SV_z_recotruthDiff.push_back(teop_gamma_z - teop_true_gamma_decay_z.at(icant));
      }

      outputtree->Fill();
    }
  }

  outputfile->cd();
  outputfile->Write();
  outputfile->Close();
}
