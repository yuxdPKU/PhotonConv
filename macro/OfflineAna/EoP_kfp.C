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
  float teop_track_pt, teop_track_eta;
  int teop_track_charge, teop_track_crossing;
  float teop_track_mass_1, teop_track_mass_2;
  float teop_track12_deta, teop_track12_dphi;
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
  outputtree->Branch("_track_pt",&teop_track_pt,"_track_pt/F");
  outputtree->Branch("_track_eta",&teop_track_eta,"_track_eta/F");
  outputtree->Branch("_track_mass_1",&teop_track_mass_1,"_track_mass_1/F");
  outputtree->Branch("_track_mass_2",&teop_track_mass_2,"_track_mass_2/F");
  outputtree->Branch("_track_phi_projemc",&teop_track_phi_projemc,"_track_phi_projemc/F");
  outputtree->Branch("_track_z_projemc",&teop_track_z_projemc,"_track_z_projemc/F");
  outputtree->Branch("_track_quality",&teop_track_quality,"_track_quality/F");
  outputtree->Branch("_track_charge",&teop_track_charge,"_track_charge/I");
  outputtree->Branch("_track_crossing",&teop_track_crossing,"_track_crossing/I");
  outputtree->Branch("_track12_deta",&teop_track12_deta,"_track12_deta/F");
  outputtree->Branch("_track12_dphi",&teop_track12_dphi,"_track12_dphi/F");
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
        //if (verbosity>0) std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail gamma/e+/e- quality cut"<<std::endl;
        //if (verbosity>0) std::cout<<"gamma quality = "<<_gamma_chi2->at(ican) / _gamma_nDoF->at(ican)<<" e+ quality = "<<_ep_chi2->at(ican) / _ep_nDoF->at(ican)<<" e- quality = "<<_em_chi2->at(ican) / _em_nDoF->at(ican)<<std::endl;
        //continue;
      }

      float gamma_r_cut = 0;
      float gamma_r = sqrt(_gamma_x->at(ican)*_gamma_x->at(ican)+_gamma_y->at(ican)*_gamma_y->at(ican));
      if (gamma_r < gamma_r_cut)
      {  
        //if (verbosity>0) std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail secondary vertex cut"<<std::endl;
        //continue;
      }

      float gamma_mass_cut = 1;
      if (_gamma_mass->at(ican) > gamma_mass_cut)
      {  
        //if (verbosity>0) std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail mass cut"<<std::endl;
        //continue;
      }

      teop_gamma_mass = _gamma_mass->at(ican);
      teop_gamma_radius = gamma_r;
      teop_gamma_quality = _gamma_chi2->at(ican) / _gamma_nDoF->at(ican);
      teop_gamma_vertex_volume = _gamma_vertex_volume->at(ican);
      teop_gamma_SV_chi2_per_nDoF = _gamma_SV_chi2_per_nDoF->at(ican);
      if (verbosity>0)
      {
        cout<<"track e+ px = "<<_ep_px->at(ican)<<" py = "<<_ep_py->at(ican)<<" pz = "<<_ep_pz->at(ican)<<" phi_emc = "<<_ep_phi_emc->at(ican)<<" z_emc = "<<_ep_z_emc->at(ican)<<endl;
        cout<<"track e- px = "<<_em_px->at(ican)<<" py = "<<_em_py->at(ican)<<" pz = "<<_em_pz->at(ican)<<" phi_emc = "<<_em_phi_emc->at(ican)<<" z_emc = "<<_em_z_emc->at(ican)<<endl;
      }

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
        float emcal_phi = atan2(_emcal_y->at(iem), _emcal_x->at(iem));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_z = radius_scale * _emcal_z->at(iem);
        float emcal_eta = asinh(_emcal_z->at(iem)/sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) ));

        float dphi = PiRange(_ep_phi_emc->at(ican) - emcal_phi);
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

      if (vec_emcal_matched_e.size() == 0)
      {
        if (verbosity>0)cout<<"no matched???"<<endl;
        continue;
      }

      teop_emcal_e = vec_emcal_matched_e.at(index);
      teop_emcal_phi = vec_emcal_matched_phi.at(index);
      teop_emcal_z = vec_emcal_matched_z.at(index);
      teop_emcal_eta = vec_emcal_matched_eta.at(index);
      teop_track_p = _ep_p->at(ican);
      teop_track_pt = _ep_pT->at(ican);
      teop_track_eta = _ep_pseudorapidity->at(ican);
      teop_track_mass_1 = _ep_mass->at(ican);
      teop_track_mass_2 = sqrt( pow(_ep_pE->at(ican),2) - pow(_ep_px->at(ican),2) - pow(_ep_py->at(ican),2) - pow(_ep_pz->at(ican),2));
      teop_track_phi_projemc = _ep_phi_emc->at(ican);
      teop_track_z_projemc = _ep_z_emc->at(ican);
      teop_track_quality = _ep_chi2->at(ican) / _ep_nDoF->at(ican);
      teop_track_charge = 1;
      teop_track_crossing = _ep_crossing->at(ican);
      teop_track12_deta = _ep_pseudorapidity->at(ican) - _em_pseudorapidity->at(ican);
      teop_track12_dphi = _ep_phi->at(ican) - _em_phi->at(ican);
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
        float emcal_phi = atan2(_emcal_y->at(iem), _emcal_x->at(iem));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_z = radius_scale * _emcal_z->at(iem);
        float emcal_eta = asinh(_emcal_z->at(iem)/sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) ));

        float dphi = PiRange(_em_phi_emc->at(ican) - emcal_phi);
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

      if (vec_emcal_matched_e.size() == 0)
      {
        if (verbosity>0)cout<<"no matched???"<<endl;
        continue;
      }

      teop_emcal_e = vec_emcal_matched_e.at(index);
      teop_emcal_phi = vec_emcal_matched_phi.at(index);
      teop_emcal_z = vec_emcal_matched_z.at(index);
      teop_emcal_eta = vec_emcal_matched_eta.at(index);
      teop_track_p = _em_p->at(ican);
      teop_track_pt = _em_pT->at(ican);
      teop_track_eta = _em_pseudorapidity->at(ican);
      teop_track_mass_1 = _em_mass->at(ican);
      teop_track_mass_2 = sqrt( pow(_em_pE->at(ican),2) - pow(_em_px->at(ican),2) - pow(_em_py->at(ican),2) - pow(_em_pz->at(ican),2));
      teop_track_phi_projemc = _em_phi_emc->at(ican);
      teop_track_z_projemc = _em_z_emc->at(ican);
      teop_track_quality = _em_chi2->at(ican) / _em_nDoF->at(ican);
      teop_track_charge = -1;
      teop_track_crossing = _em_crossing->at(ican);
      teop_track12_deta = _ep_pseudorapidity->at(ican) - _em_pseudorapidity->at(ican);
      teop_track12_dphi = _ep_phi->at(ican) - _em_phi->at(ican);
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

/*
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

  TH1* h1_nmatch = new TH1F("h1_nmatch", "h1_nmatch", 10, 0, 10);
  h1_nmatch->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_nmatch->GetXaxis()->SetTitle("Ncalo");
  h1_nmatch->GetYaxis()->SetTitle(Form("Events / %.1f",(10.-0.)/10.));
  h1_nmatch->SetMinimum(0.9);

  TH2* h2_dphi_dz_track_emcal = new TH2F("h2_dphi_dz_track_emcal", "h2_dphi_dz_track_emcal", 100, -0.2, 0.2, 100, -20, 20);
  h2_dphi_dz_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h2_dphi_dz_track_emcal->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_dz_track_emcal->GetYaxis()->SetTitle("#DeltaZ [cm]");
  h2_dphi_dz_track_emcal->GetZaxis()->SetTitle("Entries");

  TH1* h1_dphi_track_emcal = new TH1F("h1_dphi_track_emcal", "h1_dphi_track_emcal", 100, -.2, .2);
  h1_dphi_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_dphi_track_emcal->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_dphi_track_emcal->GetYaxis()->SetTitle(Form("Events / %.3f rad",0.4/100.));
  h1_dphi_track_emcal->SetMinimum(0);

  TH1* h1_dz_track_emcal = new TH1F("h1_dz_track_emcal", "h1_dz_track_emcal", 100, -20, 20);
  h1_dz_track_emcal->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
  h1_dz_track_emcal->GetXaxis()->SetTitle("#DeltaZ [cm]");
  h1_dz_track_emcal->GetYaxis()->SetTitle(Form("Events / %1f cm",40./100.));
  h1_dz_track_emcal->SetMinimum(0);

  //loop over events
  for (int index = 0; index<(matched_emcal_e.size()); index++)
  {
    for (int i = 0; i < (matched_emcal_e.at(index).size()); i++)
    {
      h2_eOp_p->Fill(matched_track_p.at(index).at(i), matched_emcal_e.at(index).at(i) / matched_track_p.at(index).at(i));
      h1_eOp->Fill(matched_emcal_e.at(index).at(i) / matched_track_p.at(index).at(i));
    }
  }

  for (int index = 0; index<(matched_number.size()); index++)
  {
    for (int i = 0; i < (matched_number.at(index).size()); i++)
    {
      h1_nmatch->Fill(matched_number.at(index).at(i));
    }
  }

  for (int index = 0; index<(matched_track_emcal_dphi.size()); index++)
  {
    for (int i = 0; i < (matched_track_emcal_dphi.at(index).size()); i++)
    {
      h2_dphi_dz_track_emcal->Fill(matched_track_emcal_dphi.at(index).at(i),matched_track_emcal_dz.at(index).at(i));
      if (matched_track_emcal_dphi.at(index).at(i)>-0.2 && matched_track_emcal_dphi.at(index).at(i)<1.0)
      {
        h1_dz_track_emcal->Fill(matched_track_emcal_dz.at(index).at(i));
      }
      if (matched_track_emcal_dz.at(index).at(i)>-5 && matched_track_emcal_dz.at(index).at(i)<5)
      {
        h1_dphi_track_emcal->Fill(matched_track_emcal_dphi.at(index).at(i));
      }
    }
  }

  TPaveText *pt = new TPaveText(.50, .72, .80, .92, "NDC");
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt->AddText("p+p #sqrt{s}=200 GeV");
  pt->AddText(Form("Run %d",runnumber));

  TPaveText *pt1 = new TPaveText(.20, .72, .40, .92, "NDC");
  pt1->SetFillColor(0);
  //pt1->SetFillStyle(0);//transparent
  pt1->SetLineColor(0);
  pt1->SetBorderSize(0);
  pt1->SetTextColor(kBlack);
  pt1->AddText("#it{#bf{sPHENIX}} Internal");
  pt1->AddText("p+p #sqrt{s}=200 GeV");
  pt1->AddText(Form("Run %d",runnumber));

  TCanvas *can1 = new TCanvas(Form("can1"), "can1", 800, 600);
  //can1->SetLeftMargin(0.10);
  can1->SetRightMargin(0.15);
  can1->cd();
  h2_eOp_p->Draw("COLZ");
  pt->Draw("same");

  TCanvas *can2 = new TCanvas(Form("can2"), "can2", 800, 600);
  //can2->SetLeftMargin(0.12);
  //can2->SetRightMargin(0.05);
  can2->cd();
  h1_eOp->Draw("");
  pt->Draw("same");

  TCanvas *can3 = new TCanvas(Form("can3"), "can3", 800, 600);
  //can3->SetLeftMargin(0.12);
  //can3->SetRightMargin(0.05);
  can3->cd();
  can3->SetLogy(1);
  h1_nmatch->Draw("");
  pt->Draw("same");

  TLine *vline_left = new TLine(-0.02, -20, -0.02, 20);
  vline_left->SetLineColor(kRed);
  vline_left->SetLineStyle(2);
  vline_left->SetLineWidth(1);

  TLine *vline_right = new TLine(0.1, -20, 0.1, 20);
  vline_right->SetLineColor(kRed);
  vline_right->SetLineStyle(2);
  vline_right->SetLineWidth(1);

  TLine *hline_up = new TLine(-0.2, 5, 0.2, 5);
  hline_up->SetLineColor(kRed);
  hline_up->SetLineStyle(2);
  hline_up->SetLineWidth(1);

  TLine *hline_bot = new TLine(-0.2, -5, 0.2, -5);
  hline_bot->SetLineColor(kRed);
  hline_bot->SetLineStyle(2);
  hline_bot->SetLineWidth(1);

  TCanvas *can4 = new TCanvas(Form("can4"), "can4", 800, 600);
  //can4->SetLeftMargin(0.10);
  can4->SetRightMargin(0.15);
  can4->cd();
  h2_dphi_dz_track_emcal->Draw("COLZ");
  vline_left->Draw("same");
  vline_right->Draw("same");
  hline_up->Draw("same");
  hline_bot->Draw("same");
  pt1->Draw("same");

  TArrow *arrow5_left = new TArrow(-5,6000,-5,0,0.03,">");
  arrow5_left->SetFillColor(4);
  arrow5_left->SetFillStyle(1001);
  arrow5_left->SetLineColor(4);
  arrow5_left->SetLineWidth(1);

  TArrow *arrow5_right = new TArrow(5,6000,5,0,0.03,">");
  arrow5_right->SetFillColor(4);
  arrow5_right->SetFillStyle(1001);
  arrow5_right->SetLineColor(4);
  arrow5_right->SetLineWidth(1);

  TCanvas *can5 = new TCanvas(Form("can5"), "can5", 800, 600);
  //can5->SetLeftMargin(0.10);
  //can5->SetRightMargin(0.15);
  can5->cd();
  h1_dz_track_emcal->Draw("");
  arrow5_left->Draw("same");
  arrow5_right->Draw("same");
  pt1->Draw("same");

  TArrow *arrow6_left = new TArrow(-0.02,2000,-0.02,0,0.03,">");
  arrow6_left->SetFillColor(4);
  arrow6_left->SetFillStyle(1001);
  arrow6_left->SetLineColor(4);
  arrow6_left->SetLineWidth(1);

  TArrow *arrow6_right = new TArrow(0.1,2000,0.1,0,0.03,">");
  arrow6_right->SetFillColor(4);
  arrow6_right->SetFillStyle(1001);
  arrow6_right->SetLineColor(4);
  arrow6_right->SetLineWidth(1);

  TCanvas *can6 = new TCanvas(Form("can6"), "can6", 800, 600);
  //can6->SetLeftMargin(0.10);
  //can6->SetRightMargin(0.15);
  can6->cd();
  h1_dphi_track_emcal->Draw("");
  arrow6_left->Draw("same");
  arrow6_right->Draw("same");
  pt1->Draw("same");

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

  can1->Update();
  can1->SaveAs(Form("figure/%d/EoP_p_run%d.pdf",runnumber,runnumber));

  can2->Update();
  can2->SaveAs(Form("figure/%d/EoP_run%d.pdf",runnumber,runnumber));

  can3->Update();
  can3->SaveAs(Form("figure/%d/ncalo_matched_run%d.pdf",runnumber,runnumber));

  can4->Update();
  can4->SaveAs(Form("figure/%d/dphi_vs_dz_run%d.pdf",runnumber,runnumber));

  can5->Update();
  can5->SaveAs(Form("figure/%d/dz_signalcut_run%d.pdf",runnumber,runnumber));

  can6->Update();
  can6->SaveAs(Form("figure/%d/dphi_signalcut_run%d.pdf",runnumber,runnumber));
*/

}
