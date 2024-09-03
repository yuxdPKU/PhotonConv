#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void EoP_kfp(int runnumber)
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  int runlist[77] = {51452,51490,51491,51493,51505,51508,51509,51510,51511,51512,51518,51520,51559,51560,51562,51563,51564,51565,51569,51570,51571,51572,51575,51576,51603,51604,51605,51606,51607,51608,51609,51610,51611,51613,51614,51615,51616,51619,51688,51710,51712,51713,51714,51715,51718,51719,51721,51725,51732,51733,51735,51736,51738,51740,51742,51762,51764,51765,51768,51771,51777,51778,51824,51831,51838,51840,51841,51842,51843,51854,51855,51856,51860,51865,51874,51878,51886};

  TChain* chain = new TChain("tree_KFP");
  //chain->Add(Form("./%d/TrackCalo_*_ana.root",runnumber));
  for (int i=0; i<77; i++) chain->Add(Form("./%d/TrackCalo_*_ana.root",runlist[i]));

  setBranch_kfp(chain);

  TFile* outputfile = new TFile(Form("./eop_%d_kfp.root",runnumber),"recreate");
  TTree* outputtree = new TTree("tree","tree with eop info");

  float teop_gamma_mass, teop_gamma_radius;
  float teop_gamma_quality, teop_track_quality;
  float teop_emcal_e, teop_emcal_phi, teop_emcal_z;
  float teop_track_p, teop_track_phi, teop_track_z;
  float teop_dphi, teop_dz;
  int teop_track_charge;
  float teop_dca;

  outputtree->Branch("_runNumber",&_runNumber,"_runNumber/I");
  outputtree->Branch("_eventNumber",&_eventNumber,"_eventNumber/I");
  outputtree->Branch("_gamma_mass",&teop_gamma_mass,"_gamma_mass/F");
  outputtree->Branch("_gamma_radius",&teop_gamma_radius,"_gamma_radius/F");
  outputtree->Branch("_gamma_quality",&teop_gamma_quality,"_gamma_quality/F");
  outputtree->Branch("_emcal_e",&teop_emcal_e,"_emcal_e/F");
  outputtree->Branch("_emcal_phi",&teop_emcal_phi,"_emcal_phi/F");
  outputtree->Branch("_emcal_z",&teop_emcal_z,"_emcal_z/F");
  outputtree->Branch("_track_p",&teop_track_p,"_track_p/F");
  outputtree->Branch("_track_phi",&teop_track_phi,"_track_phi/F");
  outputtree->Branch("_track_z",&teop_track_z,"_track_z/F");
  outputtree->Branch("_track_quality",&teop_track_quality,"_track_quality/F");
  outputtree->Branch("_track_charge",&teop_track_charge,"_track_charge/I");
  outputtree->Branch("_dphi",&teop_dphi,"_dphi/F");
  outputtree->Branch("_dz",&teop_dz,"_dz/F");
  outputtree->Branch("_dca",&teop_dca,"_dca/F");

  std::vector<float> vec_emcal_matched_e;
  std::vector<float> vec_emcal_matched_phi;
  std::vector<float> vec_emcal_matched_z;
  std::vector<float> vec_track_emcal_residual_phi;
  std::vector<float> vec_track_emcal_residual_z;

  float max_e = 0;
  int index_maxE = -1;

  int nevent  = chain->GetEntries();
  cout<<"total nevent = "<<nevent<<endl;
  for(int i = 0; i < nevent; i++)
  {
    if (i % (nevent / 10) == 0) cout << "Processing progress: " << i / (nevent / 10) << "0%" << endl;
    chain->GetEntry(i);

    for (int ican = 0; ican < _numCan; ican++)
    {
      float quality_cut = 1000;
      if ( (_gamma_chi2->at(ican) / _gamma_nDoF->at(ican))>quality_cut || (_ep_chi2->at(ican) / _ep_nDoF->at(ican))>quality_cut || (_em_chi2->at(ican) / _em_nDoF->at(ican))>quality_cut)
      {
        //std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail gamma/e+/e- quality cut"<<std::endl;
        //std::cout<<"gamma quality = "<<_gamma_chi2->at(ican) / _gamma_nDoF->at(ican)<<" e+ quality = "<<_ep_chi2->at(ican) / _ep_nDoF->at(ican)<<" e- quality = "<<_em_chi2->at(ican) / _em_nDoF->at(ican)<<std::endl;
        //continue;
      }

      float gamma_r_cut = 0;
      float gamma_r = sqrt(_gamma_x->at(ican)*_gamma_x->at(ican)+_gamma_y->at(ican)*_gamma_y->at(ican));
      if (gamma_r < gamma_r_cut)
      {  
        //std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail secondary vertex cut"<<std::endl;
        //continue;
      }

      float gamma_mass_cut = 1;
      if (_gamma_mass->at(ican) > gamma_mass_cut)
      {  
        //std::cout<<"Run "<<_runNumber<<" Event "<<_eventNumber<<" Candidate "<<ican<<": fail mass cut"<<std::endl;
        //continue;
      }

      teop_gamma_mass = _gamma_mass->at(ican);
      teop_gamma_radius = gamma_r;
      teop_gamma_quality = _gamma_chi2->at(ican) / _gamma_nDoF->at(ican);
//cout<<"track e+ px = "<<_ep_px->at(ican)<<" py = "<<_ep_py->at(ican)<<" pz = "<<_ep_pz->at(ican)<<" phi_emc = "<<_ep_phi_emc->at(ican)<<" z_emc = "<<_ep_z_emc->at(ican)<<endl;
//cout<<"track e- px = "<<_em_px->at(ican)<<" py = "<<_em_py->at(ican)<<" pz = "<<_em_pz->at(ican)<<" phi_emc = "<<_em_phi_emc->at(ican)<<" z_emc = "<<_em_z_emc->at(ican)<<endl;

      // match e+ with calo
      vec_emcal_matched_e.clear();
      vec_emcal_matched_phi.clear();
      vec_emcal_matched_z.clear();
      vec_track_emcal_residual_phi.clear();
      vec_track_emcal_residual_z.clear();
      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
//cout<<"emcal "<<iem<<" e = "<<_emcal_e->at(iem)<<" x = "<<_emcal_x->at(iem)<<" y = "<<_emcal_y->at(iem)<<" z = "<<_emcal_z->at(iem)<<endl;
        // cut good emcal clusters
        if(_emcal_e->at(iem) < min_EMCal_E) continue;
        //if(fabs(_emcal_eta->at(iem)) > 1.1) continue;
        float emcal_phi = atan2(_emcal_y->at(iem), _emcal_x->at(iem));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_z = radius_scale * _emcal_z->at(iem);

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
          vec_track_emcal_residual_phi.push_back(dphi);
          vec_track_emcal_residual_z.push_back(dz);
        }

      }

      max_e = 0;
      index_maxE = -1;
      for(unsigned int iem = 0; iem < vec_emcal_matched_e.size(); iem++)
      {
        if (vec_emcal_matched_e.at(iem)>max_e)
        {
          max_e = vec_emcal_matched_e.at(iem);
          index_maxE = iem;
        }
      }

      if (vec_emcal_matched_e.size() == 0)
      {
cout<<"no matched???"<<endl;
        continue;
      }

      teop_emcal_e = max_e;
      teop_emcal_phi = vec_emcal_matched_phi.at(index_maxE);
      teop_emcal_z = vec_emcal_matched_z.at(index_maxE);
      teop_track_p = _ep_p->at(ican);
      teop_track_phi = _ep_phi_emc->at(ican);
      teop_track_z = _ep_z_emc->at(ican);
      teop_track_quality = _ep_chi2->at(ican) / _ep_nDoF->at(ican);
      teop_track_charge = 1;
      teop_dphi = vec_track_emcal_residual_phi.at(index_maxE);
      teop_dz = vec_track_emcal_residual_z.at(index_maxE);
      teop_dca = _epem_DCA->at(ican);

      outputtree->Fill();
//cout<<"fill e+"<<endl;

      // match e- with calo
      vec_emcal_matched_e.clear();
      vec_emcal_matched_phi.clear();
      vec_emcal_matched_z.clear();
      vec_track_emcal_residual_phi.clear();
      vec_track_emcal_residual_z.clear();
      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
//cout<<"emcal "<<iem<<" e = "<<_emcal_e->at(iem)<<" x = "<<_emcal_x->at(iem)<<" y = "<<_emcal_y->at(iem)<<" z = "<<_emcal_z->at(iem)<<endl;
        // cut good emcal clusters
        if(_emcal_e->at(iem) < min_EMCal_E) continue;
        //if(fabs(_emcal_eta->at(iem)) > 1.1) continue;
        float emcal_phi = atan2(_emcal_y->at(iem), _emcal_x->at(iem));
        float radius_scale = emcal_radius / sqrt( pow(_emcal_x->at(iem),2) + pow(_emcal_y->at(iem),2) );
        float emcal_z = radius_scale * _emcal_z->at(iem);

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
          vec_track_emcal_residual_phi.push_back(dphi);
          vec_track_emcal_residual_z.push_back(dz);
        }

      }

      max_e = 0;
      index_maxE = -1;
      for(unsigned int iem = 0; iem < vec_emcal_matched_e.size(); iem++)
      {
        if (vec_emcal_matched_e.at(iem)>max_e)
        {
          max_e = vec_emcal_matched_e.at(iem);
          index_maxE = iem;
        }
      }

      if (vec_emcal_matched_e.size() == 0)
      {
cout<<"no matched???"<<endl;
        continue;
      }

      teop_emcal_e = max_e;
      teop_emcal_phi = vec_emcal_matched_phi.at(index_maxE);
      teop_emcal_z = vec_emcal_matched_z.at(index_maxE);
      teop_track_p = _em_p->at(ican);
      teop_track_phi = _em_phi_emc->at(ican);
      teop_track_z = _em_z_emc->at(ican);
      teop_track_quality = _em_chi2->at(ican) / _em_nDoF->at(ican);
      teop_track_charge = -1;
      teop_dphi = vec_track_emcal_residual_phi.at(index_maxE);
      teop_dz = vec_track_emcal_residual_z.at(index_maxE);
      teop_dca = _epem_DCA->at(ican);

      outputtree->Fill();
//cout<<"fill e-"<<endl;
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
