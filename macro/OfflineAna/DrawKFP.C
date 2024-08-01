#include <filesystem>
#include "utilities.h"

namespace fs = std::filesystem;

void DrawKFP(int runnumber)
{
  gStyle->SetOptStat(0);

  TString inputfile = Form("46730/final_%d_kfp.root",runnumber);

  TFile *file = new TFile(inputfile, "READ");
  TTree *tree = (TTree*)file->Get("DecayTree");

  float gamma_mass;
  float gamma_x;
  float gamma_y;
  float gamma_z;
  float track_1_track_2_DCA;
  float track_1_chi2;
  float track_2_chi2;
  int track_1_nDoF;
  int track_2_nDoF;

  tree->SetBranchAddress("gamma_mass", &gamma_mass);
  tree->SetBranchAddress("gamma_x", &gamma_x);
  tree->SetBranchAddress("gamma_y", &gamma_y);
  tree->SetBranchAddress("gamma_z", &gamma_z);
  tree->SetBranchAddress("track_1_track_2_DCA", &track_1_track_2_DCA);
  tree->SetBranchAddress("track_1_chi2", &track_1_chi2);
  tree->SetBranchAddress("track_2_chi2", &track_2_chi2);
  tree->SetBranchAddress("track_1_nDoF", &track_1_nDoF);
  tree->SetBranchAddress("track_2_nDoF", &track_2_nDoF);

  TH1F* h1_mass = new TH1F("h1_mass", "h1_mass", 20, 0, 1);
  h1_mass->SetTitle(Form("Run 46730"));
  h1_mass->GetXaxis()->SetTitle("#it{M}_{e^{+}e^{-}} [GeV/#it{c}^{2}]");
  h1_mass->GetYaxis()->SetTitle("Events");
  h1_mass->SetMinimum(0);

  TH1F* h1_radius = new TH1F("h1_radius", "h1_radius", 100, 0, 100);
  h1_radius->SetTitle(Form("Run 46730"));
  h1_radius->GetXaxis()->SetTitle("Radius [cm]");
  h1_radius->GetYaxis()->SetTitle("Events");
  h1_radius->SetMinimum(0);

  TH1F* h1_DCA = new TH1F("h1_DCA", "h1_DCA", 20, 0, 1);
  h1_DCA->SetTitle(Form("Run 46730"));
  h1_DCA->GetXaxis()->SetTitle("DCA_{e^{+}e^{-}} [cm]");
  h1_DCA->GetYaxis()->SetTitle("Events");
  h1_DCA->SetMinimum(0);

  TH2* h2_mass_radius = new TH2F("h2_mass_radius", "h2_mass_radius", 20, 0, 1, 100, 0, 100);
  h2_mass_radius->SetTitle("Run 46730");
  h2_mass_radius->GetXaxis()->SetTitle("#it{M}_{e^{+}e^{-}} [GeV/#it{c}^{2}]");
  h2_mass_radius->GetYaxis()->SetTitle("Radius [cm]");
  h2_mass_radius->GetYaxis()->SetTitleOffset(1.2);
  h2_mass_radius->GetZaxis()->SetTitle("Entries");
  h2_mass_radius->GetZaxis()->SetTitleOffset(1.2);

  for(int i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    if (track_1_chi2 / track_1_nDoF > 100 || track_2_chi2 / track_2_nDoF > 100) continue;

    h1_mass->Fill(gamma_mass);
    h1_radius->Fill(sqrt(gamma_x*gamma_x+gamma_y*gamma_y));
    h1_DCA->Fill(track_1_track_2_DCA);
    h2_mass_radius->Fill(gamma_mass, sqrt(gamma_x*gamma_x+gamma_y*gamma_y));
  }

  TCanvas *can1 = new TCanvas("can1", "can1", 800, 800);
  can1->SetLeftMargin(0.12);
  can1->SetRightMargin(0.05);
  can1->cd();
  h1_mass->Draw();

  TCanvas *can2 = new TCanvas("can2", "can2", 800, 800);
  can2->SetLeftMargin(0.12);
  can2->SetRightMargin(0.05);
  can2->cd();
  h1_radius->Draw();

  TCanvas *can3 = new TCanvas("can", "can3", 800, 800);
  can3->SetLeftMargin(0.12);
  can3->SetRightMargin(0.05);
  can3->cd();
  h1_DCA->Draw();

  TCanvas *can4 = new TCanvas("can4", "can4", 800, 800);
  can4->SetLeftMargin(0.12);
  can4->SetRightMargin(0.13);
  can4->cd();
  h2_mass_radius->Draw("COLZ");

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
  can1->SaveAs("figure/46730/TrackEMcal_KFP_mass_run46730.pdf");

  can2->Update();
  can2->SaveAs("figure/46730/TrackEMcal_KFP_radius_run46730.pdf");

  can3->Update();
  can3->SaveAs("figure/46730/TrackEMcal_KFP_DCA_run46730.pdf");

  can4->Update();
  can4->SaveAs("figure/46730/TrackEMcal_KFP_mass_radius_run46730.pdf");
}
