#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>
#include <TArrow.h>

void drawsPHENIXInternal(TPaveText *pt);
void drawArrow(TArrow *arrow);

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void plot_recotruth(std::string inputfile = "eop_15_kfp_v2.root")
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree");
  chain->Add("eop_15_kfp_v2.root");

// track E reco truth ratio
  TH1F* h1_ep_pE_recotruth_ratio = new TH1F("h1_ep_pE_recotruth_ratio","",100,0,2);
  h1_ep_pE_recotruth_ratio->GetXaxis()->SetTitle("Track #it{E} reco/truth ratio [GeV]");
  h1_ep_pE_recotruth_ratio->GetYaxis()->SetTitle(Form("Events / %.2f [GeV]",2./100.));
  h1_ep_pE_recotruth_ratio->SetLineColor(kRed);

  TH1F* h1_em_pE_recotruth_ratio = new TH1F("h1_em_pE_recotruth_ratio","",100,0,2);
  h1_em_pE_recotruth_ratio->GetXaxis()->SetTitle("Track #it{E} reco/truth ratio [GeV]");
  h1_em_pE_recotruth_ratio->GetYaxis()->SetTitle(Form("Events / %.2f [GeV]",3./100.));
  h1_em_pE_recotruth_ratio->SetLineColor(kBlue);

  TLegend *legend1 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend1->AddEntry(h1_ep_pE_recotruth_ratio, "e^{+}", "l");
  legend1->AddEntry(h1_em_pE_recotruth_ratio, "e^{-}", "l");

  TPaveText *pt1 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt1);

  chain->Draw("_ep_pE/_true_ep_pE>>h1_ep_pE_recotruth_ratio","fabs(_track12_deta)<0.02");
  chain->Draw("_em_pE/_true_em_pE>>h1_em_pE_recotruth_ratio","fabs(_track12_deta)<0.02");

  TCanvas *can1 = new TCanvas("can1", "", 800, 800);
  can1->cd();
  h1_em_pE_recotruth_ratio->Scale(h1_ep_pE_recotruth_ratio->Integral()/h1_em_pE_recotruth_ratio->Integral());
  h1_ep_pE_recotruth_ratio->SetMaximum(h1_ep_pE_recotruth_ratio->GetMaximum() > h1_em_pE_recotruth_ratio->GetMaximum() ? 1.1*h1_ep_pE_recotruth_ratio->GetMaximum() : 1.1*h1_em_pE_recotruth_ratio->GetMaximum());
  h1_ep_pE_recotruth_ratio->SetMinimum(0);
  h1_ep_pE_recotruth_ratio->Draw("hist");
  h1_em_pE_recotruth_ratio->Draw("same,hist");
  legend1->Draw("same");
  pt1->Draw("same");
  can1->Update();
  can1->SaveAs("./figure/track_pE_recotruth_ratio.pdf");


// track phi reco truth diff
  TH1F* h1_ep_phi_recotruth_diff = new TH1F("h1_ep_phi_recotruth_diff","",100,-1,1);
  h1_ep_phi_recotruth_diff->GetXaxis()->SetTitle("Track #phi reco - truth [rad]");
  h1_ep_phi_recotruth_diff->GetYaxis()->SetTitle(Form("Events / %.2f",2./100.));
  h1_ep_phi_recotruth_diff->SetLineColor(kRed);

  TH1F* h1_em_phi_recotruth_diff = new TH1F("h1_em_phi_recotruth_diff","",100,-1,1);
  h1_em_phi_recotruth_diff->GetXaxis()->SetTitle("Track #phi reco - truth [rad]");
  h1_em_phi_recotruth_diff->GetYaxis()->SetTitle(Form("Events / %.2f",2./100.));
  h1_em_phi_recotruth_diff->SetLineColor(kBlue);

  TLegend *legend2 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend2->AddEntry(h1_ep_phi_recotruth_diff, "e^{+}", "l");
  legend2->AddEntry(h1_em_phi_recotruth_diff, "e^{-}", "l");

  TPaveText *pt2 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt2);

  chain->Draw("_ep_phi_recotruthDiff>>h1_ep_phi_recotruth_diff","fabs(_track12_deta)<0.02");
  chain->Draw("_em_phi_recotruthDiff>>h1_em_phi_recotruth_diff","fabs(_track12_deta)<0.02");

  TCanvas *can2 = new TCanvas("can2", "", 800, 800);
  can2->cd();
  h1_em_phi_recotruth_diff->Scale(h1_ep_phi_recotruth_diff->Integral()/h1_em_phi_recotruth_diff->Integral());
  h1_ep_phi_recotruth_diff->SetMaximum(h1_ep_phi_recotruth_diff->GetMaximum() > h1_em_phi_recotruth_diff->GetMaximum() ? 1.1*h1_ep_phi_recotruth_diff->GetMaximum() : 1.1*h1_em_phi_recotruth_diff->GetMaximum());
  h1_ep_phi_recotruth_diff->SetMinimum(0);
  h1_ep_phi_recotruth_diff->Draw("hist");
  h1_em_phi_recotruth_diff->Draw("same,hist");
  legend1->Draw("same");
  pt2->Draw("same");
  can2->Update();
  can2->SaveAs("./figure/track_phi_recotruth_diff.pdf");


// track eta reco truth diff
  TH1F* h1_ep_eta_recotruth_diff = new TH1F("h1_ep_eta_recotruth_diff","",100,-0.5,0.5);
  h1_ep_eta_recotruth_diff->GetXaxis()->SetTitle("Track #eta reco - truth");
  h1_ep_eta_recotruth_diff->GetYaxis()->SetTitle(Form("Events / %.2f",1./100.));
  h1_ep_eta_recotruth_diff->SetLineColor(kRed);

  TH1F* h1_em_eta_recotruth_diff = new TH1F("h1_em_eta_recotruth_diff","",100,-0.5,0.5);
  h1_em_eta_recotruth_diff->GetXaxis()->SetTitle("Track #eta reco - truth");
  h1_em_eta_recotruth_diff->GetYaxis()->SetTitle(Form("Events / %.2f",1/100.));
  h1_em_eta_recotruth_diff->SetLineColor(kBlue);

  TLegend *legend3 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend3->AddEntry(h1_ep_eta_recotruth_diff, "e^{+}", "l");
  legend3->AddEntry(h1_em_eta_recotruth_diff, "e^{-}", "l");

  TPaveText *pt3 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt3);

  chain->Draw("_ep_eta_recotruthDiff>>h1_ep_eta_recotruth_diff","fabs(_track12_deta)<0.02");
  chain->Draw("_em_eta_recotruthDiff>>h1_em_eta_recotruth_diff","fabs(_track12_deta)<0.02");

  TCanvas *can3 = new TCanvas("can3", "", 800, 800);
  can3->cd();
  h1_em_eta_recotruth_diff->Scale(h1_ep_eta_recotruth_diff->Integral()/h1_em_eta_recotruth_diff->Integral());
  h1_ep_eta_recotruth_diff->SetMaximum(h1_ep_eta_recotruth_diff->GetMaximum() > h1_em_eta_recotruth_diff->GetMaximum() ? 1.1*h1_ep_eta_recotruth_diff->GetMaximum() : 1.1*h1_em_eta_recotruth_diff->GetMaximum());
  h1_ep_eta_recotruth_diff->SetMinimum(0);
  h1_ep_eta_recotruth_diff->Draw("hist");
  h1_em_eta_recotruth_diff->Draw("same,hist");
  legend1->Draw("same");
  pt2->Draw("same");
  can3->Update();
  can3->SaveAs("./figure/track_eta_recotruth_diff.pdf");


// calo E reco truth ratio
  TH1F* h1_epemc_E_recotruth_ratio = new TH1F("h1_epemc_E_recotruth_ratio","",100,0,2);
  h1_epemc_E_recotruth_ratio->GetXaxis()->SetTitle("Calo #it{E} reco/truth ratio [GeV]");
  h1_epemc_E_recotruth_ratio->GetYaxis()->SetTitle(Form("Events / %.2f [GeV]",2./100.));
  h1_epemc_E_recotruth_ratio->SetLineColor(kRed);

  TH1F* h1_ememc_E_recotruth_ratio = new TH1F("h1_ememc_E_recotruth_ratio","",100,0,2);
  h1_ememc_E_recotruth_ratio->GetXaxis()->SetTitle("Calo #it{E} reco/truth ratio [GeV]");
  h1_ememc_E_recotruth_ratio->GetYaxis()->SetTitle(Form("Events / %.2f [GeV]",3./100.));
  h1_ememc_E_recotruth_ratio->SetLineColor(kBlue);

  TLegend *legend4 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend4->AddEntry(h1_epemc_E_recotruth_ratio, "e^{+}", "l");
  legend4->AddEntry(h1_ememc_E_recotruth_ratio, "e^{-}", "l");

  TPaveText *pt4 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt4);

  chain->Draw("_ep_emcal_e/_true_ep_pE>>h1_epemc_E_recotruth_ratio","fabs(_track12_deta)<0.02");
  chain->Draw("_em_emcal_e/_true_em_pE>>h1_ememc_E_recotruth_ratio","fabs(_track12_deta)<0.02");

  TCanvas *can4 = new TCanvas("can4", "", 800, 800);
  can4->cd();
  h1_ememc_E_recotruth_ratio->Scale(h1_epemc_E_recotruth_ratio->Integral()/h1_ememc_E_recotruth_ratio->Integral());
  h1_epemc_E_recotruth_ratio->SetMaximum(h1_epemc_E_recotruth_ratio->GetMaximum() > h1_ememc_E_recotruth_ratio->GetMaximum() ? 1.1*h1_epemc_E_recotruth_ratio->GetMaximum() : 1.1*h1_ememc_E_recotruth_ratio->GetMaximum());
  h1_epemc_E_recotruth_ratio->SetMinimum(0);
  h1_epemc_E_recotruth_ratio->Draw("hist");
  h1_ememc_E_recotruth_ratio->Draw("same,hist");
  legend4->Draw("same");
  pt4->Draw("same");
  can4->Update();
  can4->SaveAs("./figure/calo_E_recotruth_ratio.pdf");


// SV phi reco truth diff
  TH1F* h1_sv_phi_recotruth_diff = new TH1F("h1_sv_phi_recotruth_diff","",100,-0.5,0.5);
  h1_sv_phi_recotruth_diff->GetXaxis()->SetTitle("SV #phi reco - truth [rad]");
  h1_sv_phi_recotruth_diff->GetYaxis()->SetTitle(Form("Events / %.2f",1./100.));
  h1_sv_phi_recotruth_diff->SetLineColor(kBlack);

  TPaveText *pt5 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt5);

  chain->Draw("_SV_phi_recotruthDiff>>h1_sv_phi_recotruth_diff","fabs(_track12_deta)<0.02");

  TCanvas *can5 = new TCanvas("can5", "", 800, 800);
  can5->cd();
  h1_sv_phi_recotruth_diff->SetMinimum(0);
  h1_sv_phi_recotruth_diff->Draw("hist");
  pt5->Draw("same");
  can5->Update();
  can5->SaveAs("./figure/sv_phi_recotruth_diff.pdf");


// SV phi reco truth diff
  TH1F* h1_sv_z_recotruth_diff = new TH1F("h1_sv_z_recotruth_diff","",100,-30,30);
  h1_sv_z_recotruth_diff->GetXaxis()->SetTitle("SV z reco - truth [cm]");
  h1_sv_z_recotruth_diff->GetYaxis()->SetTitle(Form("Events / %.1f [cm]",60./100.));
  h1_sv_z_recotruth_diff->SetLineColor(kBlack);

  TPaveText *pt6 = new TPaveText(.17, .70, .47, .85, "NDC");
  drawsPHENIXInternal(pt6);

  chain->Draw("_SV_z_recotruthDiff>>h1_sv_z_recotruth_diff","fabs(_track12_deta)<0.02");

  TCanvas *can6 = new TCanvas("can6", "", 800, 800);
  can6->cd();
  h1_sv_z_recotruth_diff->SetMinimum(0);
  h1_sv_z_recotruth_diff->Draw("hist");
  pt6->Draw("same");
  can6->Update();
  can6->SaveAs("./figure/sv_z_recotruth_diff.pdf");



}

void drawsPHENIXInternal(TPaveText *pt)
{
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt->AddText("p+p #sqrt{s}=200 GeV");
  pt->AddText(Form("MDC2 PhotonJet "));
  pt->AddText(Form("pTmin=5GeV with pileup"));
}

void drawArrow(TArrow *arrow)
{
  arrow->SetLineColor(2);
  arrow->SetLineWidth(2);
}
