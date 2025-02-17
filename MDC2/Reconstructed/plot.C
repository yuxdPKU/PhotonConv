#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>
#include <TArrow.h>

void drawsPHENIXInternal(TPaveText *pt);
void drawArrow(TArrow *arrow);

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void plot(std::string inputfile = "eop_15_kfp_v2.root")
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree");
  chain->Add("eop_15_kfp_v2.root");

  TH1F* h1_epemc_dphi = new TH1F("h1_epemc_dphi","",100,-0.15,0.15);
  h1_epemc_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_epemc_dphi->GetYaxis()->SetTitle(Form("Events / %.2f",3./100.));
  h1_epemc_dphi->SetLineColor(kRed);
  h1_epemc_dphi->GetXaxis()->SetLabelSize(0.05);

  TH1F* h1_ememc_dphi = new TH1F("h1_ememc_dphi","",100,-0.15,0.15);
  h1_ememc_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_ememc_dphi->GetYaxis()->SetTitle(Form("Events / %.2f",3./100.));
  h1_ememc_dphi->SetLineColor(kBlue);
  h1_ememc_dphi->GetXaxis()->SetLabelSize(0.05);

  TLegend *legend1 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend1->AddEntry(h1_epemc_dphi, "e^{+}", "l");
  legend1->AddEntry(h1_ememc_dphi, "e^{-}", "l");

  TPaveText *pt1 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt1);

  TArrow *arrow1l = new TArrow(-0.03,1000,-0.03,0,0.05,">");
  drawArrow(arrow1l);
  TArrow *arrow1r = new TArrow(0.03,1000,0.03,0,0.05,">");
  drawArrow(arrow1r);

  TH1F* h1_epemc_dz = new TH1F("h1_epemc_dz","",100,-10,10);
  h1_epemc_dz->GetXaxis()->SetTitle("#Deltaz [cm]");
  h1_epemc_dz->GetYaxis()->SetTitle(Form("Events / %.1f [cm]",20./100.));
  h1_epemc_dz->SetLineColor(kRed);

  TH1F* h1_ememc_dz = new TH1F("h1_ememc_dz","",100,-10,10);
  h1_ememc_dz->GetXaxis()->SetTitle("#Deltaz [cm]");
  h1_ememc_dz->GetYaxis()->SetTitle(Form("Events / %.1f [cm]",20./100.));
  h1_ememc_dz->SetLineColor(kBlue);

  TLegend *legend2 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend2->AddEntry(h1_epemc_dz, "e^{+}", "l");
  legend2->AddEntry(h1_ememc_dz, "e^{-}", "l");

  TPaveText *pt2 = new TPaveText(.17, .60, .47, .75, "NDC");
  drawsPHENIXInternal(pt2);

  TArrow *arrow2l = new TArrow(-5,200,-5,0,0.05,">");
  drawArrow(arrow2l);
  TArrow *arrow2r = new TArrow(5,200,5,0,0.05,">");
  drawArrow(arrow2r);

  TH1F* h1_track12_deta = new TH1F("h1_track12_deta","",100,-0.2,0.2);
  h1_track12_deta->GetXaxis()->SetTitle("#Delta#eta");
  h1_track12_deta->GetYaxis()->SetTitle(Form("Events / %.3f",0.4/100.));
  h1_track12_deta->SetLineColor(kBlack);

  TPaveText *pt3 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt3);

  TArrow *arrow3l = new TArrow(-0.02,200,-0.02,0,0.05,">");
  drawArrow(arrow3l);
  TArrow *arrow3r = new TArrow(0.02,200,0.02,0,0.05,">");
  drawArrow(arrow3r);

  TH1F* h1_track12_dphi = new TH1F("h1_track12_dphi","",100,-0.2,0.2);
  h1_track12_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_track12_dphi->GetYaxis()->SetTitle(Form("Events / %.3f",0.4/100.));
  h1_track12_dphi->SetLineColor(kBlack);

  TPaveText *pt4 = new TPaveText(.22, .70, .52, .85, "NDC");
  drawsPHENIXInternal(pt4);

  TH1F* h1_gamma_mass = new TH1F("h1_gamma_mass","",100,0,1);
  h1_gamma_mass->GetXaxis()->SetTitle("#it{M}_{e^{+}e^{-}} [GeV/#it{c}^{2}]");
  h1_gamma_mass->GetYaxis()->SetTitle(Form("Events / %.2f [GeV/#it{c}^{2}]",1/100.));
  h1_gamma_mass->SetLineColor(kBlack);

  TPaveText *pt5 = new TPaveText(.62, .72, .92, .92, "NDC");
  drawsPHENIXInternal(pt5);

  TH1F* h1_SV_radius = new TH1F("h1_SV_radius","",50,0,50);
  h1_SV_radius->GetXaxis()->SetTitle("#it{R}_{e^{+}e^{-}} [cm]");
  h1_SV_radius->GetYaxis()->SetTitle(Form("Events / %.1f [cm]",50./50.));
  h1_SV_radius->SetLineColor(kBlack);

  TPaveText *pt6 = new TPaveText(.62, .72, .92, .92, "NDC");
  drawsPHENIXInternal(pt6);

  TH1F* h1_eop_ep = new TH1F("h1_eop_ep","",40,0,2);
  h1_eop_ep->GetXaxis()->SetTitle("E/p");
  h1_eop_ep->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_ep->SetLineColor(kRed);

  TH1F* h1_eop_em = new TH1F("h1_eop_em","",40,0,2);
  h1_eop_em->GetXaxis()->SetTitle("E/p");
  h1_eop_em->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_em->SetLineColor(kBlue);

  TLegend *legend7 = new TLegend(0.17, 0.75, 0.47, 0.9);
  legend7->AddEntry(h1_eop_ep, "e^{+}", "l");
  legend7->AddEntry(h1_eop_em, "e^{-}", "l");

  TPaveText *pt7 = new TPaveText(.17, .52, .47, .72, "NDC");
  drawsPHENIXInternal(pt7);

  TH1F* h1_true_gamma_decay_radius = new TH1F("h1_true_gamma_decay_radius","",100,0,25);
  h1_true_gamma_decay_radius->GetXaxis()->SetTitle("True #it{R}_{e^{+}e^{-}} [cm]");
  h1_true_gamma_decay_radius->GetYaxis()->SetTitle(Form("Events / %.1f [cm]",50./50.));
  h1_true_gamma_decay_radius->SetLineColor(kBlack);

  TPaveText *pt8 = new TPaveText(.62, .72, .92, .92, "NDC");
  drawsPHENIXInternal(pt8);

  chain->Draw("_epemc_dphi>>h1_epemc_dphi");
  chain->Draw("_ememc_dphi>>h1_ememc_dphi");
  chain->Draw("_epemc_dz>>h1_epemc_dz");
  chain->Draw("_ememc_dz>>h1_ememc_dz");
  chain->Draw("_track12_deta>>h1_track12_deta");
  //chain->Draw("_track12_dphi>>h1_track12_dphi");
  chain->Draw("_track12_dphi>>h1_track12_dphi","fabs(_track12_deta)<0.02 && fabs(_SV_radius-_true_gamma_decay_radius)<5");
  chain->Draw("_gamma_mass>>h1_gamma_mass","fabs(_track12_deta)<0.02");
  chain->Draw("_SV_radius>>h1_SV_radius","fabs(_track12_deta)<0.02");
  chain->Draw("_true_gamma_decay_radius>>h1_true_gamma_decay_radius");
  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<1 && _SV_radius>0 && _SV_radius<50");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<1 && _SV_radius>0 && _SV_radius<50");

  TCanvas *can1 = new TCanvas("can1", "", 800, 800);
  can1->cd();
  h1_ememc_dphi->Scale(h1_epemc_dphi->Integral()/h1_ememc_dphi->Integral());
  h1_epemc_dphi->SetMaximum(h1_epemc_dphi->GetMaximum() > h1_ememc_dphi->GetMaximum() ? 1.1*h1_epemc_dphi->GetMaximum() : 1.1*h1_ememc_dphi->GetMaximum());
  h1_epemc_dphi->Draw("hist");
  h1_ememc_dphi->Draw("same,hist");
  legend1->Draw("same");
  pt1->Draw("same");
  arrow1l->Draw("same");
  arrow1r->Draw("same");
  can1->Update();
  can1->SaveAs("./figure/dphi.pdf");

  TCanvas *can2 = new TCanvas("can2", "", 800, 800);
  can2->cd();
  h1_ememc_dz->Scale(h1_epemc_dz->Integral()/h1_ememc_dz->Integral());
  h1_epemc_dz->SetMaximum(h1_epemc_dz->GetMaximum() > h1_ememc_dz->GetMaximum() ? 1.1*h1_epemc_dz->GetMaximum() : 1.1*h1_ememc_dz->GetMaximum());
  h1_epemc_dz->Draw("hist");
  h1_ememc_dz->Draw("same,hist");
  legend2->Draw("same");
  pt2->Draw("same");
  arrow2l->Draw("same");
  arrow2r->Draw("same");
  can2->Update();
  can2->SaveAs("./figure/dz.pdf");

  TCanvas *can3 = new TCanvas("can3", "", 800, 800);
  can3->cd();
  h1_track12_deta->SetMinimum(0);
  h1_track12_deta->Draw("hist");
  pt3->Draw("same");
  arrow3l->Draw("same");
  arrow3r->Draw("same");
  can3->Update();
  can3->SaveAs("./figure/track12_deta.pdf");

  TCanvas *can4 = new TCanvas("can4", "", 800, 800);
  can4->cd();
  h1_track12_dphi->SetMinimum(0);
  h1_track12_dphi->Draw("hist");
  pt4->Draw("same");
  can4->Update();
  can4->SaveAs("./figure/track12_dphi.pdf");

  TCanvas *can5 = new TCanvas("can5", "", 800, 800);
  can5->cd();
  h1_gamma_mass->SetMinimum(0);
  h1_gamma_mass->Draw("hist");
  pt5->Draw("same");
  can5->Update();
  can5->SaveAs("./figure/gamma_mass.pdf");

  TCanvas *can6 = new TCanvas("can6", "", 800, 800);
  can6->cd();
  h1_SV_radius->SetMinimum(0);
  h1_SV_radius->Draw("hist");
  pt6->Draw("same");
  can6->Update();
  can6->SaveAs("./figure/gamma_radius.pdf");

  TCanvas *can7 = new TCanvas("can7", "", 800, 800);
  can7->cd();
  h1_eop_em->Scale(h1_eop_ep->Integral()/h1_eop_em->Integral());
  h1_eop_ep->SetMaximum(h1_eop_ep->GetMaximum() > h1_eop_em->GetMaximum() ? 1.1*h1_eop_ep->GetMaximum() : 1.1*h1_eop_em->GetMaximum());
  h1_eop_ep->SetMinimum(0);
  h1_eop_ep->Draw("hist");
  h1_eop_em->Draw("same,hist");
  pt7->Draw("same");
  legend7->Draw("same");
  can7->Update();
  can7->SaveAs("./figure/eop.pdf");

  TCanvas *can8 = new TCanvas("can8", "", 800, 800);
  can8->cd();
  h1_true_gamma_decay_radius->SetMaximum(8000);
  h1_true_gamma_decay_radius->SetMinimum(0);
  h1_true_gamma_decay_radius->Draw("hist");
  pt8->Draw("same");
  can8->Update();
  can8->SaveAs("./figure/true_gamma_radius.pdf");


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
