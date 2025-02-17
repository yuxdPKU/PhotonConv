#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>
#include <TArrow.h>

void drawsPHENIXInternal(TPaveText *pt);
void drawArrow(TArrow *arrow);

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void plot(int runnumber = 53744)
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TChain* chain = new TChain("tree");
  chain->Add(Form("eop_%d_kfp_v2.root",runnumber));

  TH1F* h1_epemc_dphi = new TH1F("h1_epemc_dphi","",100,-0.15,0.15);
  h1_epemc_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_epemc_dphi->GetYaxis()->SetTitle(Form("Events / %.2f",3./100.));
  h1_epemc_dphi->SetLineColor(kRed);
  //h1_epemc_dphi->GetXaxis()->SetLabelSize(0.05);

  TH1F* h1_ememc_dphi = new TH1F("h1_ememc_dphi","",100,-0.15,0.15);
  h1_ememc_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_ememc_dphi->GetYaxis()->SetTitle(Form("Events / %.2f",3./100.));
  h1_ememc_dphi->SetLineColor(kBlue);
  //h1_ememc_dphi->GetXaxis()->SetLabelSize(0.05);

  TLegend *legend1 = new TLegend(0.22, 0.75, 0.42, 0.9);
  legend1->AddEntry(h1_epemc_dphi, "e^{+}", "l");
  legend1->AddEntry(h1_ememc_dphi, "e^{-}", "l");

  TPaveText *pt1 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt1);

  TArrow *arrow1l = new TArrow(-0.05,100,-0.05,0,0.05,">");
  drawArrow(arrow1l);
  TArrow *arrow1r = new TArrow(0.10,100,0.10,0,0.05,">");
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

  TPaveText *pt2 = new TPaveText(.20, .60, .50, .75, "NDC");
  drawsPHENIXInternal(pt2);

  TArrow *arrow2l = new TArrow(-3,100,-3,0,0.05,">");
  drawArrow(arrow2l);
  TArrow *arrow2r = new TArrow(5,100,5,0,0.05,">");
  drawArrow(arrow2r);

  TH1F* h1_track12_deta = new TH1F("h1_track12_deta","",100,-0.2,0.2);
  h1_track12_deta->GetXaxis()->SetTitle("#Delta#eta");
  h1_track12_deta->GetYaxis()->SetTitle(Form("Events / %.3f",0.4/100.));
  h1_track12_deta->SetLineColor(kBlack);

  TPaveText *pt3 = new TPaveText(.22, .60, .52, .75, "NDC");
  drawsPHENIXInternal(pt3);

  TArrow *arrow3l = new TArrow(-0.02,100,-0.02,0,0.05,">");
  drawArrow(arrow3l);
  TArrow *arrow3r = new TArrow(0.02,100,0.02,0,0.05,">");
  drawArrow(arrow3r);

  TH1F* h1_track12_dphi = new TH1F("h1_track12_dphi","",100,-0.2,0.2);
  h1_track12_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_track12_dphi->GetYaxis()->SetTitle(Form("Events / %.3f",0.4/100.));
  h1_track12_dphi->SetLineColor(kBlack);

  TPaveText *pt4 = new TPaveText(.22, .60, .52, .75, "NDC");
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

  TLegend *legend7 = new TLegend(0.70, 0.6, 1.0, 0.7);
  legend7->AddEntry(h1_eop_ep, "e^{+}", "l");
  legend7->AddEntry(h1_eop_em, "e^{-}", "l");

  TPaveText *pt7 = new TPaveText(.70, .70, 1.0, 1., "NDC");
  drawsPHENIXInternal(pt7);

  TH1F* h1_eop_ep_cons = new TH1F("h1_eop_ep_cons","",40,0,2);
  h1_eop_ep_cons->GetXaxis()->SetTitle("E/p");
  h1_eop_ep_cons->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_ep_cons->SetLineColor(kRed);

  TH1F* h1_eop_em_cons = new TH1F("h1_eop_em_cons","",40,0,2);
  h1_eop_em_cons->GetXaxis()->SetTitle("E/p");
  h1_eop_em_cons->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_em_cons->SetLineColor(kBlue);

  TLegend *legend8 = new TLegend(0.70, 0.6, 1.0, 0.7);
  legend8->AddEntry(h1_eop_ep_cons, "e^{+}", "l");
  legend8->AddEntry(h1_eop_em_cons, "e^{-}", "l");

  TPaveText *pt8 = new TPaveText(.70, .70, 1.0, 1., "NDC");
  drawsPHENIXInternal(pt8);
  pt8->AddText("Require partner E/p#in[0.8,1.3]");

  TH2F* h2_eop_p_ep = new TH2F("h2_eop_p_ep","",20,0,5,40,0,2);
  h2_eop_p_ep->GetXaxis()->SetTitle("p [GeV/#it{c}]");
  h2_eop_p_ep->GetYaxis()->SetTitle("E/p");
  h2_eop_p_ep->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53744, 3 Million");
  TPaveText *pt9 = new TPaveText(.70, .70, 1.0, 1., "NDC");
  drawsPHENIXInternal(pt9);
  pt9->AddText("e^{+}");

  TH2F* h2_eop_p_em = new TH2F("h2_eop_p_em","",20,0,5,40,0,2);
  h2_eop_p_em->GetXaxis()->SetTitle("p [GeV/#it{c}]");
  h2_eop_p_em->GetYaxis()->SetTitle("E/p");
  h2_eop_p_em->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53744, 3 Million");
  TPaveText *pt10 = new TPaveText(.70, .70, 1.0, 1., "NDC");
  drawsPHENIXInternal(pt10);
  pt10->AddText("e^{-}");

  chain->Draw("_epemc_dphi>>h1_epemc_dphi");
  chain->Draw("_ememc_dphi>>h1_ememc_dphi");
  chain->Draw("_epemc_dz>>h1_epemc_dz");
  chain->Draw("_ememc_dz>>h1_ememc_dz");
  /*
  chain->Draw("_track12_deta>>h1_track12_deta");
  chain->Draw("_track12_dphi>>h1_track12_dphi");
  chain->Draw("_gamma_mass>>h1_gamma_mass","fabs(_track12_deta)<0.02");
  chain->Draw("_SV_radius>>h1_SV_radius","fabs(_track12_deta)<0.02");
  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<1 && _SV_radius>0 && _SV_radius<50");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<1 && _SV_radius>0 && _SV_radius<50");
  chain->Draw("_eop_ep>>h1_eop_ep_cons","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<1 && _SV_radius>0 && _SV_radius<50");
  chain->Draw("_eop_em>>h1_eop_em_cons","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<1 && _SV_radius>0 && _SV_radius<50");
  */

  chain->Draw("_track12_deta>>h1_track12_deta","_ep_emcal_e/_ep_p>0.8 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.8 && _em_emcal_e/_em_p<1.3 && _SV_radius>5 && _SV_radius<30");
  chain->Draw("_track12_dphi>>h1_track12_dphi","_ep_emcal_e/_ep_p>0.8 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.8 && _em_emcal_e/_em_p<1.3 && _SV_radius>5 && _SV_radius<30");

  chain->Draw("_gamma_mass>>h1_gamma_mass","fabs(_track12_deta)<0.02 && _ep_emcal_e/_ep_p>0.8 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.8 && _em_emcal_e/_em_p<1.3 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02");
  chain->Draw("_SV_radius>>h1_SV_radius","fabs(_track12_deta)<0.02 && _ep_emcal_e/_ep_p>0.8 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.8 && _em_emcal_e/_em_p<1.3 && _gamma_mass>0 && _gamma_mass<0.2 && fabs(_track12_dphi)<0.02");
  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<0.2 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<0.2 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02");
  chain->Draw("_eop_ep>>h1_eop_ep_cons","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<0.2 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02");
  chain->Draw("_eop_em>>h1_eop_em_cons","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<0.2 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02");
  chain->Draw("_ep_emcal_e/_ep_p:_ep_p>>h2_eop_p_ep","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<0.2 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02","colz");
  chain->Draw("_em_emcal_e/_em_p:_em_p>>h2_eop_p_em","fabs(_track12_deta)<0.02 && _gamma_mass>0 && _gamma_mass<0.2 && _SV_radius>5 && _SV_radius<30 && fabs(_track12_dphi)<0.02","colz");

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
  h1_eop_em_cons->Scale(h1_eop_ep_cons->Integral()/h1_eop_em_cons->Integral());
  h1_eop_ep_cons->SetMaximum(h1_eop_ep_cons->GetMaximum() > h1_eop_em_cons->GetMaximum() ? 1.1*h1_eop_ep_cons->GetMaximum() : 1.1*h1_eop_em_cons->GetMaximum());
  h1_eop_ep_cons->SetMinimum(0);
  h1_eop_ep_cons->Draw("hist");
  h1_eop_em_cons->Draw("same,hist");
  pt8->Draw("same");
  legend8->Draw("same");
  can8->Update();
  can8->SaveAs("./figure/eop_cons.pdf");

  TCanvas *can9 = new TCanvas("can9", "", 800, 800);
  can9->cd();
  h2_eop_p_ep->Draw("colz");
  pt9->Draw("same");
  can9->Update();
  can9->SaveAs("./figure/eop_p_ep.pdf");

  TCanvas *can10 = new TCanvas("can10", "", 800, 800);
  can10->cd();
  h2_eop_p_em->Draw("colz");
  pt10->Draw("same");
  can10->Update();
  can10->SaveAs("./figure/eop_p_em.pdf");

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
  pt->AddText(Form("Run 53744, 3 Million"));
}

void drawArrow(TArrow *arrow)
{
  arrow->SetLineColor(2);
  arrow->SetLineWidth(2);
}
