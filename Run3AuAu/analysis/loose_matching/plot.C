#include <filesystem>
#include "../utilities.h"
//#include <sPhenixStyle.C>
#include <TArrow.h>

void drawsPHENIXInternal(TPaveText *pt);
void drawArrow(TArrow *arrow);
void fit_gauss(TH2* h, TString name, bool verbose=0, TPaveText *pt=nullptr);
TH2* DuplicateTH2(const TH2* originalHist);
void KeepSlopeBand(TH2* hist, double xMinAtY0 = 0.01, double xMaxAtY0 = 0.04);

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void plot()
{
  //SetsPhenixStyle();
  gStyle->SetOptStat(0);

  TGaxis::SetMaxDigits(3);

  TChain* chain = new TChain("tree");
  chain->Add("eop_kfp_unlikesign.root");
  //chain->Add("eop_kfp_likesign.root");

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

  TPaveText *pt2 = new TPaveText(.40, .30, .60, .45, "NDC");
  drawsPHENIXInternal(pt2);

  TArrow *arrow2l = new TArrow(-3,100,-3,0,0.05,">");
  drawArrow(arrow2l);
  TArrow *arrow2r = new TArrow(5,100,5,0,0.05,">");
  drawArrow(arrow2r);

  TH1F* h1_track12_deta = new TH1F("h1_track12_deta","",100,-0.2,0.2);
  h1_track12_deta->GetXaxis()->SetTitle("#Delta#eta");
  h1_track12_deta->GetYaxis()->SetTitle(Form("Events / %.3f",0.4/100.));
  h1_track12_deta->SetLineColor(kBlack);

  TPaveText *pt3 = new TPaveText(.17, .45, .47, .60, "NDC");
  drawsPHENIXInternal(pt3);

  TArrow *arrow3l = new TArrow(-0.02,8000,-0.02,0,0.05,">");
  drawArrow(arrow3l);
  TArrow *arrow3r = new TArrow(0.02,8000,0.02,0,0.05,">");
  drawArrow(arrow3r);

  TH1F* h1_track12_dphi = new TH1F("h1_track12_dphi","",100,-0.2,0.2);
  h1_track12_dphi->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h1_track12_dphi->GetYaxis()->SetTitle(Form("Events / %.3f",0.4/100.));
  h1_track12_dphi->SetLineColor(kBlack);

  TPaveText *pt4 = new TPaveText(.17, .60, .47, .75, "NDC");
  drawsPHENIXInternal(pt4);

  TArrow *arrow4l = new TArrow(-0.02,3000,-0.02,0,0.05,">");
  drawArrow(arrow4l);
  TArrow *arrow4r = new TArrow(0.02,3000,0.02,0,0.05,">");
  drawArrow(arrow4r);

  TH1F* h1_gamma_mass = new TH1F("h1_gamma_mass","",100,0,0.05);
  h1_gamma_mass->GetXaxis()->SetTitle("#it{M}_{e^{+}e^{-}} [GeV/#it{c}^{2}]");
  h1_gamma_mass->GetYaxis()->SetTitle(Form("Events / %.2f [GeV/#it{c}^{2}]",1/100.));
  h1_gamma_mass->SetLineColor(kBlack);

  //TPaveText *pt5 = new TPaveText(.62, .77, .92, .92, "NDC");
  TPaveText *pt5 = new TPaveText(.55, .57, .85, .72, "NDC");
  drawsPHENIXInternal(pt5);

  TH1F* h1_SV_radius = new TH1F("h1_SV_radius","",50,0,30);
  h1_SV_radius->GetXaxis()->SetTitle("#it{R}_{e^{+}e^{-}} [cm]");
  h1_SV_radius->GetYaxis()->SetTitle(Form("Events / %.1f [cm]",50./50.));
  h1_SV_radius->SetLineColor(kBlack);

  TPaveText *pt6 = new TPaveText(.62, .77, .92, .92, "NDC");
  drawsPHENIXInternal(pt6);

  TH1F* h1_eop_ep = new TH1F("h1_eop_ep","",40,0,2);
  h1_eop_ep->GetXaxis()->SetTitle("E/p");
  h1_eop_ep->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_ep->SetLineColor(kRed);

  TH1F* h1_eop_em = new TH1F("h1_eop_em","",40,0,2);
  h1_eop_em->GetXaxis()->SetTitle("E/p");
  h1_eop_em->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_em->SetLineColor(kBlue);

  TLegend *legend7 = new TLegend(0.75, 0.6, 0.92, 0.7);
  legend7->AddEntry(h1_eop_ep, "e^{+}", "l");
  legend7->AddEntry(h1_eop_em, "e^{-}", "l");

  TPaveText *pt7 = new TPaveText(.62, .77, .92, .92, "NDC");
  drawsPHENIXInternal(pt7);

  TH1F* h1_eop_ep_cons = new TH1F("h1_eop_ep_cons","",40,0,2);
  h1_eop_ep_cons->GetXaxis()->SetTitle("E/p");
  h1_eop_ep_cons->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_ep_cons->SetLineColor(kRed);

  TH1F* h1_eop_em_cons = new TH1F("h1_eop_em_cons","",40,0,2);
  h1_eop_em_cons->GetXaxis()->SetTitle("E/p");
  h1_eop_em_cons->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_em_cons->SetLineColor(kBlue);

  TLegend *legend8 = new TLegend(0.75, 0.6, 0.92, 0.7);
  legend8->AddEntry(h1_eop_ep_cons, "e^{+}", "l");
  legend8->AddEntry(h1_eop_em_cons, "e^{-}", "l");

  TPaveText *pt8 = new TPaveText(.62, .77, .92, .92, "NDC");
  drawsPHENIXInternal(pt8);
  //pt8->AddText("Require partner E/p#in[0.7,1.3]");

  TH2F* h2_eop_p_ep = new TH2F("h2_eop_p_ep","",20,0,5,40,0,3);
  h2_eop_p_ep->GetXaxis()->SetTitle("e^{+} p [GeV/#it{c}]");
  h2_eop_p_ep->GetYaxis()->SetTitle("e^{+} E/p");
  h2_eop_p_ep->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");
  TPaveText *pt9 = new TPaveText(.55, .77, .80, .92, "NDC");
  drawsPHENIXInternal(pt9);

  TH2F* h2_eop_p_em = new TH2F("h2_eop_p_em","",20,0,5,40,0,3);
  h2_eop_p_em->GetXaxis()->SetTitle("e^{-} p [GeV/#it{c}]");
  h2_eop_p_em->GetYaxis()->SetTitle("e^{-} E/p");
  h2_eop_p_em->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");
  TPaveText *pt10 = new TPaveText(.55, .77, .80, .92, "NDC");
  drawsPHENIXInternal(pt10);

  TH1F* h1_ep_crossing = new TH1F("h1_ep_crossing","",600,-100,400);
  h1_ep_crossing->GetXaxis()->SetTitle("Crossing");
  h1_ep_crossing->GetYaxis()->SetTitle(Form("Events"));
  h1_ep_crossing->SetLineColor(kRed);

  TH1F* h1_em_crossing = new TH1F("h1_em_crossing","",600,-100,400);
  h1_em_crossing->GetXaxis()->SetTitle("Crossing");
  h1_em_crossing->GetYaxis()->SetTitle(Form("Events"));
  h1_em_crossing->SetLineColor(kBlue);

  TLegend *legend11 = new TLegend(0.75, 0.6, 0.92, 0.7);
  legend11->AddEntry(h1_ep_crossing, "e^{+}", "l");
  legend11->AddEntry(h1_em_crossing, "e^{-}", "l");

  TPaveText *pt11 = new TPaveText(.62, .77, .92, .92, "NDC");
  drawsPHENIXInternal(pt11);

  TH1F* h1_ep_nsi = new TH1F("h1_ep_nsi","",8,0,8);
  h1_ep_nsi->GetXaxis()->SetTitle("nmvtx+nintt");
  h1_ep_nsi->GetYaxis()->SetTitle(Form("Events"));
  h1_ep_nsi->SetLineColor(kRed);

  TH1F* h1_em_nsi = new TH1F("h1_em_nsi","",8,0,8);
  h1_em_nsi->GetXaxis()->SetTitle("nmvtx+nintt");
  h1_em_nsi->GetYaxis()->SetTitle(Form("Events"));
  h1_em_nsi->SetLineColor(kBlue);

  TLegend *legend12 = new TLegend(0.75, 0.6, 0.92, 0.7);
  legend12->AddEntry(h1_ep_nsi, "e^{+}", "l");
  legend12->AddEntry(h1_em_nsi, "e^{-}", "l");

  TPaveText *pt12 = new TPaveText(.62, .77, .92, .92, "NDC");
  drawsPHENIXInternal(pt12);

  TH2F* h2_eop_qopt_ep = new TH2F("h2_eop_qopt_ep","",20,-3,3,40,0,3);
  h2_eop_qopt_ep->GetXaxis()->SetTitle("q/pT [(GeV/#it{c})^{-1}]");
  h2_eop_qopt_ep->GetYaxis()->SetTitle("E/p");
  h2_eop_qopt_ep->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");
  TPaveText *pt13 = new TPaveText(.55, .77, .80, .92, "NDC");
  drawsPHENIXInternal(pt13);

  TH2F* h2_eop_qopt_em = new TH2F("h2_eop_qopt_em","",20,-3,3,40,0,3);
  h2_eop_qopt_em->GetXaxis()->SetTitle("q/pT [(GeV/#it{c})^{-1}]");
  h2_eop_qopt_em->GetYaxis()->SetTitle("E/p");
  h2_eop_qopt_em->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");
  TPaveText *pt14 = new TPaveText(.55, .77, .80, .92, "NDC");
  drawsPHENIXInternal(pt14);

  TH2F* h2_eop_relemcalphi = new TH2F("h2_eop_relemcalphi","",16,0,2*TMath::Pi()/128.,40,0,2);
  h2_eop_relemcalphi->GetXaxis()->SetNoExponent(true);
  h2_eop_relemcalphi->GetXaxis()->SetTitle("fmod(EMCal #phi, 2#pi/128) [rad]");
  h2_eop_relemcalphi->GetYaxis()->SetTitle("E/p");
  h2_eop_relemcalphi->GetXaxis()->SetNdivisions(505);
  h2_eop_relemcalphi->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TH2F* h2_eop_reltrackphi = new TH2F("h2_eop_reltrackphi","",16,0,2*TMath::Pi()/128.,40,0,2);
  h2_eop_reltrackphi->GetXaxis()->SetNoExponent(true);
  h2_eop_reltrackphi->GetXaxis()->SetTitle("fmod(Track projected #phi, 2#pi/128) [rad]");
  h2_eop_reltrackphi->GetYaxis()->SetTitle("E/p");
  h2_eop_reltrackphi->GetXaxis()->SetNdivisions(505);
  h2_eop_reltrackphi->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TH2F* h2_eop_reltrackphi_ep = new TH2F("h2_eop_reltrackphi_ep","",16,0,2*TMath::Pi()/128.,40,0,2);
  h2_eop_reltrackphi_ep->GetXaxis()->SetNoExponent(true);
  h2_eop_reltrackphi_ep->GetXaxis()->SetTitle("fmod(e^{+} Track projected #phi, 2#pi/128) [rad]");
  h2_eop_reltrackphi_ep->GetYaxis()->SetTitle("E/p");
  h2_eop_reltrackphi_ep->GetXaxis()->SetNdivisions(505);
  h2_eop_reltrackphi_ep->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TH2F* h2_eop_reltrackphi_em = new TH2F("h2_eop_reltrackphi_em","",16,0,2*TMath::Pi()/128.,40,0,2);
  h2_eop_reltrackphi_em->GetXaxis()->SetNoExponent(true);
  h2_eop_reltrackphi_em->GetXaxis()->SetTitle("fmod(e^{-} Track projected #phi, 2#pi/128) [rad]");
  h2_eop_reltrackphi_em->GetYaxis()->SetTitle("E/p");
  h2_eop_reltrackphi_em->GetXaxis()->SetNdivisions(505);
  h2_eop_reltrackphi_em->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TH2F* h2_relemcalphi_reltrackphi = new TH2F("h2_relemcalphi_reltrackphi","",16,0,2*TMath::Pi()/128.,16,0,2*TMath::Pi()/128.);
  h2_relemcalphi_reltrackphi->GetXaxis()->SetNoExponent(true);
  h2_relemcalphi_reltrackphi->GetYaxis()->SetNoExponent(true);
  h2_relemcalphi_reltrackphi->GetXaxis()->SetTitle("fmod(Track projected #phi, 2#pi/128) [rad]");
  h2_relemcalphi_reltrackphi->GetYaxis()->SetTitle("fmod(EMCal #phi, 2#pi/128) [rad]");
  h2_relemcalphi_reltrackphi->GetXaxis()->SetNdivisions(505);
  h2_relemcalphi_reltrackphi->GetYaxis()->SetNdivisions(505);
  h2_relemcalphi_reltrackphi->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TH2F* h2_relemcalphi_reltrackphi_ep = new TH2F("h2_relemcalphi_reltrackphi_ep","",16,0,2*TMath::Pi()/128.,16,0,2*TMath::Pi()/128.);
  h2_relemcalphi_reltrackphi_ep->GetXaxis()->SetNoExponent(true);
  h2_relemcalphi_reltrackphi_ep->GetYaxis()->SetNoExponent(true);
  h2_relemcalphi_reltrackphi_ep->GetXaxis()->SetTitle("fmod(e^{+} Track projected #phi, 2#pi/128) [rad]");
  h2_relemcalphi_reltrackphi_ep->GetYaxis()->SetTitle("fmod(EMCal #phi, 2#pi/128) [rad]");
  h2_relemcalphi_reltrackphi_ep->GetXaxis()->SetNdivisions(505);
  h2_relemcalphi_reltrackphi_ep->GetYaxis()->SetNdivisions(505);
  h2_relemcalphi_reltrackphi_ep->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TH2F* h2_relemcalphi_reltrackphi_em = new TH2F("h2_relemcalphi_reltrackphi_em","",16,0,2*TMath::Pi()/128.,16,0,2*TMath::Pi()/128.);
  h2_relemcalphi_reltrackphi_em->GetXaxis()->SetNoExponent(true);
  h2_relemcalphi_reltrackphi_em->GetYaxis()->SetNoExponent(true);
  h2_relemcalphi_reltrackphi_em->GetXaxis()->SetTitle("fmod(e^{-} Track projected #phi, 2#pi/128) [rad]");
  h2_relemcalphi_reltrackphi_em->GetYaxis()->SetTitle("fmod(EMCal #phi, 2#pi/128) [rad]");
  h2_relemcalphi_reltrackphi_em->GetXaxis()->SetNdivisions(505);
  h2_relemcalphi_reltrackphi_em->GetYaxis()->SetNdivisions(505);
  h2_relemcalphi_reltrackphi_em->SetTitle("#it{#bf{sPHENIX}} Internal, p+p #sqrt{s}=200 GeV, Run 53741-53783, 43.32 Million");

  TPaveText *pt15 = new TPaveText(.2, .77, .55, .92, "NDC");
  drawsPHENIXInternal(pt15);

  TLine *line_v1 = new TLine(0, 0, 0, 2*TMath::Pi()/128.);
  line_v1->SetLineColor(kRed);
  TLine *line_v2 = new TLine(2*TMath::Pi()/128., 0, 2*TMath::Pi()/128., 2*TMath::Pi()/128.);
  line_v2->SetLineColor(kRed);
  TLine *line_h1 = new TLine(0, 0, 2*TMath::Pi()/128., 0);
  line_h1->SetLineColor(kRed);
  TLine *line_h2 = new TLine(0, 2*TMath::Pi()/128., 2*TMath::Pi()/128., 2*TMath::Pi()/128.);
  line_h2->SetLineColor(kRed);

  TH1F* h1_pT_ep = new TH1F("h1_pT_ep","",100,0,3);
  h1_pT_ep->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  h1_pT_ep->GetYaxis()->SetTitle(Form("Events / %.2f [GeV/#it{c}]",3./100.));
  h1_pT_ep->SetLineColor(kRed);

  TH1F* h1_pT_em = new TH1F("h1_pT_em","",100,0,3);
  h1_pT_em->GetXaxis()->SetTitle("p_{T} [GeV/#it{c}]");
  h1_pT_em->GetYaxis()->SetTitle(Form("Events / %.2f [GeV/#it{c}]",3./100.));
  h1_pT_em->SetLineColor(kBlue);

  TLegend *legend13 = new TLegend(0.75, 0.6, 0.92, 0.7);
  legend13->AddEntry(h1_pT_ep, "e^{+}", "l");
  legend13->AddEntry(h1_pT_em, "e^{-}", "l");

  TPaveText *pt16 = new TPaveText(.55, .77, .80, .92, "NDC");
  drawsPHENIXInternal(pt16);

  TH1F* h1_Eemc_ep = new TH1F("h1_Eemc_ep","",100,0,3);
  h1_Eemc_ep->GetXaxis()->SetTitle("E [GeV]");
  h1_Eemc_ep->GetYaxis()->SetTitle(Form("Events / %.2f [GeV]",3./100.));
  h1_Eemc_ep->SetLineColor(kRed);

  TH1F* h1_Eemc_em = new TH1F("h1_Eemc_em","",100,0,3);
  h1_Eemc_em->GetXaxis()->SetTitle("E [GeV]");
  h1_Eemc_em->GetYaxis()->SetTitle(Form("Events / %.2f [GeV]",3./100.));
  h1_Eemc_em->SetLineColor(kBlue);

  TLegend *legend14 = new TLegend(0.75, 0.6, 0.92, 0.7);
  legend14->AddEntry(h1_Eemc_ep, "e^{+}", "l");
  legend14->AddEntry(h1_Eemc_em, "e^{-}", "l");

  TPaveText *pt17 = new TPaveText(.55, .77, .80, .92, "NDC");
  drawsPHENIXInternal(pt17);

  chain->Draw("_epemc_dphi>>h1_epemc_dphi");
  chain->Draw("_ememc_dphi>>h1_ememc_dphi");
  chain->Draw("_epemc_dz>>h1_epemc_dz");
  chain->Draw("_ememc_dz>>h1_ememc_dz");

  //chain->Draw("_track12_deta>>h1_track12_deta");
  //chain->Draw("_track12_dphi>>h1_track12_dphi");
  chain->Draw("_track12_deta>>h1_track12_deta","_ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1 && fabs(_track12_dphi)<0.02");
  chain->Draw("_track12_dphi>>h1_track12_dphi","_ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1 && fabs(_track12_deta)<0.02");

  //chain->Draw("_gamma_mass>>h1_gamma_mass");
  //chain->Draw("_gamma_mass>>h1_gamma_mass","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02");
  chain->Draw("_gamma_mass>>h1_gamma_mass","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _SV_radius>3 && _SV_radius<30");

  //chain->Draw("_SV_radius>>h1_SV_radius");
  //chain->Draw("_SV_radius>>h1_SV_radius","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02");
  chain->Draw("_SV_radius>>h1_SV_radius","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _gamma_mass>0 && _gamma_mass<0.1");

  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em");
  //chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02");
  //chain->Draw("_em_emcal_e/_em_p>>h1_eop_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02");
  //chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1");
  //chain->Draw("_em_emcal_e/_em_p>>h1_eop_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1");

  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep_cons","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30 && (_em_emcal_e/_em_p)>0.7 && (_em_emcal_e/_em_p)<1.3");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em_cons","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30 && (_ep_emcal_e/_ep_p)>0.7 && (_ep_emcal_e/_ep_p)<1.3");

  chain->Draw("_ep_emcal_e/_ep_p:_ep_p>>h2_eop_p_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_emcal_e/_em_p:_em_p>>h2_eop_p_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("_ep_emcal_e/_ep_p:1/_ep_pt>>h2_eop_qopt_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_emcal_e/_em_p:-1/_em_pt>>h2_eop_qopt_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("_ep_crossing>>h1_ep_crossing","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1");
  chain->Draw("_em_crossing>>h1_em_crossing","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1");

  chain->Draw("_ep_nmvtx+_ep_nintt>>h1_ep_nsi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1");
  chain->Draw("_em_nmvtx+_em_nintt>>h1_em_nsi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _SV_radius>3 && _SV_radius<30 && _gamma_mass>0 && _gamma_mass<0.1");

  chain->Draw("_ep_emcal_e/_ep_p:fmod(_ep_emcal_phi < 0 ? _ep_emcal_phi + 2*TMath::Pi() : _ep_emcal_phi, 2*TMath::Pi()/128.)>>h2_eop_relemcalphi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_emcal_e/_em_p:fmod(_em_emcal_phi < 0 ? _em_emcal_phi + 2*TMath::Pi() : _em_emcal_phi, 2*TMath::Pi()/128.)>>+h2_eop_relemcalphi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("_ep_emcal_e/_ep_p:fmod(_ep_phi_projemc < 0 ? _ep_phi_projemc + 2*TMath::Pi() : _ep_phi_projemc, 2*TMath::Pi()/128.)>>h2_eop_reltrackphi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_emcal_e/_em_p:fmod(_em_phi_projemc < 0 ? _em_phi_projemc + 2*TMath::Pi() : _em_phi_projemc, 2*TMath::Pi()/128.)>>+h2_eop_reltrackphi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("_ep_emcal_e/_ep_p:fmod(_ep_phi_projemc < 0 ? _ep_phi_projemc + 2*TMath::Pi() : _ep_phi_projemc, 2*TMath::Pi()/128.)>>h2_eop_reltrackphi_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_emcal_e/_em_p:fmod(_em_phi_projemc < 0 ? _em_phi_projemc + 2*TMath::Pi() : _em_phi_projemc, 2*TMath::Pi()/128.)>>h2_eop_reltrackphi_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("fmod(_ep_emcal_phi < 0 ? _ep_emcal_phi + 2*TMath::Pi() : _ep_emcal_phi, 2*TMath::Pi()/128.):fmod(_ep_phi_projemc < 0 ? _ep_phi_projemc + 2*TMath::Pi() : _ep_phi_projemc, 2*TMath::Pi()/128.)>>h2_relemcalphi_reltrackphi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("fmod(_em_emcal_phi < 0 ? _em_emcal_phi + 2*TMath::Pi() : _em_emcal_phi, 2*TMath::Pi()/128.):fmod(_em_phi_projemc < 0 ? _em_phi_projemc + 2*TMath::Pi() : _em_phi_projemc, 2*TMath::Pi()/128.)>>+h2_relemcalphi_reltrackphi","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("fmod(_ep_emcal_phi < 0 ? _ep_emcal_phi + 2*TMath::Pi() : _ep_emcal_phi, 2*TMath::Pi()/128.):fmod(_ep_phi_projemc < 0 ? _ep_phi_projemc + 2*TMath::Pi() : _ep_phi_projemc, 2*TMath::Pi()/128.)>>h2_relemcalphi_reltrackphi_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("fmod(_em_emcal_phi < 0 ? _em_emcal_phi + 2*TMath::Pi() : _em_emcal_phi, 2*TMath::Pi()/128.):fmod(_em_phi_projemc < 0 ? _em_phi_projemc + 2*TMath::Pi() : _em_phi_projemc, 2*TMath::Pi()/128.)>>h2_relemcalphi_reltrackphi_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("_ep_pt>>h1_pT_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_pt>>h1_pT_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

  chain->Draw("_ep_emcal_e>>h1_Eemc_ep","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");
  chain->Draw("_em_emcal_e>>h1_Eemc_em","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30","colz");

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
  arrow4l->Draw("same");
  arrow4r->Draw("same");
  can4->Update();
  can4->SaveAs("./figure/track12_dphi.pdf");

  TCanvas *can5 = new TCanvas("can5", "", 1000, 800);
  can5->cd();
  gPad->SetRightMargin(0.11);
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
  //h1_eop_em->Scale(h1_eop_ep->Integral()/h1_eop_em->Integral());
  h1_eop_ep->SetMaximum(h1_eop_ep->GetMaximum() > h1_eop_em->GetMaximum() ? 1.1*h1_eop_ep->GetMaximum() : 1.1*h1_eop_em->GetMaximum());
  h1_eop_ep->SetMinimum(0);
  h1_eop_ep->Draw("hist");
  h1_eop_em->Draw("same,hist");
  pt7->Draw("same");
  legend7->Draw("same");
  h1_eop_ep->Draw("same,hist");
  h1_eop_em->Draw("same,hist");
  can7->Update();
  can7->SaveAs("./figure/eop.pdf");
  Int_t eop1_bin = h1_eop_ep->FindBin(0.7);
  Int_t eop2_bin = h1_eop_ep->FindBin(1.3);
  std::cout<<"Integrated 0.7<E/p<1.3 region, e+ "<<h1_eop_ep->Integral(eop1_bin,eop2_bin)<<" , e- "<<h1_eop_em->Integral(eop1_bin,eop2_bin)<<std::endl;

  TCanvas *can8 = new TCanvas("can8", "", 800, 800);
  can8->cd();
  //h1_eop_em_cons->Scale(h1_eop_ep_cons->Integral()/h1_eop_em_cons->Integral());
  h1_eop_ep_cons->SetMaximum(h1_eop_ep_cons->GetMaximum() > h1_eop_em_cons->GetMaximum() ? 1.1*h1_eop_ep_cons->GetMaximum() : 1.1*h1_eop_em_cons->GetMaximum());
  h1_eop_ep_cons->SetMinimum(0);
  h1_eop_ep_cons->Draw("hist");
  h1_eop_em_cons->Draw("same,hist");
  pt8->Draw("same");
  legend8->Draw("same");
  h1_eop_ep_cons->Draw("same,hist");
  h1_eop_em_cons->Draw("same,hist");
  can8->Update();
  can8->SaveAs("./figure/eop_cons.pdf");

  TCanvas *can9 = new TCanvas("can9", "", 1200, 800);
  gPad->SetRightMargin(0.15);
  can9->cd();
  h2_eop_p_ep->Draw("colz");
  pt9->Draw("same");
  can9->Update();
  can9->SaveAs("./figure/eop_p_ep.pdf");

  TCanvas *can10 = new TCanvas("can10", "", 1200, 800);
  gPad->SetRightMargin(0.15);
  can10->cd();
  h2_eop_p_em->Draw("colz");
  pt10->Draw("same");
  can10->Update();
  can10->SaveAs("./figure/eop_p_em.pdf");

  TCanvas *can11 = new TCanvas("can11", "", 800, 800);
  can11->cd();
  gPad->SetLogy(1);
  h1_ep_crossing->SetMaximum(h1_ep_crossing->GetMaximum() > h1_ep_crossing->GetMaximum() ? 1.1*h1_ep_crossing->GetMaximum() : 1.1*h1_ep_crossing->GetMaximum());
  h1_ep_crossing->SetMinimum(0.9);
  h1_ep_crossing->Draw("hist");
  h1_em_crossing->Draw("same,hist");
  pt11->Draw("same");
  legend11->Draw("same");
  can11->Update();
  can11->SaveAs("./figure/crossing.pdf");

  TCanvas *can12 = new TCanvas("can12", "", 800, 800);
  can12->cd();
  h1_ep_nsi->SetMaximum(h1_ep_nsi->GetMaximum() > h1_ep_nsi->GetMaximum() ? 1.1*h1_ep_nsi->GetMaximum() : 1.1*h1_ep_nsi->GetMaximum());
  h1_ep_nsi->SetMinimum(0.9);
  h1_ep_nsi->Draw("hist");
  h1_em_nsi->Draw("same,hist");
  pt12->Draw("same");
  legend12->Draw("same");
  can12->Update();
  can12->SaveAs("./figure/nsi.pdf");

  TCanvas *can13 = new TCanvas("can13", "", 1200, 800);
  gPad->SetRightMargin(0.15);
  can13->cd();
  h2_eop_qopt_ep->Draw("colz");
  h2_eop_qopt_em->Draw("colz,same");
  pt13->Draw("same");
  pt14->Draw("same");
  can13->Update();
  can13->SaveAs("./figure/eop_qopt.pdf");

  //TCanvas *can14 = new TCanvas("can14", "", 1200, 800);
  //gPad->SetRightMargin(0.15);
  //can14->cd();
  //h2_eop_relemcalphi->Draw("colz");
  //pt15->Draw("same");
  //can14->Update();
  //can14->SaveAs("./figure/eop_relemcalphi.pdf");
  fit_gauss(h2_eop_relemcalphi, Form("./figure/eop_relemcalphi.pdf"), 0, pt15);
  fit_gauss(h2_eop_reltrackphi, Form("./figure/eop_reltrackphi.pdf"), 0, pt15);
  fit_gauss(h2_eop_reltrackphi_ep, Form("./figure/eop_reltrackphi_ep.pdf"), 0, pt15);
  fit_gauss(h2_eop_reltrackphi_em, Form("./figure/eop_reltrackphi_em.pdf"), 0, pt15);

  TH2* h2_relemcalphi_reltrackphi_2x2 = DuplicateTH2(h2_relemcalphi_reltrackphi);
  TCanvas *can14 = new TCanvas("can14", "", 1200, 800);
  can14->cd();
  gPad->SetRightMargin(0.15);
  //h2_relemcalphi_reltrackphi->Draw("colz");
  h2_relemcalphi_reltrackphi_2x2->Draw("colz");
  line_v1->Draw();
  line_v2->Draw();
  line_h1->Draw();
  line_h2->Draw();
  can14->Update();
  can14->SaveAs("./figure/relemcalphi_reltrackphi.pdf");

  KeepSlopeBand(h2_relemcalphi_reltrackphi_2x2, 0, 2*TMath::Pi()/128.);
  TCanvas *can14a = new TCanvas("can14a", "", 1200, 800);
  can14a->cd();
  gPad->SetRightMargin(0.15);
  h2_relemcalphi_reltrackphi_2x2->Draw("colz");
  pt15->Draw("same");
  can14a->Update();
  can14a->SaveAs("./figure/relemcalphi_reltrackphi_cutband.pdf");

  TCanvas *can15 = new TCanvas("can15", "", 1200, 800);
  can15->cd();
  h1_pT_ep->Draw("hist");
  h1_pT_em->Draw("hist,same");
  pt16->Draw("same");
  legend13->Draw("same");
  can15->Update();
  can15->SaveAs("./figure/pt.pdf");

  TCanvas *can16 = new TCanvas("can16", "", 1200, 800);
  can16->cd();
  h1_Eemc_ep->Draw("hist");
  h1_Eemc_em->Draw("hist,same");
  pt17->Draw("same");
  legend14->Draw("same");
  can16->Update();
  can16->SaveAs("./figure/Eemc.pdf");

  TH2* h2_relemcalphi_reltrackphi_ep_2x2 = DuplicateTH2(h2_relemcalphi_reltrackphi_ep);
  TCanvas *can17 = new TCanvas("can17", "", 1200, 800);
  can17->cd();
  gPad->SetRightMargin(0.15);
  //h2_relemcalphi_reltrackphi_ep->Draw("colz");
  h2_relemcalphi_reltrackphi_ep_2x2->Draw("colz");
  line_v1->Draw();
  line_v2->Draw();
  line_h1->Draw();
  line_h2->Draw();
  can17->Update();
  can17->SaveAs("./figure/relemcalphi_reltrackphi_ep.pdf");

  KeepSlopeBand(h2_relemcalphi_reltrackphi_ep_2x2, 0, 2*TMath::Pi()/128.);
  //h2_relemcalphi_reltrackphi_ep_2x2->FitSlicesY();
  //TH1 *h2_relemcalphi_reltrackphi_ep_2x2_1 = (TH1*)gDirectory->Get("h2_relemcalphi_reltrackphi_ep_2x2_1");
  //TProfile *h2_relemcalphi_reltrackphi_ep_2x2_profileX = h2_relemcalphi_reltrackphi_ep_2x2->ProfileX();
  TCanvas *can17a = new TCanvas("can17a", "", 1200, 800);
  can17a->cd();
  gPad->SetRightMargin(0.15);
  h2_relemcalphi_reltrackphi_ep_2x2->Draw("colz");
  //h2_relemcalphi_reltrackphi_ep_2x2_1->Draw("same");
  //h2_relemcalphi_reltrackphi_ep_2x2_profileX->Draw("same");
  pt15->Draw("same");
  can17a->Update();
  can17a->SaveAs("./figure/relemcalphi_reltrackphi_ep_cutband.pdf");

  TH2* h2_relemcalphi_reltrackphi_em_2x2 = DuplicateTH2(h2_relemcalphi_reltrackphi_em);
  TCanvas *can18 = new TCanvas("can18", "", 1200, 800);
  can18->cd();
  gPad->SetRightMargin(0.15);
  //h2_relemcalphi_reltrackphi_em->Draw("colz");
  h2_relemcalphi_reltrackphi_em_2x2->Draw("colz");
  line_v1->Draw();
  line_v2->Draw();
  line_h1->Draw();
  line_h2->Draw();
  can18->Update();
  can18->SaveAs("./figure/relemcalphi_reltrackphi_em.pdf");

  KeepSlopeBand(h2_relemcalphi_reltrackphi_em_2x2, 0, 2*TMath::Pi()/128.);
  //h2_relemcalphi_reltrackphi_em_2x2->FitSlicesY();
  //TH1 *h2_relemcalphi_reltrackphi_em_2x2_1 = (TH1*)gDirectory->Get("h2_relemcalphi_reltrackphi_em_2x2_1");
  //TProfile *h2_relemcalphi_reltrackphi_em_2x2_profileX = h2_relemcalphi_reltrackphi_em_2x2->ProfileX();
  TCanvas *can18a = new TCanvas("can18a", "", 1200, 800);
  can18a->cd();
  gPad->SetRightMargin(0.15);
  h2_relemcalphi_reltrackphi_em_2x2->Draw("colz");
  //h2_relemcalphi_reltrackphi_em_2x2_profileX->Draw("same");
  //h2_relemcalphi_reltrackphi_em_2x2_1->Draw("same");
  pt15->Draw("same");
  can18a->Update();
  can18a->SaveAs("./figure/relemcalphi_reltrackphi_em_cutband.pdf");

}

void drawsPHENIXInternal(TPaveText *pt)
{
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText("#it{#bf{sPHENIX}} Internal");
  pt->AddText("Au+Au #sqrt{s_{NN}}=200 GeV");
  pt->AddText(Form("Run 66641"));
  //pt->AddText("43.32 Million");
}

void drawArrow(TArrow *arrow)
{
  arrow->SetLineColor(2);
  arrow->SetLineWidth(2);
}

void fit_gauss(TH2* h, TString name, bool verbose=0, TPaveText *pt=nullptr)
{
  TGraphErrors *graph_full = new TGraphErrors();
  TGraphErrors *graph_sub = new TGraphErrors();

  int n = 1;
  for (int i = 1; i <= h->GetNbinsX(); i+=n)
  {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d fmod(phi,2#pi/128)=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    TF1 *gausFit_full = new TF1("gausFit_full", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit_full, "Q", "", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
//    double fitmin=maxBinCenter-0.5;
//    double fitmax=maxBinCenter+0.5;
    double fitmin=0.5;
    double fitmax=1.5;
    TF1 *gausFit_sub = new TF1("gausFit_sub", "gaus", fitmin, fitmax);
    projection->Fit(gausFit_sub, "Q", "", fitmin, fitmax);

    if (verbose>0)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      gPad->SetRightMargin(0.1);
      projection->Draw("hist");
      gausFit_full->SetLineColor(kRed);
      gausFit_full->Draw("same");
      gausFit_sub->SetLineColor(kBlue);
      gausFit_sub->Draw("same");
      can->Update();
      can->SaveAs(name + Form(".gausFit_%d.pdf",i));
      delete can;
    }

    Double_t mean_full = gausFit_full->GetParameter(1);
    Double_t error_full = gausFit_full->GetParError(1);

    Double_t mean_sub = gausFit_sub->GetParameter(1);
    Double_t error_sub = gausFit_sub->GetParError(1);

    if (verbose>0)
    {
      cout<<"Bin "<<i<<": full fit mean = "<<mean_full<<" +/- "<<error_full<<" , sub fit mean = "<<mean_sub<<" +/- "<<error_sub<<endl;
    }

    double center = 0;
    double width = 0;
    for (int j=0; j<n; j++)
    {
      center += h->GetXaxis()->GetBinCenter(i+j) / n;
      width += h->GetXaxis()->GetBinWidth(i+j) / 2;
    }
    graph_full->SetPoint(i-1, center, mean_full);
    graph_full->SetPointError(i-1, width, error_full);
    graph_sub->SetPoint(i-1, center, mean_sub);
    graph_sub->SetPointError(i-1, width, error_sub);
  }

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  gPad->SetRightMargin(0.15);
  h->Draw("COLZ");

  graph_full->SetMarkerStyle(20);
  graph_full->SetMarkerColor(kRed);
  //graph_full->Draw("P");
  graph_sub->SetMarkerStyle(20);
  graph_sub->SetMarkerColor(kBlack);
  graph_sub->Draw("P");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 1, h->GetXaxis()->GetXmax(), 1);
  line->SetLineColor(kRed);
  line->Draw();

  if (pt) pt->Draw("same");
  c1->Update();
  c1->SaveAs(name);
  delete c1;
}

TH2* DuplicateTH2(const TH2* originalHist)
{
    int nx = originalHist->GetNbinsX();
    int ny = originalHist->GetNbinsY();
    float xmin = originalHist->GetXaxis()->GetXmin();
    float xmax = originalHist->GetXaxis()->GetXmax();
    float ymin = originalHist->GetYaxis()->GetXmin();
    float ymax = originalHist->GetYaxis()->GetXmax();

    TH2* newHist = new TH2F(Form("%s_2x2",originalHist->GetName()), "2x2 Duplicated Histogram", 
                             2 * nx, xmin, 2 * xmax, 2 * ny, ymin, 2 * ymax);
    newHist->GetXaxis()->SetNoExponent(true);
    newHist->GetYaxis()->SetNoExponent(true);
    newHist->GetXaxis()->SetTitle(originalHist->GetXaxis()->GetTitle());
    newHist->GetYaxis()->SetTitle(originalHist->GetYaxis()->GetTitle());
    newHist->GetXaxis()->SetNdivisions(505);
    newHist->GetYaxis()->SetNdivisions(505);
    newHist->SetTitle(newHist->GetTitle());

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double content = originalHist->GetBinContent(i, j);
            double error = originalHist->GetBinError(i, j);

            newHist->SetBinContent(i, j, content);
            newHist->SetBinError(i, j, error);

            newHist->SetBinContent(i + nx, j, content);
            newHist->SetBinError(i + nx, j, error);

            newHist->SetBinContent(i, j + ny, content);
            newHist->SetBinError(i, j + ny, error);

            newHist->SetBinContent(i + nx, j + ny, content);
            newHist->SetBinError(i + nx, j + ny, error);
        }
    }

    return newHist;
}

void KeepSlopeBand(TH2* hist, double xMinAtY0 = 0.01, double xMaxAtY0 = 0.04)
{
    int nx = hist->GetNbinsX();
    int ny = hist->GetNbinsY();

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);

            double bandMin = x - xMaxAtY0;
            double bandMax = x - xMinAtY0;

            if (y < bandMin || y > bandMax) {
                hist->SetBinContent(i, j, 0);
                hist->SetBinError(i, j, 0);
            }
        }
    }
}
