#include <filesystem>
#include "../utilities.h"
#include <sPhenixStyle.C>
#include <TArrow.h>

void drawsPHENIXInternal(TPaveText *pt);
float GetMeanpT(TChain* chain, TCut cut);

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void plot()
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TGaxis::SetMaxDigits(3);

  TChain* chain = new TChain("tree");
  chain->Add("../eop_kfp_unlikesign.root");
  //chain->Add("eop_kfp_likesign.root");

  //date timestamp
  TPaveText *ptDate;
  ptDate = new TPaveText(0.67,0.92,1.05,1.00, "NDC");
  ptDate->SetFillColor(0);
  ptDate->SetFillStyle(0);
  ptDate->SetTextFont(42);
  ptDate->SetTextSize(0.05);
  TText *pt_LaTexDate = ptDate->AddText("11/03/2025");
  ptDate->SetBorderSize(0);

  TH1F* h1_gamma_mass = new TH1F("h1_gamma_mass","",50,0,50); // unit MeV
  h1_gamma_mass->GetXaxis()->SetTitle("#it{M}_{e^{+}e^{-}} [MeV]");
  h1_gamma_mass->GetYaxis()->SetTitle(Form("Events / 1 MeV"));
  h1_gamma_mass->SetLineColor(kBlack);

  TPaveText *pt5 = new TPaveText(.40, .57, .85, .77, "NDC");
  drawsPHENIXInternal(pt5);

  TH1F* h1_SV_radius = new TH1F("h1_SV_radius","",30,0,30);
  h1_SV_radius->GetXaxis()->SetTitle("#it{R}_{e^{+}e^{-}} [cm]");
  h1_SV_radius->GetYaxis()->SetTitle(Form("Events / 1 cm"));
  h1_SV_radius->SetLineColor(kBlack);

  TPaveText *pt6 = new TPaveText(.47, .67, .92, .87, "NDC");
  drawsPHENIXInternal(pt6);

  TH1F* h1_eop_ep_cons = new TH1F("h1_eop_ep_cons","",40,0,2);
  h1_eop_ep_cons->GetXaxis()->SetTitle("E/p");
  h1_eop_ep_cons->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_ep_cons->SetLineColor(kRed);

  TH1F* h1_eop_em_cons = new TH1F("h1_eop_em_cons","",40,0,2);
  h1_eop_em_cons->GetXaxis()->SetTitle("E/p");
  h1_eop_em_cons->GetYaxis()->SetTitle(Form("Events / %.2f",2./40.));
  h1_eop_em_cons->SetLineColor(kBlue);

  TLegend *legend8 = new TLegend(0.80, 0.5, 0.92, 0.6);
  legend8->AddEntry(h1_eop_ep_cons, "e^{+}", "l");
  legend8->AddEntry(h1_eop_em_cons, "e^{-}", "l");

  TPaveText *pt8 = new TPaveText(.62, .67, .92, .87, "NDC");
  drawsPHENIXInternal(pt8);
  float meanpt = GetMeanpT(chain,Form("fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30"));
  std::cout<<"Mean pT = "<<meanpt<<std::endl;
  pt8->AddText(Form("Mean pT ~ 1 GeV"));

  chain->Draw("1000*_gamma_mass>>h1_gamma_mass","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _SV_radius>3 && _SV_radius<30");

  chain->Draw("_SV_radius>>h1_SV_radius","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _gamma_mass>0 && _gamma_mass<0.1");

  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep_cons","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30 && (_em_emcal_e/_em_p)>0.7 && (_em_emcal_e/_em_p)<1.3");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em_cons","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30 && (_ep_emcal_e/_ep_p)>0.7 && (_ep_emcal_e/_ep_p)<1.3");

  TCanvas *can5 = new TCanvas("can5", "", 800, 800);
  gPad->SetTopMargin(0.08);
  can5->cd();
  h1_gamma_mass->SetMinimum(0);
  h1_gamma_mass->Draw("hist");
  pt5->Draw("same");
  ptDate->Draw();
  gPad->Modified();
  can5->Update();
  can5->SaveAs("./figure/gamma_mass.pdf");

  TCanvas *can6 = new TCanvas("can6", "", 800, 800);
  gPad->SetTopMargin(0.08);
  can6->cd();
  h1_SV_radius->SetMinimum(0);
  h1_SV_radius->Draw("hist");
  pt6->Draw("same");
  ptDate->Draw();
  gPad->Modified();
  can6->Update();
  can6->SaveAs("./figure/gamma_radius.pdf");

  TCanvas *can8 = new TCanvas("can8", "", 800, 800);
  gPad->SetTopMargin(0.08);
  can8->cd();
  h1_eop_ep_cons->SetMaximum(h1_eop_ep_cons->GetMaximum() > h1_eop_em_cons->GetMaximum() ? 1.1*h1_eop_ep_cons->GetMaximum() : 1.1*h1_eop_em_cons->GetMaximum());
  h1_eop_ep_cons->SetMinimum(0);
  h1_eop_ep_cons->Draw("hist");
  h1_eop_em_cons->Draw("same,hist");
  pt8->Draw("same");
  legend8->Draw("same");
  h1_eop_ep_cons->Draw("same,hist");
  h1_eop_em_cons->Draw("same,hist");
  ptDate->Draw();
  gPad->Modified();
  can8->Update();
  can8->SaveAs("./figure/eop_cons.pdf");
}

void drawsPHENIXInternal(TPaveText *pt)
{
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText("#it{#bf{sPHENIX}} Preliminary");
  pt->AddText("p+p #sqrt{s}=200 GeV");
  pt->AddText(Form("Conversion Photon Candidates"));
}

float GetMeanpT(TChain* chain, TCut cut)
{
  TTree* tree = chain->CopyTree(cut);
  float _em_emcal_e, _em_p, _em_pt;
  float _ep_emcal_e, _ep_p, _ep_pt;
  tree->SetBranchAddress("_em_emcal_e",&_em_emcal_e);
  tree->SetBranchAddress("_em_p",&_em_p);
  tree->SetBranchAddress("_em_pt",&_em_pt);
  tree->SetBranchAddress("_ep_emcal_e",&_ep_emcal_e);
  tree->SetBranchAddress("_ep_p",&_ep_p);
  tree->SetBranchAddress("_ep_pt",&_ep_pt);

  float ntot=0;
  float nsum=0;
  int nevent=tree->GetEntries();
  for (int i=0; i<nevent; i++)
  {
    tree->GetEntry(i);
    if (_em_emcal_e/_em_p>0.7 && _em_emcal_e/_em_p<1.3 && _ep_emcal_e/_ep_p<2)
    {
      ntot++;
      nsum+=_ep_pt;
    }
    if (_ep_emcal_e/_ep_p>0.7 && _ep_emcal_e/_ep_p<1.3 && _em_emcal_e/_em_p<2)
    {
      ntot++;
      nsum+=_em_pt;
    }
  }
  float meanpt=nsum/ntot;

  return meanpt;
}
