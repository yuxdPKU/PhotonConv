#include <filesystem>
#include "utilities.h"
#include <sPhenixStyle.C>
#include <TArrow.h>

void drawsPHENIXInternal(TPaveText *pt);
void drawArrow(TArrow *arrow);
void fit_gauss(TH2* h, TString name, bool verbose=0, TPaveText *pt=nullptr);
TH2* DuplicateTH2(const TH2* originalHist);
void KeepSlopeBand(TH2* hist, double xMinAtY0 = 0.01, double xMaxAtY0 = 0.04);

namespace fs = std::filesystem;
TVector3 z_direction(0,0,1);

void fit()
{
  SetsPhenixStyle();
  //gStyle->SetOptStat(0);

  TGaxis::SetMaxDigits(3);

  TChain* chain = new TChain("tree");
  chain->Add("eop_kfp_unlikesign.root");

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

  chain->Draw("_ep_emcal_e/_ep_p>>h1_eop_ep_cons","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30 && (_em_emcal_e/_em_p)>0.7 && (_em_emcal_e/_em_p)<1.3");
  chain->Draw("_em_emcal_e/_em_p>>h1_eop_em_cons","fabs(_track12_deta)<0.02 && fabs(_track12_dphi)<0.02 && _gamma_mass>0 && _gamma_mass<0.1 && _SV_radius>3 && _SV_radius<30 && (_ep_emcal_e/_ep_p)>0.7 && (_ep_emcal_e/_ep_p)<1.3");

  double fitmin=0.7;
  double fitmax=1.3;
  TF1 *gausFit_ep = new TF1("gausFit_ep", "gaus", 0.7, 1.3);
  TF1 *gausFit_em = new TF1("gausFit_em", "gaus", 0.7, 1.4);
  h1_eop_ep_cons->Fit(gausFit_ep, "Q", "", 0.7, 1.3);
  h1_eop_em_cons->Fit(gausFit_em, "Q", "", 0.7, 1.4);

  std::cout << "Fit to e+ eop, mean = "<<gausFit_ep->GetParameter(1)<<"+/-"<<gausFit_ep->GetParError(1)<<" , Sigma = "<<gausFit_ep->GetParameter(2)<<"+/-"<<gausFit_ep->GetParError(2)<<std::endl;
  std::cout << "Fit to e- eop, mean = "<<gausFit_em->GetParameter(1)<<"+/-"<<gausFit_em->GetParError(1)<<" , Sigma = "<<gausFit_em->GetParameter(2)<<"+/-"<<gausFit_em->GetParError(2)<<std::endl;

  TPaveText *pt8 = new TPaveText(.65, .77, .92, .92, "NDC");
  drawsPHENIXInternal(pt8);
  //pt8->AddText("Require partner E/p#in[0.7,1.3]");
  TPaveText *pt9 = new TPaveText(.20, .77, .50, .92, "NDC");
  pt9->SetFillColor(0);
  //pt9->SetFillStyle(0);//transparent
  pt9->SetLineColor(0);
  pt9->SetBorderSize(0);
  pt9->SetTextColor(kBlack);
  pt9->AddText(Form("e^{+} #mu=%.2f#pm%.2f #sigma=%.2f#pm%.2f",gausFit_ep->GetParameter(1),gausFit_ep->GetParError(1),gausFit_ep->GetParameter(2),gausFit_ep->GetParError(2)));
  pt9->AddText(Form("e^{-} #mu=%.2f#pm%.2f #sigma=%.2f#pm%.2f",gausFit_em->GetParameter(1),gausFit_em->GetParError(1),gausFit_em->GetParameter(2),gausFit_em->GetParError(2)));

  TCanvas *can = new TCanvas("can", "", 800, 800);
  can->cd();
  h1_eop_ep_cons->SetMaximum(h1_eop_ep_cons->GetMaximum() > h1_eop_em_cons->GetMaximum() ? 1.1*h1_eop_ep_cons->GetMaximum() : 1.1*h1_eop_em_cons->GetMaximum());
  h1_eop_ep_cons->SetMinimum(0);
  h1_eop_ep_cons->Draw("hist");
  h1_eop_em_cons->Draw("same,hist");
  gausFit_ep->SetLineColor(kRed+2);
  gausFit_ep->Draw("same");
  gausFit_em->SetLineColor(kBlue+2);
  gausFit_em->Draw("same");
  pt8->Draw("same");
  pt9->Draw("same");
  legend8->Draw("same");
  h1_eop_ep_cons->Draw("same,hist");
  h1_eop_em_cons->Draw("same,hist");
  can->Update();
  can->SaveAs("./figure/eop_fit.pdf");

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
  pt->AddText(Form("Run 53741-53783"));
  pt->AddText("43.32 Million");
}

void drawArrow(TArrow *arrow)
{
  arrow->SetLineColor(2);
  arrow->SetLineWidth(2);
}
