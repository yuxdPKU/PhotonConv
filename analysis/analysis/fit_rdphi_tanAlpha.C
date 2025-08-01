std::tuple<double,double,double,double> fit_gauss_slice(TH2* h, TString name, double Rnominal, int verbose=0);

void fit_rdphi_tanAlpha()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const int n=15;
  double Radius[n] = {94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108};
  double err_Radius[n] = {0};
  double arr_slope[n];
  double arr_err_slope[n];
  double arr_intercept[n];
  double arr_err_intercept[n];

  for (int ir = 0; ir < n; ir++)
  {
  TChain* intree = new TChain("tree");
  intree->Add(Form("../DetailedGeo_EMCal%d/analysis/loose_matching/eop_kfp_unlikesign.root",(int)Radius[ir]));
  float _epemc_dphi, _ememc_dphi;
  float _ep_tanAlpha_emc, _em_tanAlpha_emc;
  intree->SetBranchAddress("_epemc_dphi",&_epemc_dphi);
  intree->SetBranchAddress("_ememc_dphi",&_ememc_dphi);
  intree->SetBranchAddress("_ep_tanAlpha_emc",&_ep_tanAlpha_emc);
  intree->SetBranchAddress("_em_tanAlpha_emc",&_em_tanAlpha_emc);

  double Rnominal = Radius[ir];
  TH2* h2_rdphi_tanAlpha = new TH2F("h2_rdphi_tanAlpha",Form("R#Delta#phi(track-cluster) vs. tan#alpha @ %.0f cm;tan#alpha;R#Delta#phi(track-cluster) (cm)",Rnominal),60,-0.6,0.6,50,-0.15*Rnominal,0.15*Rnominal);

  int nevent = intree->GetEntries();
  for (int i=0; i<nevent; i++)
  {
    if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
    intree->GetEntry(i);

    h2_rdphi_tanAlpha->Fill(_ep_tanAlpha_emc,Rnominal*_epemc_dphi);
    h2_rdphi_tanAlpha->Fill(_em_tanAlpha_emc,Rnominal*_ememc_dphi);
  }

  auto [slope, err_slope, intercept, err_intercept] = fit_gauss_slice(h2_rdphi_tanAlpha, Form("figure/h2_rdphi_tanAlpha_%d.pdf",(int)Rnominal),Rnominal,0);
  arr_slope[ir] = slope+Rnominal;
  arr_err_slope[ir] = err_slope;
  arr_intercept[ir] = intercept/Rnominal;
  arr_err_intercept[ir] = err_intercept/Rnominal;
  }

  TGraphErrors* gr_slope = new TGraphErrors(n, Radius, arr_slope, err_Radius, arr_err_slope);

  TGraphErrors* gr_intercept = new TGraphErrors(n, Radius, arr_intercept, err_Radius, arr_err_intercept);

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  c1->SetTopMargin(0.12);
  gr_slope->GetXaxis()->SetTitle("Radius (cm)");
  gr_slope->GetYaxis()->SetTitle("Corrected projection radius (cm)");
  gr_slope->SetMarkerStyle(20);
  gr_slope->SetMarkerColor(kRed);
  gr_slope->Draw("AP,same");
  c1->Update();
  c1->SaveAs("figure/RcorrFromRPhi.pdf");

  TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
  c2->SetTopMargin(0.12);
  gr_intercept->GetXaxis()->SetTitle("Radius (cm)");
  gr_intercept->GetYaxis()->SetTitle("#delta#phi (rad)");
  gr_intercept->SetMarkerStyle(20);
  gr_intercept->SetMarkerColor(kRed);
  gr_intercept->Draw("AP,same");
  c2->Update();
  c2->SaveAs("figure/dphiFromRPhi.pdf");

}

std::tuple<double,double,double,double> fit_gauss_slice(TH2* h, TString name, double Rnominal, int verbose=0)
{
  TGraphErrors *graph_full = new TGraphErrors();
  TGraphErrors *graph_sub = new TGraphErrors();

  int n = 0;
  for (int i = 1; i <= h->GetNbinsX(); i+=(n+1))
  {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d tanAlpha=%d",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    TF1 *gausFit_full = new TF1("gausFit_full", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit_full, "Q", "", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
    double fitmin=maxBinCenter-0.05;
    double fitmax=maxBinCenter+0.05;
    TF1 *gausFit_sub = new TF1("gausFit_sub", "gaus", fitmin, fitmax);
    projection->Fit(gausFit_sub, "Q", "", fitmin, fitmax);

    if (verbose>1)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      can->SetTopMargin(0.12);
      projection->Draw("hist");
      gausFit_full->SetLineColor(kRed);
      gausFit_full->Draw("same");
      gausFit_sub->SetLineColor(kBlue);
      gausFit_sub->Draw("same");
      myText(0.2,0.95,kBlack,Form("%s , Bin %d tanAlpha=%d",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));
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
    for (int j=0; j<(n+1); j++)
    {
      center += h->GetXaxis()->GetBinCenter(i+j) / (n+1);
      width += h->GetXaxis()->GetBinWidth(i+j) / 2;
    }
    graph_full->SetPoint(i-1, center, mean_full);
    graph_full->SetPointError(i-1, width, error_full);
    graph_sub->SetPoint(i-1, center, mean_sub);
    graph_sub->SetPointError(i-1, width, error_sub);
  }

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 800);
  c1->SetTopMargin(0.20);
  gPad->SetLogz(1);
  h->Draw("COLZ");

  graph_full->SetMarkerStyle(20);
  graph_full->SetMarkerColor(kRed);
  graph_full->Draw("P,same");
  graph_sub->SetMarkerStyle(20);
  graph_sub->SetMarkerColor(kBlack);
  graph_sub->Draw("P,same");

  TLine *line = new TLine(h->GetXaxis()->GetXmin(), 0, h->GetXaxis()->GetXmax(), 0);
  line->SetLineColor(kRed);
  line->Draw();

  myText(0.2,0.95,kBlack,h->GetTitle());

  TF1 *linearFit = new TF1("linearFit", "[0]+[1]*x", -0.4, 0.4);
  linearFit->SetParameter(0,3.5);
  linearFit->SetParameter(1,0);
  linearFit->SetParLimits(0,-10,10);
  linearFit->SetParLimits(1,-20,20);
  linearFit->SetLineColor(kViolet);
  graph_full->Fit("linearFit", "QR");

  double intercept = linearFit->GetParameter(0);
  double intercept_err = linearFit->GetParError(0);
  double slope = linearFit->GetParameter(1);
  double slope_err = linearFit->GetParError(1);

  double chi2 = linearFit->GetChisquare();
  int ndf = linearFit->GetNDF();
  double chi2_ndf = chi2 / ndf;

  if (verbose>0)
  {
    cout << "fit intercept = " << linearFit->GetParameter(0) << " +/- " << linearFit->GetParError(0) << endl;
    cout << "fit slope = " << linearFit->GetParameter(1) << " +/- " << linearFit->GetParError(1) << endl;
    cout<<"chi2_ndf = "<<chi2_ndf<<endl;
  }

  double Rcorr = Rnominal+slope;
  double Rcorr_err = slope_err;

  TPaveText *pt = new TPaveText(.05, .82, .95, .92, "NDC");
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText(Form("Slope: %.2f#pm%.2f, Intercept: %.2f#pm%.2f",slope,slope_err,intercept,intercept_err));
  pt->AddText(Form("Proposed R_{corr} = (%.2f#pm%.2f)cm",Rcorr,Rcorr_err));
  pt->Draw("same");

  c1->Update();
  c1->SaveAs(name);
  delete c1;

  return {slope,slope_err,intercept,intercept_err};
}
