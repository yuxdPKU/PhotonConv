void fit_dz()
{
  int verbosity = 1;
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  double Radius[10] = {94, 95, 96, 97, 98, 99, 100, 101, 102, 103};
  double err_Radius[10] = {0};
  double mean_ep[10], mean_em[10];
  double err_mean_ep[10], err_mean_em[10];
  double sigma_ep[10], sigma_em[10];
  double err_sigma_ep[10], err_sigma_em[10];
  for (int i=0; i<10; i++)
  {
    TFile* infile = new TFile(Form("../DetailedGeo_EMCal%d/analysis/loose_matching/eop_kfp_unlikesign.root",(int)Radius[i]),"");
    TTree* intree = (TTree*) infile->Get("tree");
    TH1F* h1_epemc_dz = new TH1F("h1_epemc_dz","",100,-10,10);
    TH1F* h1_ememc_dz = new TH1F("h1_ememc_dz","",100,-10,10);
//    intree->Draw("_epemc_dz>>h1_epemc_dz");
//    intree->Draw("_ememc_dz>>h1_ememc_dz");
    intree->Draw(Form("_epemc_dz/%g*100>>h1_epemc_dz",Radius[i]));
    intree->Draw(Form("_ememc_dz/%g*100>>h1_ememc_dz",Radius[i]));

    int maxBin_ep = h1_epemc_dz->GetMaximumBin();
    double maxBinCenter_ep = h1_epemc_dz->GetBinCenter(maxBin_ep);
    double fitmin_ep = maxBinCenter_ep - 1 * h1_epemc_dz->GetStdDev();
    double fitmax_ep = maxBinCenter_ep + 1 * h1_epemc_dz->GetStdDev();
    TF1 *gausFit_ep = new TF1("gausFit_ep", "gaus", fitmin_ep, fitmax_ep);
    h1_epemc_dz->Fit(gausFit_ep, "Q", "", fitmin_ep, fitmax_ep);

    int maxBin_em = h1_ememc_dz->GetMaximumBin();
    double maxBinCenter_em = h1_ememc_dz->GetBinCenter(maxBin_em);
    double fitmin_em = maxBinCenter_em - 1 * h1_ememc_dz->GetStdDev();
    double fitmax_em = maxBinCenter_em + 1 * h1_ememc_dz->GetStdDev();
    TF1 *gausFit_em = new TF1("gausFit_em", "gaus", fitmin_em, fitmax_em);
    h1_ememc_dz->Fit(gausFit_em, "Q", "", fitmin_em, fitmax_em);

    if (verbosity>0)
    {
      TCanvas *can = new TCanvas("can", "Canvas", 800, 600);
      can->SetTopMargin(0.12);
      h1_epemc_dz->GetXaxis()->SetTitle(Form("#DeltaZ (cm)"));
      h1_epemc_dz->SetLineColor(kRed);
      h1_epemc_dz->Draw("hist");
      h1_ememc_dz->SetLineColor(kBlue);
      h1_ememc_dz->Draw("hist,same");
      gausFit_ep->SetLineColor(kBlack);
      gausFit_ep->Draw("same");
      gausFit_em->SetLineColor(kBlack);
      gausFit_em->Draw("same");
      myText(0.2,0.95,kBlack,Form("EMCal projection R = %d cm",(int)Radius[i]));
      TLegend *legend = new TLegend(0.22, 0.65, 0.32, 0.8);
      legend->AddEntry(h1_epemc_dz, "e^{+}", "l");
      legend->AddEntry(h1_ememc_dz, "e^{-}", "l");
      legend->Draw("same");
      can->Update();
      can->SaveAs(Form("figure/dz_R%d_gausFit.pdf",(int)Radius[i]));
      delete can;
    }

    mean_ep[i] = gausFit_ep->GetParameter(1);
    sigma_ep[i] = gausFit_ep->GetParameter(2);
    err_mean_ep[i] = gausFit_ep->GetParError(1);
    err_sigma_ep[i] = gausFit_ep->GetParError(2);
    mean_em[i] = gausFit_em->GetParameter(1);
    sigma_em[i] = gausFit_em->GetParameter(2);
    err_mean_em[i] = gausFit_em->GetParError(1);
    err_sigma_em[i] = gausFit_em->GetParError(2);
  }

  TGraphErrors* gr_mean_ep = new TGraphErrors(10, Radius, mean_ep, err_Radius, err_mean_ep);
  TGraphErrors* gr_mean_em = new TGraphErrors(10, Radius, mean_em, err_Radius, err_mean_em);

  TGraphErrors* gr_sigma_ep = new TGraphErrors(10, Radius, sigma_ep, err_Radius, err_sigma_ep);
  TGraphErrors* gr_sigma_em = new TGraphErrors(10, Radius, sigma_em, err_Radius, err_sigma_em);

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  c1->SetTopMargin(0.12);
  gr_mean_ep->GetXaxis()->SetTitle("Radius (cm)");
  gr_mean_ep->GetYaxis()->SetTitle("#DeltaZ Gaus. mean (cm)");
  double ymin_mean = gr_mean_ep->GetHistogram()->GetMinimum() > gr_mean_em->GetHistogram()->GetMinimum() ? gr_mean_em->GetHistogram()->GetMinimum() : gr_mean_ep->GetHistogram()->GetMinimum();
  double ymax_mean = gr_mean_ep->GetHistogram()->GetMaximum() > gr_mean_em->GetHistogram()->GetMaximum() ? gr_mean_ep->GetHistogram()->GetMaximum() : gr_mean_em->GetHistogram()->GetMaximum();
  gr_mean_ep->SetMinimum(ymin_mean);
  gr_mean_ep->SetMaximum(ymax_mean);
  gr_mean_ep->SetMarkerStyle(20);
  gr_mean_ep->SetMarkerColor(kRed);
  gr_mean_ep->Draw("AP,same");
  gr_mean_em->SetMarkerStyle(20);
  gr_mean_em->SetMarkerColor(kBlue);
  gr_mean_em->Draw("P,same");
  c1->Update();
  c1->SaveAs("figure/dz_mean.pdf");

  TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
  c2->SetTopMargin(0.12);
  gr_sigma_ep->GetXaxis()->SetTitle("Radius (cm)");
  gr_sigma_ep->GetYaxis()->SetTitle("#DeltaZ Gaus. sigma (cm)");
  double ymin_sigma = gr_sigma_ep->GetHistogram()->GetMinimum() > gr_sigma_em->GetHistogram()->GetMinimum() ? gr_sigma_em->GetHistogram()->GetMinimum() : gr_sigma_ep->GetHistogram()->GetMinimum();
  double ymax_sigma = gr_sigma_ep->GetHistogram()->GetMaximum() > gr_sigma_em->GetHistogram()->GetMaximum() ? gr_sigma_ep->GetHistogram()->GetMaximum() : gr_sigma_em->GetHistogram()->GetMaximum();
  gr_sigma_ep->SetMinimum(ymin_sigma);
  gr_sigma_ep->SetMaximum(ymax_sigma);
  gr_sigma_ep->SetMarkerStyle(20);
  gr_sigma_ep->SetMarkerColor(kRed);
  gr_sigma_ep->Draw("AP,same");
  gr_sigma_em->SetMarkerStyle(20);
  gr_sigma_em->SetMarkerColor(kBlue);
  gr_sigma_em->Draw("P,same");
  c2->Update();
  c2->SaveAs("figure/dz_sigma.pdf");


}
