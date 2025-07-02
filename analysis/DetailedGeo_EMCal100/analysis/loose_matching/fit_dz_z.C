double fit_gauss_Zslice(TH2* h, TString name, double Rnominal, int verbose=0);

void fit_dz_z()
{
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TChain* intree = new TChain("tree");
  intree->Add(Form("eop_kfp_unlikesign.root"));
  float _epemc_dz, _ememc_dz;
  float _ep_emcal_z, _em_emcal_z;
  intree->SetBranchAddress("_epemc_dz",&_epemc_dz);
  intree->SetBranchAddress("_ememc_dz",&_ememc_dz);
  intree->SetBranchAddress("_ep_emcal_z",&_ep_emcal_z);
  intree->SetBranchAddress("_em_emcal_z",&_em_emcal_z);

  TH2* h2_dz_z_ep = new TH2F("h2_dz_z_ep","e^{+} #Deltaz(track-cluster) vs. EMCal Z;EMCal Z (cm);#Deltaz(track-cluster) (cm)",60,-150,150,50,-10,10);
  TH2* h2_dz_z_em = new TH2F("h2_dz_z_em","e^{-} #Deltaz(track-cluster) vs. EMCal Z;EMCal Z (cm);#Deltaz(track-cluster) (cm)",60,-150,150,50,-10,10);

  int nevent = intree->GetEntries();
  for (int i=0; i<nevent; i++)
  {
    if (i % (Long64_t)(nevent / 10) == 0) std::cout << "progress: " << i / (Long64_t)(nevent / 10) << "0%" << std::endl;
    intree->GetEntry(i);

    h2_dz_z_ep->Fill(_ep_emcal_z,_epemc_dz);
    h2_dz_z_em->Fill(_em_emcal_z,_ememc_dz);
  }

  double Rnominal = 100;
  double slope_ep = fit_gauss_Zslice(h2_dz_z_ep, Form("figure/h2_dz_z_ep.pdf"),Rnominal,0);
  double slope_em = fit_gauss_Zslice(h2_dz_z_em, Form("figure/h2_dz_z_em.pdf"),Rnominal,0);

  cout<<"slope_ep = "<<slope_ep<<endl;
  cout<<"slope_em = "<<slope_em<<endl;

  double Rcorr_ep = Rnominal*(1-slope_ep);
  double Rcorr_em = Rnominal*(1-slope_em);
  double Rcorr_avg = (Rcorr_ep+Rcorr_em)/2;
  cout<<"Rcorr_ep = "<<Rcorr_ep<<" , Rcorr_em = "<<Rcorr_em<<" , Rcorr_avg = "<<Rcorr_avg<<endl;
}

double fit_gauss_Zslice(TH2* h, TString name, double Rnominal, int verbose=0)
{
  TGraphErrors *graph_full = new TGraphErrors();
  TGraphErrors *graph_sub = new TGraphErrors();

  int n = 0;
  for (int i = 1; i <= h->GetNbinsX(); i+=(n+1))
  {
    TH1D *projection = h->ProjectionY("projection", i, i+n);
    projection->GetYaxis()->SetTitle("Counts");
    projection->SetTitle(Form("%s , Bin %d Z=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));

    TF1 *gausFit_full = new TF1("gausFit_full", "gaus", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());
    projection->Fit(gausFit_full, "Q", "", projection->GetXaxis()->GetXmin(), projection->GetXaxis()->GetXmax());

    int maxBin = projection->GetMaximumBin();
    double maxBinCenter = projection->GetBinCenter(maxBin);
    double fitmin=maxBinCenter-3;
    double fitmax=maxBinCenter+3;
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
      myText(0.2,0.95,kBlack,Form("%s , Bin %d Z=%d cm",h->GetTitle(),i,(int)h->GetXaxis()->GetBinCenter(i)));
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

  //TF1 *linearFit = new TF1("linearFit", "[0]+[1]*x", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  TF1 *linearFit = new TF1("linearFit", "[0]+[1]*x", -100, 100);
  linearFit->SetParameter(0,-0.01);
  linearFit->SetParameter(1,0);
  linearFit->SetParLimits(0,-5,5);
  linearFit->SetParLimits(1,-1,1);
  linearFit->SetLineColor(kViolet);
  graph_sub->Fit("linearFit", "QR");

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

  double Rcorr = Rnominal*(1-slope);
  double Rcorr_err = Rnominal*slope_err;

  TPaveText *pt = new TPaveText(.05, .82, .95, .92, "NDC");
  pt->SetFillColor(0);
  //pt->SetFillStyle(0);//transparent
  pt->SetLineColor(0);
  pt->SetBorderSize(0);
  pt->SetTextColor(kBlack);
  pt->AddText(Form("Slope: (%.2f#pm%.2f)#times10^{-2}, Intercept: (%.2f#pm%.2f)#times10^{-2}",100*slope,100*slope_err,100*intercept,100*intercept_err));
  pt->AddText(Form("Proposed R_{corr} = (%.2f#pm%.2f)cm",Rcorr,Rcorr_err));
  pt->Draw("same");

  c1->Update();
  c1->SaveAs(name);
  delete c1;

  return slope;

}
