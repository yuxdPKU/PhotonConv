#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"

std::vector<float> vec_peak_pos_val;
std::vector<float> vec_peak_pos_err;

void fit_dphi(const char* filename, const char* treename, const char* branchname, TString cut, TString outfile, TString title, int runnumber) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Get the TTree from the file
    TTree *tree = (TTree*)file->Get(treename);
    if (!tree) {
        std::cerr << "Error: TTree " << treename << " not found in file" << std::endl;
        file->Close();
        return;
    }

    // Create a RooRealVar for the branch we want to fit
    RooRealVar x(branchname, branchname, 50, -0.1, 0.1);
    double binwidth = 0.2/50;

    TFile* newfile = new TFile("tmp.root","recreate");
    TTree* newtree = (TTree*) tree->CopyTree(cut);

    // Create a RooDataSet from the TTree
    RooDataSet data("data", "dataset with x", newtree, x);

    // Define the parameters for the Gaussian
    RooRealVar mean("mean", "mean of gaussian", 0.01, -0.1, 0.1);
    RooRealVar sigma("sigma", "width of gaussian", 0.01, 0.001, 0.1);

    // Create the Gaussian PDF
    RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
    RooRealVar ngauss("ngauss","number of gaussian",10,0,1000000);

    // Define the parameters for the linear background
    RooRealVar a0("a0", "constant", 0, -10, 10);
    RooRealVar a1("a1", "slope", 0, -1, 1);
    RooPolynomial poly("poly", "linear background", x, RooArgList(a1, a0));
    RooRealVar npoly("npoly","number of poly",10,0,1000000);

    // Combine the Gaussian and the polynomial into a single PDF
    RooAddPdf model("model", "gauss + poly", RooArgList(gauss, poly), RooArgList(ngauss, npoly));

    // Fit the Gaussian to the data
    RooFitResult* fit_result = model.fitTo(data, RooFit::Save());

    // Print the fit result
    fit_result->Print();

    // Plot the data and the fit result
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, RooFit::Components(gauss), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    model.plotOn(frame, RooFit::Components(poly), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

    frame->SetTitle(title);
    frame->SetXTitle("#Delta#Phi [rad]");
    frame->SetYTitle(Form("Events / (%.3g)",binwidth));
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);

    double meanValue = mean.getVal();
    double meanError = mean.getError();
    double sigmaValue = sigma.getVal();
    double sigmaError = sigma.getError();

    vec_peak_pos_val.push_back(meanValue);
    vec_peak_pos_err.push_back(meanError);

    TPaveText *pt = new TPaveText(0.20, 0.7, 0.48, 0.85, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("#mu = %.4f #pm %.4f", meanValue, meanError));
    pt->AddText(Form("#sigma = %.4f #pm %.4f", sigmaValue, sigmaError));

    // Draw the plot on a canvas
    TCanvas *canvas = new TCanvas("canvas", "fit result", 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);
    frame->Draw();
    pt->Draw();
    canvas->SaveAs(outfile);

    // Clean up
    file->Close();
    delete file;
    delete canvas;

    system("rm tmp.root");
}

void fit_dz(const char* filename, const char* treename, const char* branchname, TString cut, TString outfile, TString title, int runnumber) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Get the TTree from the file
    TTree *tree = (TTree*)file->Get(treename);
    if (!tree) {
        std::cerr << "Error: TTree " << treename << " not found in file" << std::endl;
        file->Close();
        return;
    }

    // Create a RooRealVar for the branch we want to fit
    //RooRealVar x(branchname, branchname, 50, -50, 50);
    //double binwidth = 100./50.;
    RooRealVar x(branchname, branchname, 50, -20, 20);
    double binwidth = 40./50.;

    TFile* newfile = new TFile("tmp.root","recreate");
    TTree* newtree = (TTree*) tree->CopyTree(cut);

    // Create a RooDataSet from the TTree
    RooDataSet data("data", "dataset with x", newtree, x);

    // Define the parameters for the Gaussian
    RooRealVar mean("mean", "mean of gaussian", data.mean(x), -20, 20);
    RooRealVar sigma("sigma", "width of gaussian", data.sigma(x), 0.1, 20);

    // Create the Gaussian PDF
    RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
    RooRealVar ngauss("ngauss","number of gaussian",10,1,1000000);

    // Define the parameters for the linear background
    RooRealVar a0("a0", "constant", 0, -10, 10);
    RooRealVar a1("a1", "slope", 0, -1, 1);
    RooPolynomial poly("poly", "linear background", x, RooArgList(a1, a0));
    RooRealVar npoly("npoly","number of poly",10,1,1000000);

    // Combine the Gaussian and the polynomial into a single PDF
    RooAddPdf model("model", "gauss + poly", RooArgList(gauss, poly), RooArgList(ngauss, npoly));

    // Fit the Gaussian to the data
    RooFitResult* fit_result = model.fitTo(data, RooFit::Save());

    // Print the fit result
    fit_result->Print();

    // Plot the data and the fit result
    RooPlot* frame = x.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, RooFit::Components(gauss), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    model.plotOn(frame, RooFit::Components(poly), RooFit::LineColor(kGreen), RooFit::LineStyle(kDashed));

    frame->SetTitle(title);
    frame->SetXTitle("#DeltaZ [cm]");
    frame->SetYTitle(Form("Events / (%1.0g)",binwidth));
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);

    double meanValue = mean.getVal();
    double meanError = mean.getError();
    double sigmaValue = sigma.getVal();
    double sigmaError = sigma.getError();

    vec_peak_pos_val.push_back(meanValue);
    vec_peak_pos_err.push_back(meanError);

    TPaveText *pt = new TPaveText(0.20, 0.73, 0.35, 0.88, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(Form("#mu = %.2f #pm %.2f", meanValue, meanError));
    pt->AddText(Form("#sigma = %.2f #pm %.2f", sigmaValue, sigmaError));

    // Draw the plot on a canvas
    TCanvas *canvas = new TCanvas("canvas", "fit result", 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);
    frame->Draw();
    pt->Draw();
    canvas->SaveAs(outfile);

    // Clean up
    file->Close();
    delete file;
    delete canvas;

    system("rm tmp.root");
}

void fit_dphi_cuteta(TString inputrootfile, int runnumber)
{
    std::vector<float> vec_eta_val;
    std::vector<float> vec_eta_err;

    // fit dphi in different phi_tilt region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_eta_val.clear();
    vec_eta_err.clear();
    float eta_min = -1.3;
    float eta_max = 1.3;
    float eta_range = eta_max - eta_min;
    int nstep = 26;
    float eta_stepsize = eta_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float eta_cut_min = eta_min + i * eta_stepsize;
      float eta_cut_max = eta_min + (i+1) * eta_stepsize;
      vec_eta_val.push_back((eta_cut_min+eta_cut_max)/2);
      vec_eta_err.push_back((eta_cut_max-eta_cut_min)/2);
      fit_dphi(inputrootfile, "tree", "dphi", Form("fabs(dz)<20 && eta>%f && eta<%f",eta_cut_min,eta_cut_max), Form("figure/%d/dphi_fit_cuteta_%d.pdf",runnumber,i), Form("Run %d, #eta#in[%.2f,%.2f]",runnumber,eta_cut_min,eta_cut_max),runnumber);
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_eta_val.data(), vec_peak_pos_val.data(), vec_eta_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d;#eta;#Delta#Phi Peak [rad]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", -1.5, 1.5);
    fitFunc->SetParameters(0.01, -0.5);

    fitFunc->SetParLimits(0, -0.1, 0.1);
    fitFunc->SetParLimits(1, -5, 5);

    graph->Fit(fitFunc, "");

    double p0 = fitFunc->GetParameter(0);
    double p1 = fitFunc->GetParameter(1);
    double p0_err = fitFunc->GetParError(0);
    double p1_err = fitFunc->GetParError(1);

    TString printtxt = Form("y = (%.3fx+%.3f)",p1,p0);

    std::cout << printtxt << std::endl;

    graph->Draw("AP");

    TPaveText *pt = new TPaveText(0.15, 0.15, 0.6, 0.25, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(printtxt);
    pt->Draw();

    c1->Update();
    c1->SaveAs(Form("./figure/%d/dphi_fit_eta.pdf",runnumber));

    TH2* h2_track_eta_dphi_track = new TH2F("h2_track_eta_dphi_track", "h2_track_eta_dphi_track", 50, -2, 2, 50, -0.2, 0.2);
    h2_track_eta_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
    h2_track_eta_dphi_track->GetXaxis()->SetTitle("#eta");
    h2_track_eta_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
    h2_track_eta_dphi_track->GetYaxis()->SetTitleOffset(1.2);
    h2_track_eta_dphi_track->GetZaxis()->SetTitle("Entries");
    h2_track_eta_dphi_track->GetZaxis()->SetTitleOffset(1.2);

    TFile* file = new TFile(inputrootfile,"");
    TTree* tree = (TTree*) file->Get("tree");

    float eta, dphi, dz;

    tree->SetBranchAddress("eta",&eta);
    tree->SetBranchAddress("dphi",&dphi);
    tree->SetBranchAddress("dz",&dz);

    int nevent = tree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      tree->GetEntry(i);
      if (fabs(dz)>20) continue;
      h2_track_eta_dphi_track->Fill(eta,dphi);
    }

    TCanvas *can = new TCanvas("can", "can", 800, 800);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.16);
    can->cd();
    can->SetLogz(1);
    h2_track_eta_dphi_track->Draw("COLZ");

    TGraphErrors *gr = new TGraphErrors(vec_peak_pos_val.size(), vec_eta_val.data(), vec_peak_pos_val.data(), vec_eta_err.data(), vec_peak_pos_err.data());
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->Draw("P SAME");

    can->Update();
    can->SaveAs(Form("figure/%d/dphi_eta.pdf",runnumber));

}

void fit_dphi_cutp(TString inputrootfile, int runnumber)
{
    std::vector<float> vec_p_val;
    std::vector<float> vec_p_err;

    // fit dphi in different phi_tilt region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_p_val.clear();
    vec_p_err.clear();
    float p_min = 0;
    float p_max = 10;
    float p_range = p_max - p_min;
    int nstep = 10;
    float p_stepsize = p_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float p_cut_min = p_min + i * p_stepsize;
      float p_cut_max = p_min + (i+1) * p_stepsize;
      vec_p_val.push_back((p_cut_min+p_cut_max)/2);
      vec_p_err.push_back((p_cut_max-p_cut_min)/2);
      fit_dphi(inputrootfile, "tree", "dphi", Form("fabs(dz)<20 && p>%f && p<%f",p_cut_min,p_cut_max), Form("figure/%d/dphi_fit_cutp_%d.pdf",runnumber,i), Form("Run %d, p#in[%.2f,%.2f]",runnumber,p_cut_min,p_cut_max),runnumber);
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_p_val.data(), vec_peak_pos_val.data(), vec_p_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d;p [GeV/#it{c}];#Delta#Phi Peak [rad]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", 0, 10);
    fitFunc->SetParameters(0.025, 0);

    fitFunc->SetParLimits(0, 0, 0.05);
    fitFunc->SetParLimits(1, -1, 1);

    graph->Fit(fitFunc, "");

    double p0 = fitFunc->GetParameter(0);
    double p1 = fitFunc->GetParameter(1);
    double p0_err = fitFunc->GetParError(0);
    double p1_err = fitFunc->GetParError(1);

    TString printtxt = Form("y = (%.3fx+%.3f)",p1,p0);

    std::cout << printtxt << std::endl;

    graph->Draw("AP");

    TPaveText *pt = new TPaveText(0.15, 0.15, 0.6, 0.25, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(printtxt);
    pt->Draw();

    c1->Update();
    c1->SaveAs(Form("./figure/%d/dphi_fit_p.pdf",runnumber));

    TH2* h2_track_p_dphi_track = new TH2F("h2_track_p_dphi_track", "h2_track_p_dphi_track", 50, 0, 10, 50, -0.2, 0.2);
    h2_track_p_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
    h2_track_p_dphi_track->GetXaxis()->SetTitle("p [GeV/#it{c}]");
    h2_track_p_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
    h2_track_p_dphi_track->GetYaxis()->SetTitleOffset(1.2);
    h2_track_p_dphi_track->GetZaxis()->SetTitle("Entries");
    h2_track_p_dphi_track->GetZaxis()->SetTitleOffset(1.2);

    TFile* file = new TFile(inputrootfile,"");
    TTree* tree = (TTree*) file->Get("tree");

    float p, dphi, dz;

    tree->SetBranchAddress("p",&p);
    tree->SetBranchAddress("dphi",&dphi);
    tree->SetBranchAddress("dz",&dz);

    int nevent = tree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      tree->GetEntry(i);
      if (fabs(dz)>20) continue;
      h2_track_p_dphi_track->Fill(p,dphi);
    }

    TCanvas *can = new TCanvas("can", "can", 800, 800);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.16);
    can->cd();
    can->SetLogz(1);
    h2_track_p_dphi_track->Draw("COLZ");

    TGraphErrors *gr = new TGraphErrors(vec_peak_pos_val.size(), vec_p_val.data(), vec_peak_pos_val.data(), vec_p_err.data(), vec_peak_pos_err.data());
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->Draw("P SAME");

    can->Update();
    can->SaveAs(Form("figure/%d/dphi_p.pdf",runnumber));

}

void fit_dz_cutcaloz(TString inputrootfile, int runnumber)
{
    std::vector<float> vec_calo_z_val;
    std::vector<float> vec_calo_z_err;

    // fit dz in different calo_z region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_calo_z_val.clear();
    vec_calo_z_err.clear();
    float z_min = -130;
    float z_max = 130;
    float z_range = z_max - z_min;
    int nstep = 26;
    float z_stepsize = z_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float z_cut_min = z_min + i * z_stepsize;
      float z_cut_max = z_min + (i+1) * z_stepsize;
      vec_calo_z_val.push_back((z_cut_min+z_cut_max)/2);
      vec_calo_z_err.push_back((z_cut_max-z_cut_min)/2);
      fit_dz(inputrootfile, "tree", "dz", Form("fabs(dphi)<0.1 && calo_z>%f && calo_z<%f",z_cut_min,z_cut_max), Form("figure/%d/dz_fit_cutcaloz_%d.pdf",runnumber,i), Form("Run %d, Calo Z#in[%.0f,%.0f] cm",runnumber,z_cut_min,z_cut_max),runnumber);
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_calo_z_val.data(), vec_peak_pos_val.data(), vec_calo_z_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d, TPC-EMCal matching;Calo Z [cm];#DeltaZ=Z_{track}-Z_{calo} Peak [cm]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->GetYaxis()->SetRangeUser(-50,20);

    //TF1 *fitFunc1 = new TF1("fitFunc1", "[0] + [1]*x", -130, -0.5);
    TF1 *fitFunc1 = new TF1("fitFunc1", "[0] + [1]*x", -130, -10);
    //fitFunc1->SetParameters(-5, -0.03);
    fitFunc1->SetParameters(-10, 0.1);
    fitFunc1->SetParLimits(0, -20, 20);
    fitFunc1->SetParLimits(1, -1, 1);
    graph->Fit(fitFunc1, "R");

    //TF1 *fitFunc2 = new TF1("fitFunc2", "[0] + [1]*x", 0.5, 130);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0] + [1]*x", 10, 130);
    //fitFunc2->SetParameters(5, -0.03);
    fitFunc2->SetParameters(10, 0.1);
    fitFunc2->SetParLimits(0, -20, 20);
    fitFunc2->SetParLimits(1, -1, 1);
    graph->Fit(fitFunc2, "R");

    double p0_1 = fitFunc1->GetParameter(0);
    double p1_1 = fitFunc1->GetParameter(1);
    double p0_1_err = fitFunc1->GetParError(0);
    double p1_1_err = fitFunc1->GetParError(1);
    TString printtxt1 = Form("y = %.3fx + (%.3f)",p1_1,p0_1);
    std::cout << "Fit left part: " << printtxt1 << std::endl;

    double p0_2 = fitFunc2->GetParameter(0);
    double p1_2 = fitFunc2->GetParameter(1);
    double p0_2_err = fitFunc2->GetParError(0);
    double p1_2_err = fitFunc2->GetParError(1);
    TString printtxt2 = Form("y = %.3fx + (%.3f)",p1_2,p0_2);
    std::cout << "Fit right part: " << printtxt2 << std::endl;

    graph->Draw("AP");

    fitFunc1->SetLineColor(kRed);
    fitFunc1->Draw("SAME");
    fitFunc2->SetLineColor(kBlue);
    fitFunc2->Draw("SAME");

    double dz_separation_0 = -0.23335280*2*105;
    double driftvelo_pre = (8.0 / 1000) * 107.0 / 105.0 * (1+ dz_separation_0 / 2. / 105.);
cout<<"dz_separation_0 = "<<dz_separation_0<<" , driftvelo_pre = "<<driftvelo_pre<<endl;
    double dz_separation_0_new = -(p1_1+p1_2)/2. * 2 * 105;
    double dz_separation_0_err = -(p1_1_err+p1_2_err)/2. * 2 * 105;
    double driftvelo_new = driftvelo_pre * (1+ dz_separation_0_new / 2. / 105.);
cout<<"dz_separation_0_new = "<<dz_separation_0_new<<" , driftvelo_new = "<<driftvelo_new<<endl;
    double driftvelo_err = driftvelo_pre * (fabs(dz_separation_0_err) / 2. / 105.);
    TString printtxt3 = Form("DV used in reconstruction: %.5f cm/ns",driftvelo_pre);
    TString printtxt4 = Form("Calibrated DV from data: (%.5f#pm%.5f) cm/ns",driftvelo_new,driftvelo_err);

    TPaveText *pt1 = new TPaveText(0.15, 0.30, 0.45, 0.40, "NDC");
    pt1->SetFillColor(0);
    pt1->SetTextAlign(12);
    pt1->AddText(printtxt1);
    pt1->Draw("same");

    TPaveText *pt2 = new TPaveText(0.55, 0.30, 0.85, 0.40, "NDC");
    pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
    pt2->AddText(printtxt2);
    pt2->Draw("same");

    TPaveText *pt3 = new TPaveText(0.20, 0.15, 0.90, 0.30, "NDC");
    pt3->SetFillColor(0);
    pt3->SetTextAlign(12);
    pt3->AddText(printtxt3);
    pt3->AddText(printtxt4);
    pt3->Draw("same");

    TString printtxt5;
    if (runnumber==50936)
    {
      printtxt5 = Form("2024-08-09, 111x111, TPC HV: (GEMs - 3.30 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50932)
    {
      printtxt5 = Form("2024-08-09, 111x111, TPC HV: (GEMs - 3.35 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50456)
    {
      printtxt5 = Form("2024-08-04, 56x56, TPC HV: (GEMs - 3.35 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50457)
    {
      printtxt5 = Form("2024-08-04, 56x56, TPC HV: (GEMs - 3.30 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50006)
    {
      printtxt5 = Form("2024-08-02, 28x28, TPC HV: (GEMs - 3.30 kV, CM - 43.2 kV), zero mrad crossing angle");
    }
    else if (runnumber==50015)
    {
      printtxt5 = Form("2024-08-02, 28x28, TPC HV: (GEMs - 3.35 kV, CM - 43.2 kV), zero mrad crossing angle");
    }
    else if (runnumber==49707)
    {
      printtxt5 = Form("2024-07-30, 6x6, TPC HV: (GEMs - 3.30 kV, CM - 43.3 kV), zero mrad crossing angle");
    }
    else if (runnumber==49715)
    {
      printtxt5 = Form("2024-07-30, 6x6, TPC HV: (GEMs - 3.35 kV, CM - 43.3 kV), zero mrad crossing angle");
    }
    else
    {
      printtxt5 = "";
    }

    TPaveText *pt4 = new TPaveText(0.15, 0.40, 0.90, 0.50, "NDC");
    pt4->SetFillColor(0);
    pt4->SetTextAlign(12);
    pt4->AddText(printtxt5);
    pt4->Draw("same");

    c1->Update();
    c1->SaveAs(Form("./figure/%d/dz_fit_caloz.pdf",runnumber));
}

void fit_dz_cuttrkrz(TString inputrootfile, int runnumber)
{
    std::vector<float> vec_trkr_z_val;
    std::vector<float> vec_trkr_z_err;

    // fit dz in different trkr_z region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_trkr_z_val.clear();
    vec_trkr_z_err.clear();
    float z_min = -130;
    float z_max = 130;
    float z_range = z_max - z_min;
    int nstep = 26;
    float z_stepsize = z_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float z_cut_min = z_min + i * z_stepsize;
      float z_cut_max = z_min + (i+1) * z_stepsize;
      vec_trkr_z_val.push_back((z_cut_min+z_cut_max)/2);
      vec_trkr_z_err.push_back((z_cut_max-z_cut_min)/2);
      fit_dz(inputrootfile, "tree", "dz", Form("fabs(dphi)<0.1 && track_z>%f && track_z<%f",z_cut_min,z_cut_max), Form("figure/%d/dz_fit_cuttrkrz_%d.pdf",runnumber,i), Form("Run %d, Track Z#in[%.0f,%.0f] cm",runnumber,z_cut_min,z_cut_max),runnumber);
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_trkr_z_val.data(), vec_peak_pos_val.data(), vec_trkr_z_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d, TPC-EMCal matching;Track Z [cm];#DeltaZ=Z_{track}-Z_{calo} Peak [cm]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    graph->GetYaxis()->SetRangeUser(-50,20);

    //TF1 *fitFunc1 = new TF1("fitFunc1", "[0] + [1]*x", -130, -0.5);
    TF1 *fitFunc1 = new TF1("fitFunc1", "[0] + [1]*x", -130, -10);
    //fitFunc1->SetParameters(-5, -0.03);
    fitFunc1->SetParameters(-10, 0.1);
    fitFunc1->SetParLimits(0, -20, 20);
    fitFunc1->SetParLimits(1, -1, 1);
    graph->Fit(fitFunc1, "R");

    //TF1 *fitFunc2 = new TF1("fitFunc2", "[0] + [1]*x", 0.5, 130);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0] + [1]*x", 10, 130);
    //fitFunc2->SetParameters(5, -0.03);
    fitFunc2->SetParameters(10, 0.11);
    fitFunc2->SetParLimits(0, -20, 20);
    fitFunc2->SetParLimits(1, -1, 1);
    graph->Fit(fitFunc2, "R");

    double p0_1 = fitFunc1->GetParameter(0);
    double p1_1 = fitFunc1->GetParameter(1);
    double p0_1_err = fitFunc1->GetParError(0);
    double p1_1_err = fitFunc1->GetParError(1);
    TString printtxt1 = Form("y = %.3fx + (%.3f)",p1_1,p0_1);
    std::cout << "Fit left part: " << printtxt1 << std::endl;

    double p0_2 = fitFunc2->GetParameter(0);
    double p1_2 = fitFunc2->GetParameter(1);
    double p0_2_err = fitFunc2->GetParError(0);
    double p1_2_err = fitFunc2->GetParError(1);
    TString printtxt2 = Form("y = %.3fx + (%.3f)",p1_2,p0_2);
    std::cout << "Fit right part: " << printtxt2 << std::endl;

    graph->Draw("AP");

    fitFunc1->SetLineColor(kRed);
    fitFunc1->Draw("SAME");
    fitFunc2->SetLineColor(kBlue);
    fitFunc2->Draw("SAME");

    double dz_separation_0 = -0.23335280*2*105;
    double driftvelo_pre = (8.0 / 1000) * 107.0 / 105.0 * (1+ dz_separation_0 / 2. / 105.);
    double dz_separation_0_new = -(p1_1+p1_2)/2. * 2 * 105;
    double dz_separation_0_err = -(p1_1_err+p1_2_err)/2. * 2 * 105;
    double driftvelo_new = driftvelo_pre * (1+ dz_separation_0_new / 2. / 105.);
    double driftvelo_err = driftvelo_pre * (fabs(dz_separation_0_err) / 2. / 105.);
    TString printtxt3 = Form("DV used in reconstruction: %.5f cm/ns",driftvelo_pre);
    TString printtxt4 = Form("Calibrated DV from data: (%.5f#pm%.5f) cm/ns",driftvelo_new,driftvelo_err);
    cout<<printtxt3<<endl;
    cout<<printtxt4<<endl;

    TPaveText *pt1 = new TPaveText(0.15, 0.30, 0.45, 0.40, "NDC");
    pt1->SetFillColor(0);
    pt1->SetTextAlign(12);
    pt1->AddText(printtxt1);
    pt1->Draw("same");

    TPaveText *pt2 = new TPaveText(0.55, 0.30, 0.85, 0.40, "NDC");
    pt2->SetFillColor(0);
    pt2->SetTextAlign(12);
    pt2->AddText(printtxt2);
    pt2->Draw("same");

    TPaveText *pt3 = new TPaveText(0.20, 0.15, 0.90, 0.30, "NDC");
    pt3->SetFillColor(0);
    pt3->SetTextAlign(12);
    pt3->AddText(printtxt3);
    pt3->AddText(printtxt4);
    pt3->Draw("same");

    TString printtxt5;
    if (runnumber==50936)
    {
      printtxt5 = Form("2024-08-09, 111x111, TPC HV: (GEMs - 3.30 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50932)
    {
      printtxt5 = Form("2024-08-09, 111x111, TPC HV: (GEMs - 3.35 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50456)
    {
      printtxt5 = Form("2024-08-04, 56x56, TPC HV: (GEMs - 3.35 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50457)
    {
      printtxt5 = Form("2024-08-04, 56x56, TPC HV: (GEMs - 3.30 kV, CM - 42.3 kV), 1.5 mrad crossing angle");
    }
    else if (runnumber==50006)
    {
      printtxt5 = Form("2024-08-02, 28x28, TPC HV: (GEMs - 3.30 kV, CM - 43.2 kV), zero mrad crossing angle");
    }
    else if (runnumber==50015)
    {
      printtxt5 = Form("2024-08-02, 28x28, TPC HV: (GEMs - 3.35 kV, CM - 43.2 kV), zero mrad crossing angle");
    }
    else if (runnumber==49707)
    {
      printtxt5 = Form("2024-07-30, 6x6, TPC HV: (GEMs - 3.30 kV, CM - 43.3 kV), zero mrad crossing angle");
    }
    else if (runnumber==49715)
    {
      printtxt5 = Form("2024-07-30, 6x6, TPC HV: (GEMs - 3.35 kV, CM - 43.3 kV), zero mrad crossing angle");
    }
    else
    {
      printtxt5 = "";
    }

    TPaveText *pt4 = new TPaveText(0.15, 0.40, 0.90, 0.50, "NDC");
    pt4->SetFillColor(0);
    pt4->SetTextAlign(12);
    pt4->AddText(printtxt5);
    pt4->Draw("same");

    c1->Update();
    c1->SaveAs(Form("./figure/%d/dz_fit_trkrz.pdf",runnumber));
}

void fit_dphi_cute(TString inputrootfile, int runnumber)
{
    std::vector<float> vec_e_val;
    std::vector<float> vec_e_err;

    // fit dphi in different phi_tilt region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_e_val.clear();
    vec_e_err.clear();
    float e_min = 0;
    float e_max = 10;
    float e_range = e_max - e_min;
    int nstep = 10;
    float e_stepsize = e_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float e_cut_min = e_min + i * e_stepsize;
      float e_cut_max = e_min + (i+1) * e_stepsize;
      vec_e_val.push_back((e_cut_min+e_cut_max)/2);
      vec_e_err.push_back((e_cut_max-e_cut_min)/2);
      fit_dphi(inputrootfile, "tree", "dphi", Form("fabs(dz)<20 && e>%f && e<%f",e_cut_min,e_cut_max), Form("figure/%d/dphi_fit_cute_%d.pdf",runnumber,i), Form("Run %d, E#in[%.2f,%.2f]",runnumber,e_cut_min,e_cut_max),runnumber);
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_e_val.data(), vec_peak_pos_val.data(), vec_e_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d;E [GeV];#Delta#Phi Peak [rad]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", 0, 6);
    fitFunc->SetParameters(0.025, 0);

    fitFunc->SetParLimits(0, 0, 0.05);
    fitFunc->SetParLimits(1, -1, 1);

    graph->Fit(fitFunc, "");

    double p0 = fitFunc->GetParameter(0);
    double p1 = fitFunc->GetParameter(1);
    double p0_err = fitFunc->GetParError(0);
    double p1_err = fitFunc->GetParError(1);

    TString printtxt = Form("y = (%.3fx+%.3f)",p1,p0);

    std::cout << printtxt << std::endl;

    graph->Draw("AP");

    TPaveText *pt = new TPaveText(0.15, 0.15, 0.6, 0.25, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(printtxt);
    pt->Draw();

    c1->Update();
    c1->SaveAs(Form("./figure/%d/dphi_fit_e.pdf",runnumber));

    TH2* h2_track_e_dphi_track = new TH2F("h2_track_e_dphi_track", "h2_track_e_dphi_track", 50, 0, 10, 50, -0.2, 0.2);
    h2_track_e_dphi_track->SetTitle(Form("sPHENIX Internal, Run %d",runnumber));
    h2_track_e_dphi_track->GetXaxis()->SetTitle("E [GeV]");
    h2_track_e_dphi_track->GetYaxis()->SetTitle("#Delta#Phi [rad]");
    h2_track_e_dphi_track->GetYaxis()->SetTitleOffset(1.2);
    h2_track_e_dphi_track->GetZaxis()->SetTitle("Entries");
    h2_track_e_dphi_track->GetZaxis()->SetTitleOffset(1.2);

    TFile* file = new TFile(inputrootfile,"");
    TTree* tree = (TTree*) file->Get("tree");

    float e, dphi, dz;

    tree->SetBranchAddress("e",&e);
    tree->SetBranchAddress("dphi",&dphi);
    tree->SetBranchAddress("dz",&dz);

    int nevent = tree->GetEntries();
    for (int i=0; i<nevent; i++)
    {
      tree->GetEntry(i);
      if (fabs(dz)>20) continue;
      h2_track_e_dphi_track->Fill(e,dphi);
    }

    TCanvas *can = new TCanvas("can", "can", 800, 800);
    can->SetLeftMargin(0.10);
    can->SetRightMargin(0.16);
    can->cd();
    can->SetLogz(1);
    h2_track_e_dphi_track->Draw("COLZ");

    TGraphErrors *gr = new TGraphErrors(vec_peak_pos_val.size(), vec_e_val.data(), vec_peak_pos_val.data(), vec_e_err.data(), vec_peak_pos_err.data());
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->SetLineColor(kRed);
    gr->Draw("P SAME");

    can->Update();
    can->SaveAs(Form("figure/%d/dphi_e.pdf",runnumber));
}

void fit_dphi_cutphitilt(TString inputrootfile, int runnumber){
    std::vector<float> vec_phi_tilt_val;
    std::vector<float> vec_phi_tilt_err;

    // fit dphi in different phi_tilt region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_phi_tilt_val.clear();
    vec_phi_tilt_err.clear();
    float phi_tilt_min = -0.5;
    //float phi_tilt_min = -0.4;
    //float phi_tilt_min = 0;
    float phi_tilt_max = 0.45;
    //float phi_tilt_max = 0.0;
    float phi_tilt_range = phi_tilt_max - phi_tilt_min;
    int nstep = 20;
    float phi_tilt_stepsize = phi_tilt_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float phi_tilt_cut_min = phi_tilt_min + i * phi_tilt_stepsize;
      float phi_tilt_cut_max = phi_tilt_min + (i+1) * phi_tilt_stepsize;
      vec_phi_tilt_val.push_back((phi_tilt_cut_min+phi_tilt_cut_max)/2);
      vec_phi_tilt_err.push_back((phi_tilt_cut_max-phi_tilt_cut_min)/2);
      fit_dphi(inputrootfile, "tree", "dphi", Form("fabs(dz)<20 && phi_tilt>%f && phi_tilt<%f",phi_tilt_cut_min,phi_tilt_cut_max), Form("figure/%d/dphi_fit_cutphitilt_%d.pdf",runnumber,i), Form("Run %d, phi tilt#in[%.2f,%.2f]",runnumber,phi_tilt_cut_min,phi_tilt_cut_max),runnumber);
      //fit_dphi(inputrootfile, "tree", "dphi", Form("charge==1 && fabs(dz)<20 && phi_tilt>%f && phi_tilt<%f",phi_tilt_cut_min,phi_tilt_cut_max), Form("figure/%d/dphi_fit_cutphitilt_%d.pdf",runnumber,i), Form("Run %d, phi tilt#in[%.2f,%.2f]",runnumber,phi_tilt_cut_min,phi_tilt_cut_max),runnumber);
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_phi_tilt_val.data(), vec_peak_pos_val.data(), vec_phi_tilt_err.data(), vec_peak_pos_err.data());
    graph->SetTitle(Form("Run %d;#vec{p}#bullet (#vec{R}#times#vec{z})/#vec{p}#bullet #vec{R};#Delta#Phi Peak [rad]",runnumber));
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

/*
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", -0.5, 0.5);
    fitFunc->SetParameters(0.01, -0.5);

    fitFunc->SetParLimits(0, -0.1, 0.1);
    fitFunc->SetParLimits(1, -5, 5);

    graph->Fit(fitFunc, "");

    double p0 = fitFunc->GetParameter(0);
    double p1 = fitFunc->GetParameter(1);
    double p0_err = fitFunc->GetParError(0);
    double p1_err = fitFunc->GetParError(1);

    TString printtxt = Form("y = (%.3fx+%.3f)",p1,p0);

    std::cout << printtxt << std::endl;
*/

    graph->Draw("AP");

/*
    TPaveText *pt = new TPaveText(0.15, 0.15, 0.6, 0.25, "NDC");
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->AddText(printtxt);
    pt->Draw();
*/

    c1->Update();
    c1->SaveAs(Form("./figure/%d/dphi_fit_phitilt.pdf",runnumber));
}

void fit(int runnumber) {
    gStyle->SetOptStat(0);

    //fit_dphi(Form("var_%d.root",runnumber), "tree", "dphi", Form("charge==1 && fabs(dz)<20"), Form("figure/%d/dphi_pos_fit.pdf",runnumber),Form("Run %d, Positive charge",runnumber), runnumber);
    //fit_dphi(Form("var_%d.root",runnumber), "tree", "dphi", Form("charge==-1 && fabs(dz)<20"), Form("figure/%d/dphi_neg_fit.pdf",runnumber), Form("Run %d, Negative charge",runnumber), runnumber);

    //fit_dz(Form("var_%d.root",runnumber), "tree", "dz", Form("charge==1 && fabs(dphi)<0.1"), Form("figure/%d/dz_pos_fit.pdf",runnumber), Form("Run %d, Positive charge",runnumber), runnumber);
    //fit_dz(Form("var_%d.root",runnumber), "tree", "dz", Form("charge==-1 && fabs(dphi)<0.1"), Form("figure/%d/dz_neg_fit.pdf",runnumber), Form("Run %d, Negative charge",runnumber), runnumber);

    //fit_dphi_cuteta(Form("var_%d.root",runnumber), runnumber);

    //fit_dphi_cutp(Form("var_%d.root",runnumber), runnumber);

    //fit_dphi_cute(Form("var_%d.root",runnumber), runnumber);

    //fit_dphi_cutphitilt(Form("var_%d.root",runnumber), runnumber);

    //fit_dz_cutcaloz(Form("var_%d.root",runnumber), runnumber);

    fit_dz_cuttrkrz(Form("var_%d.root",runnumber), runnumber);
}
