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
std::vector<float> vec_phi_tilt_val;
std::vector<float> vec_phi_tilt_err;

void fit_dphi(const char* filename, const char* treename, const char* branchname, TString cut, TString outfile, TString title) {
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

void fit_dz(const char* filename, const char* treename, const char* branchname, TString cut, const char* outfile, TString title) {
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
    RooRealVar x(branchname, branchname, 50, -50, 50);
    double binwidth = 100/50;

    TFile* newfile = new TFile("tmp.root","recreate");
    TTree* newtree = (TTree*) tree->CopyTree(cut);

    // Create a RooDataSet from the TTree
    RooDataSet data("data", "dataset with x", newtree, x);

    // Define the parameters for the Gaussian
    RooRealVar mean("mean", "mean of gaussian", 0, -20, 20);
    RooRealVar sigma("sigma", "width of gaussian", 5, 0.1, 20);

    // Create the Gaussian PDF
    RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
    RooRealVar ngauss("ngauss","number of gaussian",6000,100,1000000);

    // Define the parameters for the linear background
    RooRealVar a0("a0", "constant", 0, -10, 10);
    RooRealVar a1("a1", "slope", 0, -1, 1);
    RooPolynomial poly("poly", "linear background", x, RooArgList(a1, a0));
    RooRealVar npoly("npoly","number of poly",3000,100,1000000);

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

    TPaveText *pt = new TPaveText(0.20, 0.7, 0.48, 0.85, "NDC");
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

int fit() {
    //fit_dphi("var.root", "tree", "dphi", Form("charge==1 && fabs(dz)<20"), "figure/dphi_pos_fit.pdf",Form("Run 46730, Positive charge"));
    //fit_dphi("var.root", "tree", "dphi", Form("charge==-1 && fabs(dz)<20"), "figure/dphi_neg_fit.pdf", Form("Run 46730, Negative charge"));

    //fit_dz("var.root", "tree", "dz", Form("charge==1 && fabs(dphi)<0.1"), "figure/dz_pos_fit.pdf", Form("Run 46730, Positive charge"));
    //fit_dz("var.root", "tree", "dz", Form("charge==-1 && fabs(dphi)<0.1"), "figure/dz_neg_fit.pdf", Form("Run 46730, Negative charge"));

    // fit dphi in different phi_tilt region
    vec_peak_pos_val.clear();
    vec_peak_pos_err.clear();
    vec_phi_tilt_val.clear();
    vec_phi_tilt_err.clear();
    float phi_tilt_min = -0.5;
    float phi_tilt_max = 0.45;
    float phi_tilt_range = phi_tilt_max - phi_tilt_min;
    int nstep = 20;
    float phi_tilt_stepsize = phi_tilt_range / nstep;
    for (int i=0; i<nstep; i++)
    {
      float phi_tilt_cut_min = phi_tilt_min + i * phi_tilt_stepsize;
      float phi_tilt_cut_max = phi_tilt_min + (i+1) * phi_tilt_stepsize;
      vec_phi_tilt_val.push_back((phi_tilt_cut_min+phi_tilt_cut_max)/2);
      vec_phi_tilt_err.push_back((phi_tilt_cut_max-phi_tilt_cut_min)/2);
      fit_dphi("var.root", "tree", "dphi", Form("fabs(dz)<20 && phi_tilt>%f && phi_tilt<%f",phi_tilt_cut_min,phi_tilt_cut_max), Form("figure/dphi_fit_cutphitilt_%d.pdf",i), Form("Run 46730, phi tilt#in[%.2f,%.2f]",phi_tilt_cut_min,phi_tilt_cut_max));
    }

    TCanvas *c1 = new TCanvas("c1", "Fitting Example", 800, 600);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    TGraphErrors *graph = new TGraphErrors(vec_peak_pos_val.size(), vec_phi_tilt_val.data(), vec_peak_pos_val.data(), vec_phi_tilt_err.data(), vec_peak_pos_err.data());
    graph->SetTitle("Run 46730;#vec{p}#bullet (#vec{R}#times#vec{z})/#vec{p}#bullet #vec{R};#Delta#Phi Peak [cm]");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x", 0, 6);
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
    c1->SaveAs("./figure/dphi_fit_phitilt.pdf");

    return 0;
}

