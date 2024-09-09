#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TAxis.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TNtuple.h"
#include "TPostScript.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooBreitWigner.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooArgusBG.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooKeysPdf.h"
#include "TString.h"
#include "RooPlot.h"
#include <sPhenixStyle.C>

using namespace std;
using namespace RooFit;

typedef std::vector<double> Vdouble;
typedef std::vector<int> Vint;

void fit(TTree *t_data, int runnumber);
double cal_err(double nsig, double nbkg, double esig, double ebkg, double cov, double fs, double fb);
double cal_err2(double nsig, double nbkg, double esig, double ebkg, double cov, double fs, double fb);

void fit_sw_kfp(int runnumber)
{
    SetsPhenixStyle();

    TFile* infile = new TFile(Form("eop_%d_kfp.root",runnumber));
    TTree* tree = (TTree*) infile->Get("tree");

    TCut cut = Form("_gamma_mass<0.1 && fabs(_track12_deta)<0.5");

    TFile *f_data_tmp = new TFile("./temp_data.root","recreate");
    TTree *t_data = tree->CopyTree(cut);

    fit(t_data,runnumber);

    delete t_data;
    f_data_tmp->Close();
    gSystem->Exec("rm -f temp*.root");
}

void fit(TTree *t_data, int runnumber){

	Double_t xmin=-0.5, xmax=0.5;  Double_t xbins=100;

	//data set
	int N = t_data->GetEntries();
	RooRealVar _track12_deta("_track12_deta","_track12_deta", xmin, xmax);
	RooDataSet* dataset=new RooDataSet("dataset","dataset",t_data,_track12_deta);

	//signal pdf
	RooRealVar mean("mean","mean", 0, -0.5, 0.5);
	RooRealVar sigma("sigma","sigma", 0.01, 0.0, 0.5);
	RooGaussModel sig("sig","gauss",_track12_deta, mean, sigma);

	//bkg pdf
        RooRealVar a0("a0", "constant", 0, -10, 10);
        RooRealVar a1("a1", "slope", 0, -1, 1);
        RooPolynomial bkg("bkg", "linear background", _track12_deta, RooArgList(a1, a0));

	RooRealVar nsig("nsig","sigal number",N*0.1,0,N);
	RooRealVar nbkg("nbkg","bkg number",N*0.9,0,N);

	// --- Signal + Background Pdf ---
	RooAddPdf *model = new RooAddPdf("model", "model", RooArgList(sig, bkg), RooArgList(nsig, nbkg));

	// --- Fit ---
	RooFitResult *m_fitres(0);
	m_fitres = model->fitTo(*dataset, Save(kTRUE), Extended(kTRUE));
	m_fitres->Print("v");

	TCanvas *can=new TCanvas("can","_track12_deta",800,600);
	// --- Draw ---
	RooPlot* frame = _track12_deta.frame(xmin,xmax,xbins);
	dataset->plotOn(frame, RooFit::Name("Data"));
	model->plotOn(frame, RooFit::Name("Model"), LineColor(kRed));
	model->plotOn(frame, RooFit::Name("Sig"), Components("sig"), LineStyle(kDashed), LineColor(kGreen));
	model->plotOn(frame, RooFit::Name("Bkg"), Components("bkg"), LineStyle(kDashed), LineColor(kBlue));
	dataset->plotOn(frame, RooFit::Name("Data"));

	TString a("Events / "); char b[20];  sprintf(b, "(%.3f",(xmax-xmin)/xbins); TString c(")");
	TString ytitle = a+b+c;
	frame->SetXTitle("#Delta#eta_{e^+e^-}");
	frame->SetYTitle(ytitle);

	TPaveText *pt1 = new TPaveText(.20, .72, .40, .92, "NDC");
	pt1->SetFillColor(0);
	//pt1->SetFillStyle(0);//transparent
	pt1->SetLineColor(0);
	pt1->SetBorderSize(0);
	pt1->SetTextColor(kBlack);
	pt1->AddText("#it{#bf{sPHENIX}} Internal");
	pt1->AddText("p+p #sqrt{s}=200 GeV");
	//pt1->AddText(Form("Run %d",runnumber));

	TLegend *legend = new TLegend(0.65,0.60,0.85,0.90,NULL,"brNDC");
//	legend->SetHeader("#Lambda_{c}^{+}#rightarrow#Sigma^{+}K_{S}^{0}");
	legend->AddEntry(frame->findObject("Data"),"Data","lep");
	legend->AddEntry(frame->findObject("Model"), "Fit result","l");
	legend->AddEntry(frame->findObject("Sig"), "Signal","l");
	legend->AddEntry(frame->findObject("Bkg"), "BKG","l");
        legend->SetBorderSize(0);
        legend->SetFillColor(0);
        legend->SetTextFont(22);
        legend->SetTextSize(0.04);
        legend->SetFillStyle(1001);

	can->cd();
        frame->Draw();
        pt1->Draw();
	legend->Draw();
        frame->Draw("same");

/*
	TArrow *arrow = new TArrow(sig_winlo,frame->GetMaximum(),sig_winlo,frame->GetMinimum(),0.01,"|");
	arrow->SetFillColor(4);
	arrow->SetFillStyle(1001);
	arrow->SetLineColor(4);
	arrow->SetLineWidth(1);
//	arrow->Draw();
	arrow = new TArrow(sig_winhi,frame->GetMaximum(),sig_winhi,frame->GetMinimum(),0.01,"|");
	arrow->SetFillColor(4);
	arrow->SetFillStyle(1001);
	arrow->SetLineColor(4);
	arrow->SetLineWidth(1);
//	arrow->Draw();
*/

	can->Update();
	can->Print(TString::Format("./figure/%d/fit_track12_deta.eps",runnumber));

	// ====================================
	// Now we use the SPlot class to add SWeights to our data set
	// based on our model and our yield variables

	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *dataset, model, RooArgList(nsig, nbkg));
	dataset->Print("v");
	sData->Print("v");

	Float_t track12_deta;
	double sig_sw, bkg_sw;
	TString txt_f_out = Form("eop_%d_sw.root",runnumber);
	TFile *f_out = new TFile(txt_f_out,"recreate");
	TTree *t_out = t_data->CloneTree(0);
	t_data->SetBranchAddress("_track12_deta",&track12_deta);

	t_out->Branch("sig_sw",&sig_sw,"sig_sw/D");
	t_out->Branch("bkg_sw",&bkg_sw,"bkg_sw/D");

	int nevent_sweight=0;
	int nevent_raw = t_data->GetEntries();
	for(int i=0;i<nevent_raw;i++){
		t_data->GetEntry(i);
		if(track12_deta<xmin)continue; if(track12_deta>xmax)continue;
		nevent_sweight++;
		const RooArgSet *row = dataset->get(nevent_sweight-1);
		RooRealVar *Mass_row = (RooRealVar*)row->find("_track12_deta");
		RooRealVar *mass_nsig_sw_row = (RooRealVar*)row->find("nsig_sw");
		RooRealVar *mass_nbkg_sw_row = (RooRealVar*)row->find("nbkg_sw");

		if(track12_deta!= Mass_row->getVal()) cout <<__LINE__<<"   !!!!ERROR!!!   "<< track12_deta << "  " << Mass_row->getVal() << endl;
		if(track12_deta!= Mass_row->getVal()) break;
		sig_sw = mass_nsig_sw_row->getVal();
		bkg_sw = mass_nbkg_sw_row->getVal();
		t_out->Fill();
	}
	f_out->cd();
	t_out->Write();
	f_out->Close();

	delete can;
}

double cal_err(double nsig, double nbkg, double esig, double ebkg, double cov, double fs, double fb){
	double ret = 0;
	ret += pow(-(fb*fb*nbkg/pow(fb*nbkg + fs*nsig,2)) + fb/(fb*nbkg + fs*nsig),2) * ebkg * ebkg;
	ret += pow(-(fb*fs*nbkg/pow(fb*nbkg + fs*nsig,2)),2) * esig * esig;
	ret += 2 * (-(fb*fb*nbkg/pow(fb*nbkg + fs*nsig,2)) + fb/(fb*nbkg + fs*nsig)) * (-(fb*fs*nbkg/pow(fb*nbkg + fs*nsig,2))) * cov;
	ret = sqrt(ret);
	return ret;
}

double cal_err2(double nsig, double nbkg, double esig, double ebkg, double cov, double fs, double fb){
	double ret = 0;
	ret += pow(-(fs*fs*nsig/pow(fb*nbkg + fs*nsig,2)) + fs/(fb*nbkg + fs*nsig),2) * esig * esig;
	ret += pow(-(fb*fs*nsig/pow(fb*nbkg + fs*nsig,2)),2) * ebkg * ebkg;
	ret += 2 * (-(fs*fs*nsig/pow(fb*nbkg + fs*nsig,2)) + fs/(fb*nbkg + fs*nsig)) * (-(fb*fs*nsig/pow(fb*nbkg + fs*nsig,2))) * cov;
	ret = sqrt(ret);
	return ret;
}


