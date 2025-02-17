#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#ifndef __CINT__
#include <RooGlobalFunc.h>
#endif
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
//#include <RooDoubleCB.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooStats/SPlot.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TChain.h>
#include <TMath.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TBranch.h>

using namespace RooFit;
using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 2)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

void K_Short_fitter( bool ppSig = false )
{

  string path, inputFile, plotTitle, saveName, xAxisTitle, branch, dataType;
  double minMass = 0.300;
  double maxMass = 0.800;

  path = "/sphenix/u/xyu3/workarea/PhotonConv/ana441_trackdst_ana/analysis/";
  inputFile = path + "allkfp.root";  
  dataType = "pp";

  //branch = "K_S0_mass"; 
  branch = "gamma_mass"; 
  plotTitle = "";
  saveName = path + "K_S0"; 
  xAxisTitle = "M(#pi^{+} #pi^{-}) [GeV/#it{c}^{2}]"; 

  stringstream cutStream;
  double motherpT = 2;
  cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && sqrt(gamma_x^2 + gamma_y^2) < 3";// && sqrt(K_S0_x^2 + K_S0_y^2) < 80 " << " && K_S0_pT >= 0.0 " << " && track_1_chi2/track_1_nDoF < 5 && track_2_chi2/track_2_nDoF < 5";// && K_S0_pT >=1";
  ///cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && sqrt(K_S0_x^2 + K_S0_y^2) < 3";// && sqrt(K_S0_x^2 + K_S0_y^2) < 80 " << " && K_S0_pT >= 0.0 " << " && track_1_chi2/track_1_nDoF < 5 && track_2_chi2/track_2_nDoF < 5";// && K_S0_pT >=1";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && 2 <= K_S0_decayLength && K_S0_decayLength <= 10";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && K_S0_decayLength <= 10 && K_S0_decayLength >= 5";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && K_S0_DIRA >= 0.9";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && K_S0_DIRA <= 0.9";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && K_S0_DIRA >= 0.7 && K_S0_DIRA <= 0.9";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && K_S0_DIRA >= 0 && K_S0_DIRA <= 0.7";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && max(abs(track_1_IP_xy), abs(track_2_IP_xy)) <= 0.7";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && 0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && 0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && track_1_bunch_crossing!=0";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && 0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && track_1_bunch_crossing==0";
  //cutStream << minMass << " <= " << branch << " && " << branch << " <= " << maxMass << " && 2 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7";
  TCut masscut = cutStream.str().c_str();
 
  /*
   * Get files and data sets
   */

  TFile* dataFile = new TFile(inputFile.c_str());
  TTree* dataTree = (TTree*)dataFile->Get("DecayTree");

  //TChain* dataTree = new TChain("DecayTree");
  //dataTree->Add("/sphenix/u/xyu3/workarea/PhotonConv/KsNew/vertex/root/*.root");
  
  string datasWeight =  inputFile.substr(0, inputFile.size()-5) + "_weighted.root";
  TFile* sWeightedDataFile  = new TFile(datasWeight.c_str(), "RECREATE");
  TTree* dataSWTree = dataTree->CopyTree(masscut);
  TTree* sWeightedDataTree = dataSWTree->CloneTree(-1);

  RooRealVar mass(branch.c_str(), "mass", minMass, maxMass);      
  RooDataSet dataSet(branch.c_str(), "data", mass, Import(*sWeightedDataTree));

  /*
   * Signal Model
   */

  //calculate secondary vertex from K_S0_x and K_S0_y
  //K_S0_mass
  RooRealVar  meanK_S0("K_S0_mean", "mean", 0.517, 0.467, 0.527);
  
  RooRealVar  sigmaGaussOne("sigmaGaussOne", "sigmaGaussOne", 0.0006, 1e-4, 0.050);
  RooGaussian K_S0GaussOne("K_S0GaussOne", "K_S0GaussOne", mass, meanK_S0, sigmaGaussOne);

  //end of K_S0 Signal model


  //RooRealVar fSigGaussOne("fGaussOne", "signalD01", 0.5, 0., 1);
  RooRealVar fSig("fSig", "fSig", 1, 0., 1);

  /*
   * Background Model
   */

  RooRealVar expConst("expConst", "expConst", -1, -1e2, 1.);         //creating the decay constant
  RooExponential background("background", "background", mass, expConst);          //creating the exponential background signal
 
  /*
   * Fitting to the data
   */

  RooArgList fitModelList(K_S0GaussOne), fitFracList(fSig);
  fitModelList.add(RooArgList(background));

  RooAddPdf model("model", "model", fitModelList, fitFracList); 
  model.fitTo(dataSet);

  //Print signal yield
  string fittedSigYield = to_string_with_precision(dataSet.numEntries() * fSig.getValV(), 0);
  string fittedSigError = to_string_with_precision(dataSet.numEntries() * fSig.getError(), 0);
  int nEvents = sWeightedDataTree->GetEntries();
  cout << "*\n*\n* nEvents = " << nEvents << "\n* nSig = " << fittedSigYield << " +/- " << fittedSigError << "\n*\n*" <<endl;

  RooPlot* frame = mass.frame(Title(plotTitle.c_str()));       //creating the frame

  RooBinning bins(minMass, maxMass);
  int nBins = 25;
  bins.addUniform(nBins, minMass, maxMass);

  dataSet.plotOn(frame, Binning(bins), XErrorSize(0), DataError(RooAbsData::SumW2));          //plotting the raw unbinned data on the fram
  
  model.plotOn(frame, Components(background), LineColor(kGray), DrawOption("F"), FillColor(kGray));
  //model.plotOn(frame, Components(RooArgSet(D0GaussTwo_AuAu, background)), LineColor(kViolet), DrawOption("F"), FillColor(kViolet), MoveToBack());
  model.plotOn(frame, Components(RooArgSet(K_S0GaussOne, background)), LineColor(kAzure+8), DrawOption("F"), FillColor(kAzure+8), MoveToBack());
  
  model.plotOn(frame, LineColor(kBlack));        //plotting the fit onto the frame
  dataSet.plotOn(frame, DrawOption("PE1"), Binning(bins),XErrorSize(0), DataError(RooAbsData::SumW2));          //plotting the raw unbinned data on the fram

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooHist* pull = frame->pullHist();
  RooPlot* frame2 = mass.frame(Title(""));
  frame2->addPlotable(pull,"PE1");

  TCanvas* c = new TCanvas("massFitCanvas", "massFitCanvas",800, 600);  

  TPad mainPad("mainPad", "mainPad", 0., 0.3, 1., 1.);
  mainPad.SetBottomMargin(0);
  mainPad.Draw();
  TPad pullPad("pullPad", "pullPad", 0., 0.0, 1., 0.3);
  pullPad.SetBottomMargin(0.5);
  pullPad.SetTopMargin(0);

  pullPad.Draw();

  mainPad.cd();

  frame->SetMarkerStyle(kCircle);
  frame->SetMarkerSize(0.02);
  frame->SetLineWidth(1);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitleFont(42);
  frame->GetYaxis()->SetLabelFont(42);
  //frame->GetYaxis()->SetRangeUser(1, 100); //For 2 GeV
  //frame->GetYaxis()->SetRangeUser(1, 29); //For 4 GeV
  float binWidth = 1000.*(maxMass - minMass) / nBins;
  string yAxisTitle = "Candidates / (" + to_string_with_precision(binWidth) + " MeV)";
  frame->GetYaxis()->SetTitle(yAxisTitle.c_str());
  frame->Draw();
  frame->Print();
  
  string fittedSigYieldString = "Yield = " + to_string_with_precision(dataSet.numEntries() * fSig.getValV(), 0) + " #pm " + fittedSigError;

  string fittedK_S0_mean = to_string_with_precision(meanK_S0.getValV() * 1000, 0);
  string fittedK_S0_meanError = to_string_with_precision(meanK_S0.getError() * 1000, 0);

  string fittedK_S0_meanString = "Mean = " + fittedK_S0_mean + " #pm " + fittedK_S0_meanError + " MeV";


  string fittedsigmaGaussOne = to_string_with_precision(sigmaGaussOne.getValV() * 1000, 0);
  string fittedsigmaGaussOneError = to_string_with_precision(sigmaGaussOne.getError() * 1000, 0);

  string fittedsigmaGaussOneString = "Width = " + fittedsigmaGaussOne + " #pm " + fittedsigmaGaussOneError + " MeV";


  TPaveText *pt;
  pt = new TPaveText(0.05,0.69,0.50,0.89, "NDC"); //For 4 GeV
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  TText *pt_LaTex = pt->AddText("#it{#bf{sPHENIX}} Internal");

  //pt_LaTex = pt->AddText("#sqrt{s} = 200 GeV, ~210k events");
  //pt_LaTex = pt->AddText("461k trigger frames");
  //pt_LaTex = pt->AddText("#it{p}+#it{p}, ~55k events");
  //pt_LaTex = pt->AddText("Uncalibrated TPC-only");
  //pt_LaTex = pt->AddText("3 cm #geq r(SV)");
  //pt_LaTex = pt->AddText("Ana 418 Build");
  //pt_LaTex = pt->AddText("Track 1&2 #Chi^{2}/DoF < 5");
  //pt_LaTex = pt->AddText("1 GeV #leq K_{S}^{0} p_{T} #leq 2 GeV" );
  
  //mean, yeild, width
  pt_LaTex = pt->AddText(fittedSigYieldString.c_str());
  pt_LaTex = pt->AddText(fittedK_S0_meanString.c_str());
  pt_LaTex = pt->AddText(fittedsigmaGaussOneString.c_str());
  //pt_LaTex = pt->AddText("pp #rightarrow c#bar{c}, PYTHIA8");
  //pt_LaTex = pt->AddText("Signal Only");

  pt->SetBorderSize(0);
  pt->DrawClone();


  

  //water mark
  /*
  TPaveText *ptMark;
  ptMark = new TPaveText(0.15,0.05,0.95,0.25, "NDC"); //For 4 GeV
  ptMark->SetFillColor(0);
  ptMark->SetFillStyle(0);
  ptMark->SetTextFont(42);
  TText *ptMark_LaTex = ptMark->AddText("#it{#bf{sPHENIX}} Internal");
  ptMark_LaTex->SetTextAngle(5);
  */

  //date timestamp
   TPaveText *ptDate;
  ptDate = new TPaveText(0.67,0.92,1.05,1.00, "NDC");
  ptDate->SetFillColor(0);
  ptDate->SetFillStyle(0);
  ptDate->SetTextFont(42);
  TText *pt_LaTexDate = ptDate->AddText("10/08/2024");
  ptDate->SetBorderSize(0);
  ptDate->Draw();
  gPad->Modified();

  gPad->SetTopMargin(0.08);
  //pt_LaTex = pt->AddText("#sqrt{s} = 200 GeV");
  //pt_LaTex = pt->AddText("TPC-only, r(SV) #geq 3 cm");
  //pt_LaTex = pt->AddText("pp #rightarrow c#bar{c}, PYTHIA8");
  //pt_LaTex = pt->AddText("Signal Only");

  //pdg line
  float PDG_mass = 0.497611;

  TLine* l1 = new TLine(PDG_mass, 0, PDG_mass, 1000);
  l1->SetLineStyle(9);
  l1->SetLineColorAlpha(kRed, 1);
  l1->SetLineWidth(2);
  l1->Draw("SAME");
  

  TLegend *legend = new TLegend(0.65,0.65,0.89,0.90);
  legend->AddEntry(frame->findObject("h_K_S0_mass"),"Data","PE2");
  legend->AddEntry(l1, "K^{0}_{S} PDG Mass", "L");
  legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]"),"Fit","L");
  legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]_Comp[K_S0GaussOne,background]"),"K_{S}^{0}#rightarrow #pi^{+} #pi^{-}","f");
  //legend->AddEntry(frame->findObject("model_Norm[D0_mass]_Comp[D0GaussTwo,background]"),"mis-ID D^{0}#rightarrow K^{-} #pi^{+}","f");
  legend->AddEntry(frame->findObject("model_Norm[K_S0_mass]_Comp[background]"),"Comb. Bkg.","f");
  
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->Draw();

  //ptMark->SetBorderSize(0);
  //ptMark->DrawClone();

  gPad->Modified();
  pullPad.cd();

  frame2->SetMarkerStyle(kCircle);
  frame2->SetMarkerSize(0.02);
  frame2->SetTitle("");
  frame2->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame2->GetXaxis()->SetTitleOffset(0.9);
  frame2->GetXaxis()->SetTitleFont(42);
  frame2->GetXaxis()->SetTitleSize(0.12);
  frame2->GetXaxis()->SetLabelFont(42);
  frame2->GetXaxis()->SetLabelSize(0.12);
  frame2->GetYaxis()->SetTitle("Pull");
  frame2->GetYaxis()->SetTitleOffset(0.40);
  frame2->GetYaxis()->SetTitleFont(42);
  frame2->GetYaxis()->SetTitleSize(0.12);
  frame2->GetYaxis()->SetLabelFont(42);
  frame2->GetYaxis()->SetLabelSize(0.12);
  frame2->GetYaxis()->SetRangeUser(-6, 6);
  frame2->GetYaxis()->SetNdivisions(5);
  TF1 *plusThreeLine = new TF1("plusThreeLine",  "pol1", 0, 10000);
  plusThreeLine->SetParameters(3, 0);
  plusThreeLine->SetLineColor(1);
  plusThreeLine->SetLineStyle(2);
  plusThreeLine->SetLineWidth(2);
  TF1 *zeroLine = new TF1("zeroLine",  "pol1", 0, 10000);
  zeroLine->SetParameters(0, 0);
  zeroLine->SetLineColor(1);
  zeroLine->SetLineStyle(2);
  zeroLine->SetLineWidth(2);
  TF1 *minusThreeLine = new TF1("minusThreeLine", "pol1", 0, 10000);
  minusThreeLine->SetParameters(-3, 0);
  minusThreeLine->SetLineColor(1);
  minusThreeLine->SetLineStyle(2);
  minusThreeLine->SetLineWidth(2);
  frame2->Draw();  
  plusThreeLine->Draw("same");
  zeroLine->Draw("same");
  minusThreeLine->Draw("same");

  

  vector<string> extensions = {".C", ".pdf", ".png"};
  for (auto &extension : extensions)
  {
    string outputFile = saveName + extension;
    //string outputFile = saveName + "_" + dataType + "_ML_MLP_greater_0p99i9" + extension;
    c->SaveAs(outputFile.c_str());
  }
 
  sWeightedDataFile->Write();
  sWeightedDataFile->Close();
}
