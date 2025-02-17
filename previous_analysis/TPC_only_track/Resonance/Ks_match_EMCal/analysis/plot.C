void makecanvas1d(TH1* h1, TString name, TString xtitle, TString ytitle);
void makecanvas2d(TH2* h2, TString name, TString xtitle, TString ytitle);
void plot()
{
    TChain* chain = new TChain("DecayTree");
    chain->Add("/sphenix/u/xyu3/workarea/PhotonConv/TPC_only_track/Resonance/Ks_match_EMCal/analysis/allkfp.root");

    TH1* h1_x = new TH1F("h1_x","",100,-10,10);
    chain->Draw("K_S0_x>>h1_x");
    makecanvas1d(h1_x,Form("figure/K_S0_x.pdf"),Form("Ks PV x [cm]"),Form("Counts"));

    TH1* h1_y = new TH1F("h1_y","",100,-10,10);
    chain->Draw("K_S0_y>>h1_y");
    makecanvas1d(h1_y,Form("figure/K_S0_y.pdf"),Form("Ks PV y [cm]"),Form("Counts"));

    TH1* h1_z = new TH1F("h1_z","",100,-30,30);
    chain->Draw("K_S0_z>>h1_z");
    makecanvas1d(h1_z,Form("figure/K_S0_z.pdf"),Form("Ks PV z [cm]"),Form("Counts"));

    TH1* h1_rxy = new TH1F("h1_rxy","",100,0,100);
    chain->Draw("sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>>h1_rxy");
    makecanvas1d(h1_rxy,Form("figure/K_S0_rxy.pdf"),Form("Ks rxy(SV) [cm]"),Form("Counts"));

    TH2* h2_mass_rxy = new TH2F("h2_mass_rxy","",30,0.3,1,50,0,10);
    chain->Draw("sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y):K_S0_mass>>h2_mass_rxy");
    makecanvas2d(h2_mass_rxy,Form("figure/K_S0_mass_rxy.pdf"),Form("Ks mass [GeV]"),Form("Ks rxy(SV) [cm]"));

    TH1* h1_p = new TH1F("h1_p","",100,0.0,2.0);
    chain->Draw("K_S0_p>>h1_p","");
    makecanvas1d(h1_p,Form("figure/K_S0_p.pdf"),Form("Ks p [GeV/c]"),Form("Counts"));

    TH2* h2_mass_p = new TH2F("h2_mass_p","",30,0.3,1,50,0,2);
    chain->Draw("K_S0_p:K_S0_mass>>h2_mass_p");
    makecanvas2d(h2_mass_p,Form("figure/K_S0_mass_p.pdf"),Form("Ks mass [GeV]"),Form("Ks p [GeV]"));

    TH1* h1_pt = new TH1F("h1_pt","",100,0.0,2.0);
    chain->Draw("K_S0_pT>>h1_pt","");
    makecanvas1d(h1_pt,Form("figure/K_S0_pt.pdf"),Form("Ks pT [GeV/c]"),Form("Counts"));

    TH2* h2_mass_pt = new TH2F("h2_mass_pt","",30,0.3,1,50,0,2);
    chain->Draw("K_S0_pT:K_S0_mass>>h2_mass_pt","");
    makecanvas2d(h2_mass_pt,Form("figure/K_S0_mass_pt.pdf"),Form("Ks mass [GeV]"),Form("Ks pT [GeV/c]"));

    TH1* h1_eta = new TH1F("h1_eta","",100,-1.1,1.1);
    chain->Draw("K_S0_pseudorapidity>>h1_eta","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3");
    makecanvas1d(h1_eta,Form("figure/K_S0_pseudorapidity.pdf"),Form("Ks #eta"),Form("Counts"));

    TH2* h2_mass_eta = new TH2F("h2_mass_eta","",30,0.3,1,50,-1.1,1.1);
    chain->Draw("K_S0_pseudorapidity:K_S0_mass>>h2_mass_eta");
    makecanvas2d(h2_mass_eta,Form("figure/K_S0_mass_pseudorapidity.pdf"),Form("Ks mass [GeV]"),Form("Ks #eta"));

    TH1* h1_phi = new TH1F("h1_phi","",100,-3.15,3.15);
    chain->Draw("K_S0_phi>>h1_phi","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3");
    makecanvas1d(h1_phi,Form("figure/K_S0_phi.pdf"),Form("Ks #phi [rad]"),Form("Counts"));

    TH2* h2_mass_phi = new TH2F("h2_mass_phi","",30,0.3,1,50,-3.15,3.15);
    chain->Draw("K_S0_phi:K_S0_mass>>h2_mass_phi");
    makecanvas2d(h2_mass_phi,Form("figure/K_S0_mass_phi.pdf"),Form("Ks mass [GeV]"),Form("Ks #phi [rad]"));

    TH1* h1_dl = new TH1F("h1_dl","",100,0,20);
    chain->Draw("K_S0_decayLength>>h1_dl","");
    makecanvas1d(h1_dl,Form("figure/K_S0_dl.pdf"),Form("Ks dL [cm]"),Form("Counts"));

    TH1* h1_dldle = new TH1F("h1_dldle","",100,0,20);
    chain->Draw("K_S0_decayLength/K_S0_decayLengthErr>>h1_dldle","");
    makecanvas1d(h1_dldle,Form("figure/K_S0_dldle.pdf"),Form("Ks dL/dLErr"),Form("Counts"));

    TH2* h2_dl_rxy = new TH2F("h2_dl_rxy","",30,0,4,30,0,4);
    chain->Draw("sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y):K_S0_decayLength>>h2_dl_rxy");
    makecanvas2d(h2_dl_rxy,Form("figure/K_S0_dl_rxy.pdf"),Form("Ks dL [cm]"),Form("Ks SV rxy [cm]"));

    TH2* h2_mass_dl = new TH2F("h2_mass_dl","",30,0.3,1,50,0,20);
    chain->Draw("K_S0_decayLength:K_S0_mass>>h2_mass_dl");
    makecanvas2d(h2_mass_dl,Form("figure/K_S0_mass_dl.pdf"),Form("Ks mass [GeV]"),Form("Ks dL [cm]"));

    TH2* h2_mass_dldle = new TH2F("h2_mass_dldle","",30,0.3,1,50,0,20);
    chain->Draw("K_S0_decayLength/K_S0_decayLengthErr:K_S0_mass>>h2_mass_dldle");
    makecanvas2d(h2_mass_dldle,Form("figure/K_S0_mass_dldle.pdf"),Form("Ks mass [GeV]"),Form("Ks dL/dLErr"));

    TH1* h1_FDchi2= new TH1F("h1_FDchi2","",100,0,100);
    chain->Draw("K_S0_FDchi2>>h1_FDchi2","");
    makecanvas1d(h1_FDchi2,Form("figure/K_S0_FDchi2.pdf"),Form("Ks FDchi2"),Form("Counts"));

    TH2* h2_mass_FDchi2 = new TH2F("h2_mass_FDchi2","",30,0.3,1,50,0,100);
    chain->Draw("K_S0_FDchi2:K_S0_mass>>h2_mass_FDchi2");
    makecanvas2d(h2_mass_FDchi2,Form("figure/K_S0_mass_FDchi2.pdf"),Form("Ks mass [GeV]"),Form("Ks FDchi2"));

    TH1* h1_Ks_IP= new TH1F("h1_Ks_IP","",100,0,20);
    chain->Draw("K_S0_IP>>h1_Ks_IP","");
    makecanvas1d(h1_Ks_IP,Form("figure/K_S0_IP.pdf"),Form("Ks IP [cm]"),Form("Counts"));

    TH2* h2_mass_IP = new TH2F("h2_mass_IP","",30,0.3,1,50,0,5);
    chain->Draw("K_S0_IP:K_S0_mass>>h2_mass_IP");
    makecanvas2d(h2_mass_IP,Form("figure/K_S0_mass_IP.pdf"),Form("Ks mass [GeV]"),Form("Ks IP [cm]"));

    TH1* h1_Ks_IPxy= new TH1F("h1_Ks_IPxy","",100,0,10);
    chain->Draw("K_S0_IP_xy>>h1_Ks_IPxy","");
    makecanvas1d(h1_Ks_IPxy,Form("figure/K_S0_IPxy.pdf"),Form("Ks IPxy [cm]"),Form("Counts"));

    TH2* h2_mass_IPxy = new TH2F("h2_mass_IPxy","",30,0.3,1,50,0,10);
    chain->Draw("K_S0_IP_xy:K_S0_mass>>h2_mass_IPxy");
    makecanvas2d(h2_mass_IPxy,Form("figure/K_S0_mass_IPxy.pdf"),Form("Ks mass [GeV]"),Form("Ks IPxy [cm]"));

    TH1* h1_Ks_IPchi2= new TH1F("h1_Ks_IPchi2","",100,0,100);
    chain->Draw("K_S0_IPchi2>>h1_Ks_IPchi2","");
    makecanvas1d(h1_Ks_IPchi2,Form("figure/K_S0_IPchi2.pdf"),Form("Ks IP chi2"),Form("Counts"));

    TH2* h2_mass_IPchi2 = new TH2F("h2_mass_IPchi2","",30,0.3,1,50,0,0.5);
    chain->Draw("K_S0_IPchi2:K_S0_mass>>h2_mass_IPchi2");
    makecanvas2d(h2_mass_IPchi2,Form("figure/K_S0_mass_IPchi2.pdf"),Form("Ks mass [GeV]"),Form("Ks IPchi2"));

    TH1* h1_Ks_track1_IPxy= new TH1F("h1_Ks_track1_IPxy","",100,0,5);
    chain->Draw("track_1_IP_xy>>h1_Ks_track1_IPxy","");
    makecanvas1d(h1_Ks_track1_IPxy,Form("figure/K_S0_track1_IPxy.pdf"),Form("Ks track1 IPxy [cm]"),Form("Counts"));

    TH2* h2_mass_track1_IPxy = new TH2F("h2_mass_track1_IPxy","",30,0.3,1,50,0,5);
    chain->Draw("track_1_IP_xy:K_S0_mass>>h2_mass_track1_IPxy");
    makecanvas2d(h2_mass_track1_IPxy,Form("figure/K_S0_mass_track1_IPxy.pdf"),Form("Ks mass [GeV]"),Form("Ks track1 IPxy [cm]"));

    TH1* h1_Ks_track2_IPxy= new TH1F("h1_Ks_track2_IPxy","",100,0,5);
    chain->Draw("track_2_IP_xy>>h1_Ks_track2_IPxy","");
    makecanvas1d(h1_Ks_track2_IPxy,Form("figure/K_S0_track2_IPxy.pdf"),Form("Ks track2 IPxy [cm]"),Form("Counts"));

    TH2* h2_mass_track2_IPxy = new TH2F("h2_mass_track2_IPxy","",30,0.3,1,50,0,5);
    chain->Draw("track_2_IP_xy:K_S0_mass>>h2_mass_track2_IPxy");
    makecanvas2d(h2_mass_track2_IPxy,Form("figure/K_S0_mass_track2_IPxy.pdf"),Form("Ks mass [GeV]"),Form("Ks track2 IPxy [cm]"));

    TH1* h1_DIRA = new TH1F("h1_DIRA","",100,-1,1);
    chain->Draw("K_S0_DIRA>>h1_DIRA","");
    makecanvas1d(h1_DIRA,Form("figure/K_S0_DIRA.pdf"),Form("Ks DIRA"),Form("Counts"));

    TH2* h2_mass_DIRA = new TH2F("h2_mass_DIRA","",30,0.3,1,20,0.5,1);
    chain->Draw("K_S0_DIRA:K_S0_mass>>h2_mass_DIRA","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<4");
    makecanvas2d(h2_mass_DIRA,Form("figure/K_S0_mass_DIRA.pdf"),Form("Ks mass [GeV]"),Form("Ks DIRA"));

    TH1* h1_track1_IP= new TH1F("h1_track1_IP","",100,0,20);
    chain->Draw("track_1_IP>>h1_track1_IP","");
    makecanvas1d(h1_track1_IP,Form("figure/track1_IP.pdf"),Form("Ks track1 IP [cm]"),Form("Counts"));

    TH2* h2_mass_track1IP = new TH2F("h2_mass_track1IP","",30,0.3,1,50,0,5);
    chain->Draw("track_1_IP:K_S0_mass>>h2_mass_track1IP");
    makecanvas2d(h2_mass_track1IP,Form("figure/K_S0_mass_track1IP.pdf"),Form("Ks mass [GeV]"),Form("Ks track1 IP [cm]"));

    TH1* h1_track2_IP= new TH1F("h1_track2_IP","",100,0,20);
    chain->Draw("track_2_IP>>h1_track2_IP","");
    makecanvas1d(h1_track2_IP,Form("figure/track2_IP.pdf"),Form("Ks track2 IP [cm]"),Form("Counts"));

    TH2* h2_mass_track2IP = new TH2F("h2_mass_track2IP","",30,0.3,1,50,0,5);
    chain->Draw("track_2_IP:K_S0_mass>>h2_mass_track2IP");
    makecanvas2d(h2_mass_track2IP,Form("figure/K_S0_mass_track2IP.pdf"),Form("Ks mass [GeV]"),Form("Ks track2 IP [cm]"));

    TH1* h1_DCA = new TH1F("h1_DCA","",100,0,0.5);
    chain->Draw("track_1_track_2_DCA>>h1_DCA","");
    makecanvas1d(h1_DCA,Form("figure/K_S0_DCA.pdf"),Form("Ks track12 DCA [cm]"),Form("Counts"));

    TH2* h2_mass_DCA = new TH2F("h2_mass_DCA","",30,0.3,1,50,0,0.05);
    chain->Draw("track_1_track_2_DCA:K_S0_mass>>h2_mass_DCA");
    makecanvas2d(h2_mass_DCA,Form("figure/K_S0_mass_DCA.pdf"),Form("Ks mass [GeV]"),Form("Ks track12 DCA [cm]"));

    TH1* h1_DCAxy = new TH1F("h1_DCAxy","",100,0,0.5);
    chain->Draw("track_1_track_2_DCA_xy>>h1_DCAxy","");
    makecanvas1d(h1_DCAxy,Form("figure/K_S0_DCAxy.pdf"),Form("Ks track12 DCAxy [cm]"),Form("Counts"));

    TH2* h2_mass_DCAxy = new TH2F("h2_mass_DCAxy","",30,0.3,1,50,0,0.05);
    chain->Draw("track_1_track_2_DCA_xy:K_S0_mass>>h2_mass_DCAxy");
    makecanvas2d(h2_mass_DCAxy,Form("figure/K_S0_mass_DCAxy.pdf"),Form("Ks mass [GeV]"),Form("Ks track12 DCAxy [cm]"));

    TH1* h1_track1_crossing = new TH1F("h1_track1_crossing","",100,-100,400);
    chain->Draw("track_1_bunch_crossing>>h1_track1_crossing","");
    makecanvas1d(h1_track1_crossing,Form("figure/K_S0_track1_crossing.pdf"),Form("Ks track1 crossing"),Form("Counts"));

    TH1* h1_diff_crossing = new TH1F("h1_diff_crossing","",100,-1000,1000);
    chain->Draw("track_1_bunch_crossing-track_2_bunch_crossing>>h1_diff_crossing","");
    makecanvas1d(h1_diff_crossing,Form("figure/K_S0_diff_crossing.pdf"),Form("Ks track1 - track2 crossing"),Form("Counts"));

    //TH2* h2_rxy_DIRA = new TH2F("h2_rxy_DIRA","",20,0,150,50,0.8,1);
    TH2* h2_rxy_DIRA = new TH2F("h2_rxy_DIRA","",20,0,150,50,0,1);
    chain->Draw("K_S0_DIRA:sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>>h2_rxy_DIRA","");
    makecanvas2d(h2_rxy_DIRA,Form("figure/K_S0_rxy_DIRA.pdf"),Form("Ks rxy(SV) [cm]"),Form("Ks DIRA"));

    TH2* h2_rxy_DIRAxy = new TH2F("h2_rxy_DIRAxy","",20,0,150,50,0,1);
    chain->Draw("K_S0_DIRA_xy:sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>>h2_rxy_DIRAxy","");
    makecanvas2d(h2_rxy_DIRAxy,Form("figure/K_S0_rxy_DIRAxy.pdf"),Form("Ks rxy(SV) [cm]"),Form("Ks DIRAxy"));

    TH1* h1_track1_quality = new TH1F("h1_track1_quality","",100,0.0,200);
    chain->Draw("track_1_chi2/track_1_nDoF>>h1_track1_quality","");
    makecanvas1d(h1_track1_quality,Form("figure/K_S0_track1_quality.pdf"),Form("Ks track1 quality"),Form("Counts"));

    TH1* h1_track2_quality = new TH1F("h1_track2_quality","",100,0.0,200);
    chain->Draw("track_2_chi2/track_2_nDoF>>h1_track2_quality","");
    makecanvas1d(h1_track2_quality,Form("figure/K_S0_track2_quality.pdf"),Form("Ks track2 quality"),Form("Counts"));

    TH1* h1_track1_pt = new TH1F("h1_track1_pt","",100,0.0,2.0);
    chain->Draw("track_1_pT>>h1_track1_pt","");
    makecanvas1d(h1_track1_pt,Form("figure/K_S0_track1_pt.pdf"),Form("Ks track1 pT [GeV/c]"),Form("Counts"));

    TH1* h1_track2_pt = new TH1F("h1_track2_pt","",100,0.0,2.0);
    chain->Draw("track_2_pT>>h1_track2_pt","");
    makecanvas1d(h1_track2_pt,Form("figure/K_S0_track2_pt.pdf"),Form("Ks track2 pT [GeV/c]"),Form("Counts"));

    TH1* h1_track1_INTT_nHits = new TH1F("h1_track1_INTT_nHits","",3,0.0,3.0);
    chain->Draw("track_1_INTT_nHits>>h1_track1_INTT_nHits","");
    makecanvas1d(h1_track1_INTT_nHits,Form("figure/K_S0_track1_INTT_nHits.pdf"),Form("Ks track1 INTT nHits"),Form("Counts"));

    TH1* h1_track2_INTT_nHits = new TH1F("h1_track2_INTT_nHits","",3,0.0,3.0);
    chain->Draw("track_2_INTT_nHits>>h1_track2_INTT_nHits","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3");
    makecanvas1d(h1_track2_INTT_nHits,Form("figure/K_S0_track2_INTT_nHits.pdf"),Form("Ks track2 INTT nHits"),Form("Counts"));

    TH1* h1_track1_MVTX_nHits = new TH1F("h1_track1_MVTX_nHits","",5,0.0,5.0);
    chain->Draw("track_1_MVTX_nHits>>h1_track1_MVTX_nHits","");
    makecanvas1d(h1_track1_MVTX_nHits,Form("figure/K_S0_track1_MVTX_nHits.pdf"),Form("Ks track1 MVTX nHits"),Form("Counts"));

    TH1* h1_track2_MVTX_nHits = new TH1F("h1_track2_MVTX_nHits","",5,0.0,5.0);
    chain->Draw("track_2_MVTX_nHits>>h1_track2_MVTX_nHits","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3");
    makecanvas1d(h1_track2_MVTX_nHits,Form("figure/K_S0_track2_MVTX_nHits.pdf"),Form("Ks track2 MVTX nHits"),Form("Counts"));

    TH1* h1_track1_TPC_nHits = new TH1F("h1_track1_TPC_nHits","",50,0.0,50.0);
    chain->Draw("track_1_TPC_nHits>>h1_track1_TPC_nHits","");
    makecanvas1d(h1_track1_TPC_nHits,Form("figure/K_S0_track1_TPC_nHits.pdf"),Form("Ks track1 TPC nHits"),Form("Counts"));

    TH1* h1_track2_TPC_nHits = new TH1F("h1_track2_TPC_nHits","",50,0.0,50.0);
    chain->Draw("track_2_TPC_nHits>>h1_track2_TPC_nHits","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3");
    makecanvas1d(h1_track2_TPC_nHits,Form("figure/K_S0_track2_TPC_nHits.pdf"),Form("Ks track2 TPC nHits"),Form("Counts"));

    TH1* h1_track1_eta = new TH1F("h1_track1_eta","",100,-1.1,1.1);
    chain->Draw("track_1_pseudorapidity>>h1_track1_eta","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1");
    makecanvas1d(h1_track1_eta,Form("figure/K_S0_track1_pseudorapidity.pdf"),Form("Ks track1 #eta"),Form("Counts"));

    TH1* h1_track1_phi = new TH1F("h1_track1_phi","",100,-3.15,3.15);
    chain->Draw("track_1_phi>>h1_track1_phi","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1");
    makecanvas1d(h1_track1_phi,Form("figure/K_S0_track1_phi.pdf"),Form("Ks track1 #phi [rad]"),Form("Counts"));

    TH1* h1_track2_eta = new TH1F("h1_track2_eta","",100,-1.1,1.1);
    chain->Draw("track_2_pseudorapidity>>h1_track2_eta","");
    makecanvas1d(h1_track2_eta,Form("figure/K_S0_track2_pseudorapidity.pdf"),Form("Ks track2 #eta"),Form("Counts"));

    TH1* h1_track2_phi = new TH1F("h1_track2_phi","",100,-3.15,3.15);
    chain->Draw("track_2_phi>>h1_track2_phi","");
    makecanvas1d(h1_track2_phi,Form("figure/K_S0_track2_phi.pdf"),Form("Ks track2 #phi [rad]"),Form("Counts"));

    TH2* h2_track1_phi_eta = new TH2F("h2_track1_phi_eta","",50,-3.15,3.15,50,-1.1,1.1);
    chain->Draw("track_1_pseudorapidity:track_1_phi>>h2_track1_phi_eta","");
    makecanvas2d(h2_track1_phi_eta,Form("figure/K_S0_track1_phi_eta.pdf"),Form("Ks track1 phi [rad]"),Form("Ks track1 eta"));

    TH2* h2_track2_phi_eta = new TH2F("h2_track2_phi_eta","",50,-3.15,3.15,50,-1.1,1.1);
    chain->Draw("track_2_pseudorapidity:track_2_phi>>h2_track2_phi_eta","");
    makecanvas2d(h2_track2_phi_eta,Form("figure/K_S0_track2_phi_eta.pdf"),Form("Ks track2 phi [rad]"),Form("Ks track2 eta"));

    TH2* h2_track12_phi = new TH2F("h2_track12_phi","",50,-3.15,3.15,50,-3.15,3.15);
    chain->Draw("track_2_phi:track_1_phi>>h2_track12_phi","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1");
    makecanvas2d(h2_track12_phi,Form("figure/K_S0_track12_phi.pdf"),Form("Ks track1 phi [rad]"),Form("Ks track2 phi [rad]"));

    TH2* h2_track12_eta = new TH2F("h2_track12_eta","",50,-1.1,1.1,50,-1.1,1.1);
    chain->Draw("track_2_pseudorapidity:track_1_pseudorapidity>>h2_track12_eta","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1");
    makecanvas2d(h2_track12_eta,Form("figure/K_S0_track12_pseudorapidity.pdf"),Form("Ks track1 eta"),Form("Ks track2 eta"));

    TH2* h2_track1_mass_eta = new TH2F("h2_track1_mass_eta","",50,0.3,1.0,50,-1.1,1.1);
    chain->Draw("track_1_pseudorapidity:K_S0_mass>>h2_track1_mass_eta","");
    makecanvas2d(h2_track1_mass_eta,Form("figure/K_S0_track1_mass_eta.pdf"),Form("Ks mass [GeV]"),Form("Ks track1 eta"));

    TH2* h2_mass_DIRAxy = new TH2F("h2_mass_DIRAxy","",50,0.3,1.0,50,0.5,1);
    chain->Draw("K_S0_DIRA_xy:K_S0_mass>>h2_mass_DIRAxy","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<4");
    makecanvas2d(h2_mass_DIRAxy,Form("figure/K_S0_mass_DIRAxy.pdf"),Form("Ks mass [GeV]"),Form("Ks DIRAxy"));

    TH2* h2_DIRA_DIRAxy = new TH2F("h2_DIRA_DIRAxy","",10,0.5,1.0,10,0.5,1);
    chain->Draw("K_S0_DIRA_xy:K_S0_DIRA>>h2_DIRA_DIRAxy","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<4");
    makecanvas2d(h2_DIRA_DIRAxy,Form("figure/K_S0_DRIA_DIRAxy.pdf"),Form("Ks DIRA"),Form("Ks DIRAxy"));

    TH1* h1_ntracks = new TH1F("h1_ntracks","",100,0,100);
    chain->Draw("nEventTracks>>h1_ntracks","");
    makecanvas1d(h1_ntracks,Form("figure/ntracks.pdf"),Form("ntracks"),Form("Counts"));


    TH1* h1_mass = new TH1F("h1_mass","",50,0.3,1.0);
    chain->Draw("K_S0_mass>>h1_mass","");
    //chain->Draw("K_S0_mass>>h1_mass","(K_S0_pseudorapidity>-0.3 && K_S0_pseudorapidity<0)");
    //chain->Draw("K_S0_mass>>h1_mass","(K_S0_pseudorapidity>-0.3 && K_S0_pseudorapidity<0) && K_S0_pT<0.3");
    //chain->Draw("K_S0_mass>>h1_mass","K_S0_pT<0.3");
    //chain->Draw("K_S0_mass>>h1_mass","K_S0_pT>0.3");
    //chain->Draw("K_S0_mass>>h1_mass","track_1_bunch_crossing!=0");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<4");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<4 && K_S0_pT < 1");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1 && track_1_track_2_DCA_xy<0.005 && track_1_pT>0.5 && track_2_pT>0.5 && track_1_chi2/track_1_nDoF<100 && track_2_chi2/track_2_nDoF<100 && (K_S0_pseudorapidity>0.1 || K_S0_pseudorapidity<-0.1)");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3 && track_1_INTT_nHits>0 && track_1_MVTX_nHits>0 && track_2_INTT_nHits>0 && track_2_MVTX_nHits>0");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3 && (track_1_INTT_nHits==0 || track_1_MVTX_nHits==0 || track_2_INTT_nHits==0 || track_2_MVTX_nHits==0)");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3 && K_S0_pT<0.3");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) <= 0.7");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<3 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) <= 0.7 && K_S0_pT<0.3");
    //chain->Draw("K_S0_mass>>h1_mass","2 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7");
    //chain->Draw("K_S0_mass>>h1_mass","0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && track_1_bunch_crossing==0");
    //chain->Draw("K_S0_mass>>h1_mass","0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && track_1_bunch_crossing!=0");
    //chain->Draw("K_S0_mass>>h1_mass","0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && min(track_1_pT,track_2_pT)>0.5");
    //chain->Draw("K_S0_mass>>h1_mass","0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1 && K_S0_DIRA_xy >= 0.995 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0");
    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1 && nEventTracks<100");

    //chain->Draw("K_S0_mass>>h1_mass","sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)<8 && sqrt(K_S0_x*K_S0_x+K_S0_y*K_S0_y)>1 && K_S0_DIRA >= 0.9 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 1 && K_S0_decayLength<10 && K_S0_IP_xy>2 && K_S0_IP_xy<7 && nEventTracks<100");


    //chain->Draw("K_S0_mass>>h1_mass","0.5 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA <= 0.9 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && runNumber==52911");
    //chain->Draw("K_S0_mass>>h1_mass","2 <= K_S0_decayLength && K_S0_decayLength <= 10 && K_S0_DIRA >= 0.7 && min(abs(track_1_IP_xy), abs(track_2_IP_xy)) >= 0.7 && K_S0_pT < .3");
    makecanvas1d(h1_mass,Form("figure/K_S0_mass.pdf"),Form("Ks mass [GeV]"),Form("Counts"));
}

void makecanvas1d(TH1* h1, TString name, TString xtitle, TString ytitle="Counts")
{
    TCanvas *can = new TCanvas("can","",800,600);
    can->cd();
    h1->SetMinimum(0);
    h1->GetXaxis()->SetTitle(xtitle);
    h1->GetYaxis()->SetTitle(ytitle);
    h1->Draw();
    can->SaveAs(name);
    delete can;
}

void makecanvas2d(TH2* h2, TString name, TString xtitle, TString ytitle)
{
    TCanvas *can = new TCanvas("can","",800,600);
    can->cd();
    h2->GetXaxis()->SetTitle(xtitle);
    h2->GetYaxis()->SetTitle(ytitle);
    h2->Draw("colz");
    can->SaveAs(name);
    delete can;
}
