void combine(TString infile, TString branchname, TString outfile, TCut cut)
{
  std::cout<<"Processing "<<infile<<" branch "<<branchname<<" -- outfile "<<outfile<<endl;
  TChain chain(branchname);
  chain.Add(infile);

  TFile ofile(outfile, "recreate");
  TTree* otree = chain.CopyTree(cut);

  ofile.Write();
  ofile.Close();
}

void hadd()
{
  int run=53744;

  // KshortReconstruction
  TString ksreco_infile = Form("./Reconstructed/%d/*%d*ks_ksreco.root",run,run);
  TString ksreco_outfile = Form("all%d_ksreco.root",run);
  TCut ksreco_cut = "crossing1 == crossing2 && sqrt(px1*px1+py1*py1) > 0.20 && sqrt(px2*px2+py2*py2) > 0.20 && cosThetaReco>0.5 && fabs(projected_pair_dca)<1";

  // KFparticle Kshort
  TString kfp_ks_unlikesign_infile = Form("./Reconstructed/%d/*%d*ks_kfp_unlikesign.root",run,run);
  TString kfp_ks_unlikesign_outfile = Form("all%d_kfp_ks_unlikesign.root",run);
  TString kfp_ks_likesign_infile = Form("./Reconstructed/%d/*%d*ks_kfp_likesign.root",run,run);
  TString kfp_ks_likesign_outfile = Form("all%d_kfp_ks_likesign.root",run);
  TCut kfp_ks_cut = "track_1_MVTX_nHits>0 && track_1_INTT_nHits>0 && track_2_MVTX_nHits>0 && track_2_INTT_nHits>0 && track_1_pT>0.2 && track_2_pT>0.2 && K_S0_DIRA>0.5 && fabs(track_1_track_2_DCA)<1";

  // KFparticle D0
  TString kfp_D0_unlikesign_infile = Form("./Reconstructed/%d/*%d*D0_kfp_unlikesign.root",run,run);
  TString kfp_D0_unlikesign_outfile = Form("all%d_kfp_D0_unlikesign.root",run);
  TString kfp_D0_likesign_infile = Form("./Reconstructed/%d/*%d*D0_kfp_likesign.root",run,run);
  TString kfp_D0_likesign_outfile = Form("all%d_kfp_D0_likesign.root",run);
  TCut kfp_D0_cut = "track_1_MVTX_nHits>0 && track_1_INTT_nHits>0 && track_2_MVTX_nHits>0 && track_2_INTT_nHits>0 && track_1_pT>0.2 && track_2_pT>0.2 && fabs(track_1_track_2_DCA)<1";

  combine(ksreco_infile, "ntp_reco_info", ksreco_outfile, ksreco_cut);
  combine(kfp_ks_unlikesign_infile, "DecayTree", kfp_ks_unlikesign_outfile, kfp_ks_cut);
  combine(kfp_ks_likesign_infile, "DecayTree", kfp_ks_likesign_outfile, kfp_ks_cut);
  combine(kfp_D0_unlikesign_infile, "DecayTree", kfp_D0_unlikesign_outfile, kfp_D0_cut);
  combine(kfp_D0_likesign_infile, "DecayTree", kfp_D0_likesign_outfile, kfp_D0_cut);

  std::cout << "Done!" << std::endl;
}
