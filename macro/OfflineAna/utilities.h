float PiRange(float deltaPhi);
bool IsClusterClose(float projPhi, float projEta, float clusterPhi, float clusterEta, float minR);
bool HasOverlap(float projPhi, float projEta, std::vector<float> towerPhi, std::vector<float> towerEta, float minDelta);
void DrawEllipseWithTGraph(double x_center, double y_center, double radius, int color);
void DrawVLineWithTGraph(double x_center, double y_min, double y_max, int color);
void DrawHLineWithTGraph(double y_center, double x_min, double x_max, int color);

int _EMCAL_NETA = 96;
int _EMCAL_NPHI = 256;
int _HCAL_NETA = 24;
int _HCAL_NPHI = 64;

float min_Track_Pt = 0.5;
float min_EMCal_E = 0.2;
float min_HCal_E = 0.2;

float emcal_radius = 100.70;//(1-(-0.077))*93.5
//float emcal_radius = 99.1;//(1-(-0.060))*93.5
float hcal_radius = 177.423;

// branch vars
  int _runNumber = 0;
  int _eventNumber = 0;
  std::vector<int> *_vertex_id = 0;
  std::vector<int> *_vertex_crossing = 0;
  std::vector<int> *_vertex_ntracks = 0;
  std::vector<float> *_vertex_x = 0;
  std::vector<float> *_vertex_y = 0;
  std::vector<float> *_vertex_z = 0;
  std::vector<float> *_cluster_x = 0;
  std::vector<float> *_cluster_y = 0;
  std::vector<float> *_cluster_z = 0;
  std::vector<int> *_track_id = 0;
  std::vector<float> *_track_quality = 0;
  std::vector<float> *_track_dcaxy = 0;
  std::vector<float> *_track_dcaz = 0;
  std::vector<int> *_track_nc_tpc = 0;
  std::vector<int> *_track_nc_mvtx = 0;
  std::vector<int> *_track_nc_intt = 0;
  std::vector<int> *_track_bc = 0;
  std::vector<int> *_track_crossing = 0;
  std::vector<float> *_track_phi = 0;
  std::vector<float> *_track_eta = 0;
  std::vector<float> *_track_pcax = 0;
  std::vector<float> *_track_pcay = 0;
  std::vector<float> *_track_pcaz = 0;
  std::vector<float> *_track_vx = 0;
  std::vector<float> *_track_vy = 0;
  std::vector<float> *_track_vz = 0;
  std::vector<float> *_track_ptq = 0;
  std::vector<float> *_track_px = 0;
  std::vector<float> *_track_py = 0;
  std::vector<float> *_track_pz = 0;
  std::vector<float> *_track_phi_origin = 0;
  std::vector<float> *_track_eta_origin = 0;
  std::vector<float> *_track_px_origin = 0;
  std::vector<float> *_track_py_origin = 0;
  std::vector<float> *_track_pz_origin = 0;
  std::vector<float> *_track_x_origin = 0;
  std::vector<float> *_track_y_origin = 0;
  std::vector<float> *_track_z_origin = 0;
  std::vector<float> *_track_phi_emc = 0;
  std::vector<float> *_track_eta_emc = 0;
  std::vector<float> *_track_px_emc = 0;
  std::vector<float> *_track_py_emc = 0;
  std::vector<float> *_track_pz_emc = 0;
  std::vector<float> *_track_x_emc = 0;
  std::vector<float> *_track_y_emc = 0;
  std::vector<float> *_track_z_emc = 0;
  std::vector<float> *_track_phi_ihc = 0;
  std::vector<float> *_track_eta_ihc = 0;
  std::vector<float> *_track_px_ihc = 0;
  std::vector<float> *_track_py_ihc = 0;
  std::vector<float> *_track_pz_ihc = 0;
  std::vector<float> *_track_x_ihc = 0;
  std::vector<float> *_track_y_ihc = 0;
  std::vector<float> *_track_z_ihc = 0;
  std::vector<float> *_track_phi_ohc = 0;
  std::vector<float> *_track_eta_ohc = 0;
  std::vector<float> *_track_px_ohc = 0;
  std::vector<float> *_track_py_ohc = 0;
  std::vector<float> *_track_pz_ohc = 0;
  std::vector<float> *_track_x_ohc = 0;
  std::vector<float> *_track_y_ohc = 0;
  std::vector<float> *_track_z_ohc = 0;
  std::vector<int> *_trClus_track_id = 0;
  std::vector<int> *_trClus_type = 0;
  std::vector<float> *_trClus_x = 0;
  std::vector<float> *_trClus_y = 0;
  std::vector<float> *_trClus_z = 0;
  std::vector<int> *_emcal_id = 0;
  std::vector<float> *_emcal_phi = 0;
  std::vector<float> *_emcal_eta = 0;
  std::vector<float> *_emcal_x = 0;
  std::vector<float> *_emcal_y = 0;
  std::vector<float> *_emcal_z = 0;
  std::vector<float> *_emcal_e = 0;
  std::vector<float> *_emcal_ecore = 0;
  std::vector<float> *_emcal_chi2 = 0;
  std::vector<float> *_emcal_prob = 0;
  std::vector<int> *_emcal_tower_cluster_id = 0;
  std::vector<float> *_emcal_tower_e = 0;
  std::vector<float> *_emcal_tower_phi = 0;
  std::vector<float> *_emcal_tower_eta = 0;
  std::vector<int> *_emcal_tower_status = 0;
  std::vector<int> *_hcal_id = 0;
  std::vector<float> *_hcal_phi = 0;
  std::vector<float> *_hcal_eta = 0;
  std::vector<float> *_hcal_x = 0;
  std::vector<float> *_hcal_y = 0;
  std::vector<float> *_hcal_z = 0;
  std::vector<float> *_hcal_e = 0;
  std::vector<int> *_hcal_tower_cluster_id = 0;
  std::vector<float> *_hcal_tower_e = 0;
  std::vector<float> *_hcal_tower_phi = 0;
  std::vector<float> *_hcal_tower_eta = 0;
  std::vector<int> *_hcal_tower_status = 0;
  std::vector<int> *_hcal_tower_io = 0;
  std::vector<int> *_ntracks = 0;
  std::vector<float> *_mbd_x = 0;
  std::vector<float> *_mbd_y = 0;
  std::vector<float> *_mbd_z = 0;
  std::vector<int> *_triggers = 0;

  int _numCan = 0;
  std::vector<float> *_gamma_mass = 0;
  std::vector<float> *_gamma_massErr = 0;
  std::vector<float> *_gamma_x = 0;
  std::vector<float> *_gamma_y = 0;
  std::vector<float> *_gamma_z = 0;
  std::vector<float> *_gamma_px = 0;
  std::vector<float> *_gamma_py = 0;
  std::vector<float> *_gamma_pz = 0;
  std::vector<float> *_gamma_pE = 0;
  std::vector<float> *_gamma_pT = 0;
  std::vector<float> *_gamma_pTErr = 0;
  std::vector<float> *_gamma_p = 0;
  std::vector<float> *_gamma_pErr = 0;
  std::vector<float> *_gamma_pseudorapidity = 0;
  std::vector<float> *_gamma_rapidity = 0;
  std::vector<float> *_gamma_theta = 0;
  std::vector<float> *_gamma_phi = 0;
  std::vector<float> *_gamma_chi2 = 0;
  std::vector<float> *_gamma_nDoF = 0;
  std::vector<float> *_gamma_vertex_volume = 0;
  std::vector<float> *_gamma_SV_chi2_per_nDoF = 0;
  std::vector<float> *_ep_mass = 0;
  std::vector<float> *_ep_x = 0;
  std::vector<float> *_ep_y = 0;
  std::vector<float> *_ep_z = 0;
  std::vector<float> *_ep_px = 0;
  std::vector<float> *_ep_py = 0;
  std::vector<float> *_ep_pz = 0;
  std::vector<float> *_ep_pE = 0;
  std::vector<float> *_ep_pE_unmoved = 0;
  std::vector<float> *_ep_pT = 0;
  std::vector<float> *_ep_pTErr = 0;
  std::vector<float> *_ep_pT_raw = 0;
  std::vector<float> *_ep_pT_unmoved = 0;
  std::vector<float> *_ep_p = 0;
  std::vector<float> *_ep_pErr = 0;
  std::vector<float> *_ep_p_raw = 0;
  std::vector<float> *_ep_p_unmoved = 0;
  std::vector<float> *_ep_pseudorapidity = 0;
  std::vector<float> *_ep_rapidity = 0;
  std::vector<float> *_ep_theta = 0;
  std::vector<float> *_ep_phi = 0;
  std::vector<float> *_ep_chi2 = 0;
  std::vector<float> *_ep_nDoF = 0;
  std::vector<int> *_ep_crossing = 0;
  std::vector<int> *_ep_clus_ican = 0;
  std::vector<float> *_ep_clus_x = 0;
  std::vector<float> *_ep_clus_y = 0;
  std::vector<float> *_ep_clus_z = 0;
  std::vector<float> *_ep_phi_emc = 0;
  std::vector<float> *_ep_eta_emc = 0;
  std::vector<float> *_ep_px_emc = 0;
  std::vector<float> *_ep_py_emc = 0;
  std::vector<float> *_ep_pz_emc = 0;
  std::vector<float> *_ep_x_emc = 0;
  std::vector<float> *_ep_y_emc = 0;
  std::vector<float> *_ep_z_emc = 0;
  std::vector<float> *_em_mass = 0;
  std::vector<float> *_em_x = 0;
  std::vector<float> *_em_y = 0;
  std::vector<float> *_em_z = 0;
  std::vector<float> *_em_px = 0;
  std::vector<float> *_em_py = 0;
  std::vector<float> *_em_pz = 0;
  std::vector<float> *_em_pE = 0;
  std::vector<float> *_em_pE_unmoved = 0;
  std::vector<float> *_em_pT = 0;
  std::vector<float> *_em_pTErr = 0;
  std::vector<float> *_em_pT_raw = 0;
  std::vector<float> *_em_pT_unmoved = 0;
  std::vector<float> *_em_p = 0;
  std::vector<float> *_em_pErr = 0;
  std::vector<float> *_em_p_raw = 0;
  std::vector<float> *_em_p_unmoved = 0;
  std::vector<float> *_em_pseudorapidity = 0;
  std::vector<float> *_em_rapidity = 0;
  std::vector<float> *_em_theta = 0;
  std::vector<float> *_em_phi = 0;
  std::vector<float> *_em_chi2 = 0;
  std::vector<float> *_em_nDoF = 0;
  std::vector<int> *_em_crossing = 0;
  std::vector<int> *_em_clus_ican = 0;
  std::vector<float> *_em_clus_x = 0;
  std::vector<float> *_em_clus_y = 0;
  std::vector<float> *_em_clus_z = 0;
  std::vector<float> *_em_phi_emc = 0;
  std::vector<float> *_em_eta_emc = 0;
  std::vector<float> *_em_px_emc = 0;
  std::vector<float> *_em_py_emc = 0;
  std::vector<float> *_em_pz_emc = 0;
  std::vector<float> *_em_x_emc = 0;
  std::vector<float> *_em_y_emc = 0;
  std::vector<float> *_em_z_emc = 0;
  std::vector<float> *_epem_DCA_2d = 0;
  std::vector<float> *_epem_DCA_3d = 0;

void setBranch(TTree* tree)
{
  tree->SetBranchAddress("_runNumber", &_runNumber);
  tree->SetBranchAddress("_eventNumber", &_eventNumber);
  tree->SetBranchAddress("_vertex_id", &_vertex_id);
  tree->SetBranchAddress("_vertex_crossing", &_vertex_crossing);
  tree->SetBranchAddress("_vertex_ntracks", &_vertex_ntracks);
  tree->SetBranchAddress("_vertex_x", &_vertex_x);
  tree->SetBranchAddress("_vertex_y", &_vertex_y);
  tree->SetBranchAddress("_vertex_z", &_vertex_z);
  tree->SetBranchAddress("_cluster_x", &_cluster_x);
  tree->SetBranchAddress("_cluster_y", &_cluster_y);
  tree->SetBranchAddress("_cluster_z", &_cluster_z);
  tree->SetBranchAddress("_track_id", &_track_id);
  tree->SetBranchAddress("_track_quality", &_track_quality);
  tree->SetBranchAddress("_track_dcaxy", &_track_dcaxy);
  tree->SetBranchAddress("_track_dcaz", &_track_dcaz);
  tree->SetBranchAddress("_track_nc_tpc", &_track_nc_tpc);
  tree->SetBranchAddress("_track_nc_mvtx", &_track_nc_mvtx);
  tree->SetBranchAddress("_track_nc_intt", &_track_nc_intt);
  tree->SetBranchAddress("_track_bc", &_track_bc);
  tree->SetBranchAddress("_track_phi", &_track_phi);
  tree->SetBranchAddress("_track_eta", &_track_eta);
  tree->SetBranchAddress("_track_pcax", &_track_pcax);
  tree->SetBranchAddress("_track_pcay", &_track_pcay);
  tree->SetBranchAddress("_track_pcaz", &_track_pcaz);
  tree->SetBranchAddress("_track_crossing", &_track_crossing);
  tree->SetBranchAddress("_track_vx", &_track_vx);
  tree->SetBranchAddress("_track_vy", &_track_vy);
  tree->SetBranchAddress("_track_vz", &_track_vz);
  tree->SetBranchAddress("_track_ptq", &_track_ptq);
  tree->SetBranchAddress("_track_px", &_track_px);
  tree->SetBranchAddress("_track_py", &_track_py);
  tree->SetBranchAddress("_track_pz", &_track_pz);
  tree->SetBranchAddress("_track_phi_origin", &_track_phi_origin);
  tree->SetBranchAddress("_track_eta_origin", &_track_eta_origin);
  tree->SetBranchAddress("_track_px_origin", &_track_px_origin);
  tree->SetBranchAddress("_track_py_origin", &_track_py_origin);
  tree->SetBranchAddress("_track_pz_origin", &_track_pz_origin);
  tree->SetBranchAddress("_track_x_origin", &_track_x_origin);
  tree->SetBranchAddress("_track_y_origin", &_track_y_origin);
  tree->SetBranchAddress("_track_z_origin", &_track_z_origin);
  tree->SetBranchAddress("_track_phi_emc", &_track_phi_emc);
  tree->SetBranchAddress("_track_eta_emc", &_track_eta_emc);
  tree->SetBranchAddress("_track_px_emc", &_track_px_emc);
  tree->SetBranchAddress("_track_py_emc", &_track_py_emc);
  tree->SetBranchAddress("_track_pz_emc", &_track_pz_emc);
  tree->SetBranchAddress("_track_x_emc", &_track_x_emc);
  tree->SetBranchAddress("_track_y_emc", &_track_y_emc);
  tree->SetBranchAddress("_track_z_emc", &_track_z_emc);
  tree->SetBranchAddress("_track_phi_ihc", &_track_phi_ihc);
  tree->SetBranchAddress("_track_eta_ihc", &_track_eta_ihc);
  tree->SetBranchAddress("_track_px_ihc", &_track_px_ihc);
  tree->SetBranchAddress("_track_py_ihc", &_track_py_ihc);
  tree->SetBranchAddress("_track_pz_ihc", &_track_pz_ihc);
  tree->SetBranchAddress("_track_x_ihc", &_track_x_ihc);
  tree->SetBranchAddress("_track_y_ihc", &_track_y_ihc);
  tree->SetBranchAddress("_track_z_ihc", &_track_z_ihc);
  tree->SetBranchAddress("_track_phi_ohc", &_track_phi_ohc);
  tree->SetBranchAddress("_track_eta_ohc", &_track_eta_ohc);
  tree->SetBranchAddress("_track_px_ohc", &_track_px_ohc);
  tree->SetBranchAddress("_track_py_ohc", &_track_py_ohc);
  tree->SetBranchAddress("_track_pz_ohc", &_track_pz_ohc);
  tree->SetBranchAddress("_track_x_ohc", &_track_x_ohc);
  tree->SetBranchAddress("_track_y_ohc", &_track_y_ohc);
  tree->SetBranchAddress("_track_z_ohc", &_track_z_ohc);
  tree->SetBranchAddress("_trClus_track_id", &_trClus_track_id);
  tree->SetBranchAddress("_trClus_type", &_trClus_type);
  tree->SetBranchAddress("_trClus_x", &_trClus_x);
  tree->SetBranchAddress("_trClus_y", &_trClus_y);
  tree->SetBranchAddress("_trClus_z", &_trClus_z);
  tree->SetBranchAddress("_emcal_id", &_emcal_id);
  tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
  tree->SetBranchAddress("_emcal_eta", &_emcal_eta);
  tree->SetBranchAddress("_emcal_x", &_emcal_x);
  tree->SetBranchAddress("_emcal_y", &_emcal_y);
  tree->SetBranchAddress("_emcal_z", &_emcal_z);
  tree->SetBranchAddress("_emcal_e", &_emcal_e);
  tree->SetBranchAddress("_emcal_ecore", &_emcal_ecore);
  tree->SetBranchAddress("_emcal_chi2", &_emcal_chi2);
  tree->SetBranchAddress("_emcal_prob", &_emcal_prob);
  tree->SetBranchAddress("_emcal_tower_cluster_id", &_emcal_tower_cluster_id);
  tree->SetBranchAddress("_emcal_tower_e", &_emcal_tower_e);
  tree->SetBranchAddress("_emcal_tower_phi", &_emcal_tower_phi);
  tree->SetBranchAddress("_emcal_tower_eta", &_emcal_tower_eta);
  tree->SetBranchAddress("_emcal_tower_status", &_emcal_tower_status);
  tree->SetBranchAddress("_hcal_id", &_hcal_id);
  tree->SetBranchAddress("_hcal_phi", &_hcal_phi);
  tree->SetBranchAddress("_hcal_eta", &_hcal_eta);
  tree->SetBranchAddress("_hcal_x", &_hcal_x);
  tree->SetBranchAddress("_hcal_y", &_hcal_y);
  tree->SetBranchAddress("_hcal_z", &_hcal_z);
  tree->SetBranchAddress("_hcal_e", &_hcal_e);
  tree->SetBranchAddress("_hcal_tower_cluster_id", &_hcal_tower_cluster_id);
  tree->SetBranchAddress("_hcal_tower_e", &_hcal_tower_e);
  tree->SetBranchAddress("_hcal_tower_phi", &_hcal_tower_phi);
  tree->SetBranchAddress("_hcal_tower_eta", &_hcal_tower_eta);
  tree->SetBranchAddress("_hcal_tower_status", &_hcal_tower_status);
  tree->SetBranchAddress("_hcal_tower_io", &_hcal_tower_io);
  tree->SetBranchAddress("_mbd_x", &_mbd_x);
  tree->SetBranchAddress("_mbd_y", &_mbd_y);
  tree->SetBranchAddress("_mbd_z", &_mbd_z);
  tree->SetBranchAddress("_triggers", &_triggers);
  tree->SetBranchAddress("_ntracks", &_ntracks);
}

void setBranch(TChain* tree)
{
  tree->SetBranchAddress("_runNumber", &_runNumber);
  tree->SetBranchAddress("_eventNumber", &_eventNumber);
  tree->SetBranchAddress("_vertex_id", &_vertex_id);
  tree->SetBranchAddress("_vertex_crossing", &_vertex_crossing);
  tree->SetBranchAddress("_vertex_ntracks", &_vertex_ntracks);
  tree->SetBranchAddress("_vertex_x", &_vertex_x);
  tree->SetBranchAddress("_vertex_y", &_vertex_y);
  tree->SetBranchAddress("_vertex_z", &_vertex_z);
  tree->SetBranchAddress("_cluster_x", &_cluster_x);
  tree->SetBranchAddress("_cluster_y", &_cluster_y);
  tree->SetBranchAddress("_cluster_z", &_cluster_z);
  tree->SetBranchAddress("_track_id", &_track_id);
  tree->SetBranchAddress("_track_quality", &_track_quality);
  tree->SetBranchAddress("_track_dcaxy", &_track_dcaxy);
  tree->SetBranchAddress("_track_dcaz", &_track_dcaz);
  tree->SetBranchAddress("_track_nc_tpc", &_track_nc_tpc);
  tree->SetBranchAddress("_track_nc_mvtx", &_track_nc_mvtx);
  tree->SetBranchAddress("_track_nc_intt", &_track_nc_intt);
  tree->SetBranchAddress("_track_bc", &_track_bc);
  tree->SetBranchAddress("_track_phi", &_track_phi);
  tree->SetBranchAddress("_track_eta", &_track_eta);
  tree->SetBranchAddress("_track_pcax", &_track_pcax);
  tree->SetBranchAddress("_track_pcay", &_track_pcay);
  tree->SetBranchAddress("_track_pcaz", &_track_pcaz);
  tree->SetBranchAddress("_track_vx", &_track_vx);
  tree->SetBranchAddress("_track_vy", &_track_vy);
  tree->SetBranchAddress("_track_vz", &_track_vz);
  tree->SetBranchAddress("_track_ptq", &_track_ptq);
  tree->SetBranchAddress("_track_px", &_track_px);
  tree->SetBranchAddress("_track_py", &_track_py);
  tree->SetBranchAddress("_track_pz", &_track_pz);
  tree->SetBranchAddress("_track_phi_origin", &_track_phi_origin);
  tree->SetBranchAddress("_track_eta_origin", &_track_eta_origin);
  tree->SetBranchAddress("_track_x_origin", &_track_x_origin);
  tree->SetBranchAddress("_track_y_origin", &_track_y_origin);
  tree->SetBranchAddress("_track_z_origin", &_track_z_origin);
  tree->SetBranchAddress("_track_phi_emc", &_track_phi_emc);
  tree->SetBranchAddress("_track_eta_emc", &_track_eta_emc);
  tree->SetBranchAddress("_track_x_emc", &_track_x_emc);
  tree->SetBranchAddress("_track_y_emc", &_track_y_emc);
  tree->SetBranchAddress("_track_z_emc", &_track_z_emc);
  tree->SetBranchAddress("_track_phi_ihc", &_track_phi_ihc);
  tree->SetBranchAddress("_track_eta_ihc", &_track_eta_ihc);
  tree->SetBranchAddress("_track_x_ihc", &_track_x_ihc);
  tree->SetBranchAddress("_track_y_ihc", &_track_y_ihc);
  tree->SetBranchAddress("_track_z_ihc", &_track_z_ihc);
  tree->SetBranchAddress("_track_phi_ohc", &_track_phi_ohc);
  tree->SetBranchAddress("_track_eta_ohc", &_track_eta_ohc);
  tree->SetBranchAddress("_track_x_ohc", &_track_x_ohc);
  tree->SetBranchAddress("_track_y_ohc", &_track_y_ohc);
  tree->SetBranchAddress("_track_z_ohc", &_track_z_ohc);
  tree->SetBranchAddress("_trClus_track_id", &_trClus_track_id);
  tree->SetBranchAddress("_trClus_type", &_trClus_type);
  tree->SetBranchAddress("_trClus_x", &_trClus_x);
  tree->SetBranchAddress("_trClus_y", &_trClus_y);
  tree->SetBranchAddress("_trClus_z", &_trClus_z);
  tree->SetBranchAddress("_emcal_id", &_emcal_id);
  tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
  tree->SetBranchAddress("_emcal_eta", &_emcal_eta);
  tree->SetBranchAddress("_emcal_x", &_emcal_x);
  tree->SetBranchAddress("_emcal_y", &_emcal_y);
  tree->SetBranchAddress("_emcal_z", &_emcal_z);
  tree->SetBranchAddress("_emcal_e", &_emcal_e);
  tree->SetBranchAddress("_emcal_ecore", &_emcal_ecore);
  tree->SetBranchAddress("_emcal_chi2", &_emcal_chi2);
  tree->SetBranchAddress("_emcal_prob", &_emcal_prob);
  tree->SetBranchAddress("_emcal_tower_cluster_id", &_emcal_tower_cluster_id);
  tree->SetBranchAddress("_emcal_tower_e", &_emcal_tower_e);
  tree->SetBranchAddress("_emcal_tower_phi", &_emcal_tower_phi);
  tree->SetBranchAddress("_emcal_tower_eta", &_emcal_tower_eta);
  tree->SetBranchAddress("_emcal_tower_status", &_emcal_tower_status);
  tree->SetBranchAddress("_hcal_id", &_hcal_id);
  tree->SetBranchAddress("_hcal_phi", &_hcal_phi);
  tree->SetBranchAddress("_hcal_eta", &_hcal_eta);
  tree->SetBranchAddress("_hcal_x", &_hcal_x);
  tree->SetBranchAddress("_hcal_y", &_hcal_y);
  tree->SetBranchAddress("_hcal_z", &_hcal_z);
  tree->SetBranchAddress("_hcal_e", &_hcal_e);
  tree->SetBranchAddress("_hcal_tower_cluster_id", &_hcal_tower_cluster_id);
  tree->SetBranchAddress("_hcal_tower_e", &_hcal_tower_e);
  tree->SetBranchAddress("_hcal_tower_phi", &_hcal_tower_phi);
  tree->SetBranchAddress("_hcal_tower_eta", &_hcal_tower_eta);
  tree->SetBranchAddress("_hcal_tower_status", &_hcal_tower_status);
  tree->SetBranchAddress("_hcal_tower_io", &_hcal_tower_io);
  tree->SetBranchAddress("_mbd_z", &_mbd_z);
  tree->SetBranchAddress("_triggers", &_triggers);
  tree->SetBranchAddress("_ntracks", &_ntracks);
}

void setBranch_kfp(TChain* tree)
{
  tree->SetBranchAddress("_runNumber", &_runNumber);
  tree->SetBranchAddress("_eventNumber", &_eventNumber);
  tree->SetBranchAddress("_numCan", &_numCan);
  tree->SetBranchAddress("_gamma_mass", &_gamma_mass);
  tree->SetBranchAddress("_gamma_massErr", &_gamma_massErr);
  tree->SetBranchAddress("_gamma_x", &_gamma_x);
  tree->SetBranchAddress("_gamma_y", &_gamma_y);
  tree->SetBranchAddress("_gamma_z", &_gamma_z);
  tree->SetBranchAddress("_gamma_px", &_gamma_px);
  tree->SetBranchAddress("_gamma_py", &_gamma_py);
  tree->SetBranchAddress("_gamma_pz", &_gamma_pz);
  tree->SetBranchAddress("_gamma_pE", &_gamma_pE);
  tree->SetBranchAddress("_gamma_pT", &_gamma_pT);
  tree->SetBranchAddress("_gamma_pTErr", &_gamma_pTErr);
  tree->SetBranchAddress("_gamma_p", &_gamma_p);
  tree->SetBranchAddress("_gamma_pErr", &_gamma_pErr);
  tree->SetBranchAddress("_gamma_pseudorapidity", &_gamma_pseudorapidity);
  tree->SetBranchAddress("_gamma_rapidity", &_gamma_rapidity);
  tree->SetBranchAddress("_gamma_theta", &_gamma_theta);
  tree->SetBranchAddress("_gamma_phi", &_gamma_phi);
  tree->SetBranchAddress("_gamma_chi2", &_gamma_chi2);
  tree->SetBranchAddress("_gamma_nDoF", &_gamma_nDoF);
  tree->SetBranchAddress("_gamma_vertex_volume", &_gamma_vertex_volume);
  tree->SetBranchAddress("_gamma_SV_chi2_per_nDoF", &_gamma_SV_chi2_per_nDoF);
  tree->SetBranchAddress("_ep_mass", &_ep_mass);
  tree->SetBranchAddress("_ep_x", &_ep_x);
  tree->SetBranchAddress("_ep_y", &_ep_y);
  tree->SetBranchAddress("_ep_z", &_ep_z);
  tree->SetBranchAddress("_ep_px", &_ep_px);
  tree->SetBranchAddress("_ep_py", &_ep_py);
  tree->SetBranchAddress("_ep_pz", &_ep_pz);
  tree->SetBranchAddress("_ep_pE", &_ep_pE);
  tree->SetBranchAddress("_ep_pE_unmoved", &_ep_pE_unmoved);
  tree->SetBranchAddress("_ep_pT", &_ep_pT);
  tree->SetBranchAddress("_ep_pTErr", &_ep_pTErr);
  tree->SetBranchAddress("_ep_pT_raw", &_ep_pT_raw);
  tree->SetBranchAddress("_ep_pT_unmoved", &_ep_pT_unmoved);
  tree->SetBranchAddress("_ep_p", &_ep_p);
  tree->SetBranchAddress("_ep_pErr", &_ep_pErr);
  tree->SetBranchAddress("_ep_p_raw", &_ep_p_raw);
  tree->SetBranchAddress("_ep_p_unmoved", &_ep_p_unmoved);
  tree->SetBranchAddress("_ep_pseudorapidity", &_ep_pseudorapidity);
  tree->SetBranchAddress("_ep_rapidity", &_ep_rapidity);
  tree->SetBranchAddress("_ep_theta", &_ep_theta);
  tree->SetBranchAddress("_ep_phi", &_ep_phi);
  tree->SetBranchAddress("_ep_chi2", &_ep_chi2);
  tree->SetBranchAddress("_ep_nDoF", &_ep_nDoF);
  tree->SetBranchAddress("_ep_crossing", &_ep_crossing);
  tree->SetBranchAddress("_ep_clus_ican", &_ep_clus_ican);
  tree->SetBranchAddress("_ep_clus_x", &_ep_clus_x);
  tree->SetBranchAddress("_ep_clus_y", &_ep_clus_y);
  tree->SetBranchAddress("_ep_clus_z", &_ep_clus_z);
  tree->SetBranchAddress("_ep_phi_emc", &_ep_phi_emc);
  tree->SetBranchAddress("_ep_eta_emc", &_ep_eta_emc);
  tree->SetBranchAddress("_ep_px_emc", &_ep_px_emc);
  tree->SetBranchAddress("_ep_py_emc", &_ep_py_emc);
  tree->SetBranchAddress("_ep_pz_emc", &_ep_pz_emc);
  tree->SetBranchAddress("_ep_x_emc", &_ep_x_emc);
  tree->SetBranchAddress("_ep_y_emc", &_ep_y_emc);
  tree->SetBranchAddress("_ep_z_emc", &_ep_z_emc);
  tree->SetBranchAddress("_em_mass", &_em_mass);
  tree->SetBranchAddress("_em_x", &_em_x);
  tree->SetBranchAddress("_em_y", &_em_y);
  tree->SetBranchAddress("_em_z", &_em_z);
  tree->SetBranchAddress("_em_px", &_em_px);
  tree->SetBranchAddress("_em_py", &_em_py);
  tree->SetBranchAddress("_em_pz", &_em_pz);
  tree->SetBranchAddress("_em_pE", &_em_pE);
  tree->SetBranchAddress("_em_pE_unmoved", &_em_pE_unmoved);
  tree->SetBranchAddress("_em_pT", &_em_pT);
  tree->SetBranchAddress("_em_pTErr", &_em_pTErr);
  tree->SetBranchAddress("_em_pT_raw", &_em_pT_raw);
  tree->SetBranchAddress("_em_pT_unmoved", &_em_pT_unmoved);
  tree->SetBranchAddress("_em_p", &_em_p);
  tree->SetBranchAddress("_em_pErr", &_em_pErr);
  tree->SetBranchAddress("_em_p_raw", &_em_p_raw);
  tree->SetBranchAddress("_em_p_unmoved", &_em_p_unmoved);
  tree->SetBranchAddress("_em_pseudorapidity", &_em_pseudorapidity);
  tree->SetBranchAddress("_em_rapidity", &_em_rapidity);
  tree->SetBranchAddress("_em_theta", &_em_theta);
  tree->SetBranchAddress("_em_phi", &_em_phi);
  tree->SetBranchAddress("_em_chi2", &_em_chi2);
  tree->SetBranchAddress("_em_nDoF", &_em_nDoF);
  tree->SetBranchAddress("_em_crossing", &_em_crossing);
  tree->SetBranchAddress("_em_clus_ican", &_em_clus_ican);
  tree->SetBranchAddress("_em_clus_x", &_em_clus_x);
  tree->SetBranchAddress("_em_clus_y", &_em_clus_y);
  tree->SetBranchAddress("_em_clus_z", &_em_clus_z);
  tree->SetBranchAddress("_em_phi_emc", &_em_phi_emc);
  tree->SetBranchAddress("_em_eta_emc", &_em_eta_emc);
  tree->SetBranchAddress("_em_px_emc", &_em_px_emc);
  tree->SetBranchAddress("_em_py_emc", &_em_py_emc);
  tree->SetBranchAddress("_em_pz_emc", &_em_pz_emc);
  tree->SetBranchAddress("_em_x_emc", &_em_x_emc);
  tree->SetBranchAddress("_em_y_emc", &_em_y_emc);
  tree->SetBranchAddress("_em_z_emc", &_em_z_emc);
  tree->SetBranchAddress("_emcal_phi", &_emcal_phi);
  tree->SetBranchAddress("_emcal_eta", &_emcal_eta);
  tree->SetBranchAddress("_emcal_x", &_emcal_x);
  tree->SetBranchAddress("_emcal_y", &_emcal_y);
  tree->SetBranchAddress("_emcal_z", &_emcal_z);
  tree->SetBranchAddress("_emcal_e", &_emcal_e);
  tree->SetBranchAddress("_epem_DCA_2d", &_epem_DCA_2d);
  tree->SetBranchAddress("_epem_DCA_3d", &_epem_DCA_3d);
}

float PiRange(float deltaPhi)
{
  if(deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  if(deltaPhi < -M_PI) deltaPhi += 2*M_PI;

  return deltaPhi;
}
bool IsClusterClose(float projPhi, float projEta, float clusterPhi, float clusterEta, float minR)
{
  if(isnan(projPhi) || isnan(projEta))
  {
    return false;
  }
  float deltaEta = projEta - clusterEta;
  float deltaPhi = projPhi - clusterPhi;
  deltaPhi = PiRange(deltaPhi);
  float deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

  if(deltaR < minR)
  {
    //std::cout << "Cluster deltaPhi: " << deltaPhi << " deltaEta: " << deltaEta << " deltaR: " << deltaR << std::endl;
    return true;
  }

  return false;
}
bool HasOverlap(float projPhi, float projEta, std::vector<float> towerPhi, std::vector<float> towerEta, float minDelta)
{
  for(unsigned int i = 0; i < towerPhi.size(); i++)
  {
    float deltaEta = std::fabs(projEta - towerEta.at(i));
    float deltaPhi = projPhi - towerPhi.at(i);
    deltaPhi = std::fabs(PiRange(deltaPhi));
    if((deltaPhi < minDelta) && (deltaEta < minDelta))
    {
      return true;
    }
  }

  return false;
}

void DrawEllipseWithTGraph(double x_center, double y_center, double radius, int color) {
    const int n_points = 1000; // Number of points to approximate the ellipse
    TGraph *gr_ellipse = new TGraph(n_points);
    for (int i = 0; i < n_points; ++i) {
        double theta = 2 * TMath::Pi() * i / n_points;
        double x = x_center + radius * TMath::Cos(theta);
        double y = y_center + radius * TMath::Sin(theta);
        gr_ellipse->SetPoint(i, x, y);
    }
    gr_ellipse->SetLineColor(color);
    gr_ellipse->SetLineWidth(2);
    gr_ellipse->Draw("L same");
}

void DrawVLineWithTGraph(double x_center, double y_min, double y_max, int color) {
    const int n_points = 1000; // Number of points to approximate the line
    TGraph *gr_line = new TGraph(n_points);
    for (int i = 0; i < n_points; ++i) {
        double x = x_center;
        double y = y_min + i * (y_max - y_min) / (double) n_points;
        gr_line->SetPoint(i, x, y);
    }
    gr_line->SetLineColor(color);
    gr_line->SetLineWidth(2);
    gr_line->Draw("L same");
}
void DrawHLineWithTGraph(double y_center, double x_min, double x_max, int color) {
    const int n_points = 1000; // Number of points to approximate the line
    TGraph *gr_line = new TGraph(n_points);
    for (int i = 0; i < n_points; ++i) {
        double x = x_min + i * (x_max - x_min) / (double) n_points;
        double y = y_center;
        gr_line->SetPoint(i, x, y);
    }
    gr_line->SetLineColor(color);
    gr_line->SetLineWidth(2);
    gr_line->Draw("L same");
}

float cal_phi(float x, float y) {
  float phi = atan2(y,x);
  return phi;
}

float cal_eta(float x, float y, float z) {
  float eta = asinh( z / sqrt( x*x + y*y ) );
  return eta;
}

double customsqrt(double x) {
    if (x < 0) {
        return -std::sqrt(-x);
    } else {
        return std::sqrt(x);
    }
}
