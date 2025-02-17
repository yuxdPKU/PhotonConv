void read()
{
  TFile* file = new TFile("trackanalysis.root");
  TTree* tree = (TTree*) file->Get("tree");

  int _runNumber=0;
  int _eventNumber=0;
  std::vector<float> *_track_ptq=0;
  std::vector<float> *_track_px=0;
  std::vector<float> *_track_py=0;
  std::vector<float> *_track_pz=0;
  std::vector<float> *_trClus_track_id=0;
  std::vector<float> *_trClus_type=0;
  std::vector<float> *_trClus_x=0;
  std::vector<float> *_trClus_y=0;
  std::vector<float> *_trClus_z=0;

  tree->SetBranchAddress("_runNumber",&_runNumber);
  tree->SetBranchAddress("_eventNumber",&_eventNumber);
  tree->SetBranchAddress("_track_ptq",&_track_ptq);
  tree->SetBranchAddress("_track_px",&_track_px);
  tree->SetBranchAddress("_track_py",&_track_py);
  tree->SetBranchAddress("_track_pz",&_track_pz);
  tree->SetBranchAddress("_trClus_track_id",&_trClus_track_id);
  tree->SetBranchAddress("_trClus_type",&_trClus_type);
  tree->SetBranchAddress("_trClus_x",&_trClus_x);
  tree->SetBranchAddress("_trClus_y",&_trClus_y);
  tree->SetBranchAddress("_trClus_z",&_trClus_z);

  int nevent = tree->GetEntries();
  for (int i=0; i<nevent; i++)
  {
    tree->GetEntry(i);
    cout<<"run "<<_runNumber<<" , event "<<_eventNumber<<endl;
    int ntrack = _track_ptq->size();
    cout<<"Number of Tracks = "<<ntrack<<endl;
    cout<<"Number of clusters = "<<_trClus_track_id->size()<<endl;
    for (int j=0; j<ntrack; j++)
    {
      cout<<"track "<<j<<", ptq = "<<_track_ptq->at(j)<<", id = "<<endl;
    }
    cout<<endl;
  }

}
