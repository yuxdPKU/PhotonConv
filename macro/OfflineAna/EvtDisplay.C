#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>

#include "TChain.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

#include <boost/format.hpp>
#include <boost/math/special_functions/sign.hpp>

void event_file_start(std::ofstream &jason_file_header, std::string date, int runid, int evtid);

void EvtDisplay()
{
    int runid = 46730;
    //int evtid = 11489;
    //int evtid = 7586;
    //int evtid = 16063;
    int evtid = 6747;

    std::ostringstream oss;
    oss << "ConversionElectron_run" << runid << "_event" << evtid << ".json";

    std::string outfile = oss.str();
    std::cout << outfile << std::endl;

    std::ofstream outdata;
    outdata.open(outfile);

    std::string date = "2024-06-25";
    event_file_start(outdata, date, runid, evtid);

    TFile * fin = new TFile("/sphenix/u/xyu3/hftg01/PhotonConv/macro/Reconstructed/46730/final_46730_ana.root");
    fin->cd();
    TTree * tree = (TTree *) fin->Get("tree");

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
    std::vector<float> *_track_phi = 0;
    std::vector<float> *_track_eta = 0;
    std::vector<float> *_track_x = 0;
    std::vector<float> *_track_y = 0;
    std::vector<float> *_track_z = 0;
    std::vector<float> *_track_ptq = 0;
    std::vector<float> *_track_px = 0;
    std::vector<float> *_track_py = 0;
    std::vector<float> *_track_pz = 0;
    std::vector<float> *_track_phi_origin = 0;
    std::vector<float> *_track_eta_origin = 0;
    std::vector<float> *_track_x_origin = 0;
    std::vector<float> *_track_y_origin = 0;
    std::vector<float> *_track_z_origin = 0;
    std::vector<float> *_track_phi_emc = 0;
    std::vector<float> *_track_eta_emc = 0;
    std::vector<float> *_track_x_emc = 0;
    std::vector<float> *_track_y_emc = 0;
    std::vector<float> *_track_z_emc = 0;
    std::vector<float> *_track_phi_ihc = 0;
    std::vector<float> *_track_eta_ihc = 0;
    std::vector<float> *_track_x_ihc = 0;
    std::vector<float> *_track_y_ihc = 0;
    std::vector<float> *_track_z_ihc = 0;
    std::vector<float> *_track_phi_ohc = 0;
    std::vector<float> *_track_eta_ohc = 0;
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
    std::vector<int> *_ntracks = 0;
    std::vector<float> *_mbd_z = 0;
    std::vector<int> *_triggers = 0;

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
    tree->SetBranchAddress("_track_x", &_track_x);
    tree->SetBranchAddress("_track_y", &_track_y);
    tree->SetBranchAddress("_track_z", &_track_z);
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
    tree->SetBranchAddress("_mbd_z", &_mbd_z);
    tree->SetBranchAddress("_triggers", &_triggers);
    tree->SetBranchAddress("_ntracks", &_ntracks);

    outdata << "     \"HITS\": {\n   \n     \"CEMC\":  [";

	tree->GetEntry(evtid);

	int ClusSize =  _cluster_x->size();
	int EMCALSize = _emcal_tower_e->size();
	
	bool firstHCALHits = true;
	bool firstEMCALHits = true;
	bool firstHits = true;

	for (int i = 0; i < EMCALSize; i++)
	{
		std::ostringstream spts;
		if (firstEMCALHits)
			firstEMCALHits = false;
		else
			spts << ",";

		spts << "{ \"eta\": ";
		spts << _emcal_tower_eta->at(i);
		spts << ", \"phi\": ";
		spts << _emcal_tower_phi->at(i);
		spts << ", \"e\": ";
		spts << _emcal_tower_e->at(i);
		spts << ", \"event\": ";
		spts << 0; // adcMax?
		spts << "}";

		outdata << (boost::format("%1%") % spts.str());
		spts.clear();
		spts.str("");
	}
	outdata << "],\n " << endl;

	outdata <<  "\"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n" << endl;

	outdata << "    \"TRACKHITS\": [\n\n ";

	for (int i = 0; i < ClusSize; i++)
	{
		std::ostringstream spts;
		if (firstHits)
			firstHits = false;
		else
			spts << ",";

		spts << "{ \"x\": ";
		spts << _cluster_x->at(i);
		spts << ", \"y\": ";
		spts << _cluster_y->at(i);
		spts << ", \"z\": ";
		spts << _cluster_z->at(i);
		spts << ", \"e\": ";
		spts << 0; // adcMax?
		spts << "}";

		outdata << (boost::format("%1%") % spts.str());
		spts.clear();
		spts.str("");
	}

    std::vector<float> vec_track_cluster_x;
    std::vector<float> vec_track_cluster_y;
    std::vector<float> vec_track_cluster_z;
    vec_track_cluster_x.clear();
    vec_track_cluster_y.clear();
    vec_track_cluster_z.clear();

    for(unsigned int itrack = 0; itrack < _track_ptq->size(); itrack++)
    {
      // cut good tpc tracks
      if(_track_nc_tpc->at(itrack) < 22) continue;
      float track_p = sqrt(pow(_track_px->at(itrack),2) + pow(_track_py->at(itrack),2) + pow(_track_pz->at(itrack),2));
      if(track_p < 1.5) continue;

      if(isnan(_track_phi_emc->at(itrack)))
      {
        continue;
      }

      // project tpc tracks to emcal and hcal
      std::pair<float, float> TrackProjsEMCal;
      TrackProjsEMCal = std::make_pair(_track_phi_emc->at(itrack), _track_z_emc->at(itrack));

      // loop all emcal clusters to match with tpc tracks
      for(unsigned int iem = 0; iem < _emcal_e->size(); iem++)
      {
        if (_emcal_e->at(iem)<1)
        {
          continue;
        }
        std::pair<float, float> EMCalPos;
        EMCalPos = std::make_pair(_emcal_phi->at(iem), _emcal_z->at(iem));

        float dphi = TrackProjsEMCal.first - EMCalPos.first;
        float dz = TrackProjsEMCal.second - EMCalPos.second;
        if(fabs(dphi)<0.1 && fabs(dz)<20)
        {
          for(unsigned int ic = 0; ic < _trClus_track_id->size(); ic++)
          {
            if(_track_id->at(itrack) != _trClus_track_id->at(ic)) continue;
            vec_track_cluster_x.push_back(_trClus_x->at(ic));
            vec_track_cluster_y.push_back(_trClus_y->at(ic));
            vec_track_cluster_z.push_back(_trClus_z->at(ic));
          }
        }
      }
    }

    int MatchedClusSize = vec_track_cluster_x.size();

    outdata << "],\n    \"JETS\": [\n         ],\n\n" << endl;

    outdata << "    \"INNERTRACKER\": [\n\n ";

	firstHits = true;
	for (int i = 0; i < MatchedClusSize; i++)
	{

		std::ostringstream spts;
		if (firstHits)
			firstHits = false;
		else
			spts << ",";

		spts << "{ \"x\": ";
		spts << vec_track_cluster_x.at(i);
		spts << ", \"y\": ";
		spts << vec_track_cluster_y.at(i);
		spts << ", \"z\": ";
		spts << vec_track_cluster_z.at(i);
		spts << ", \"e\": ";
		spts << 0; // adcMax?
		spts << "}";

		outdata << (boost::format("%1%") % spts.str());
		spts.clear();
		spts.str("");
	}

    outdata << "]" << endl;
    outdata << "}" << endl;

    outdata << "}" << endl;

    outdata.close();
}

void event_file_start(std::ofstream &jason_file_header, std::string date, int runid, int evtid)
{
    jason_file_header << "{\n    \"EVENT\": {\n        \"runid\": " << runid << ", \n        \"evtid\": " << evtid << ", \n        \"time\": 0, \n        \"type\": \"Collision\", \n        \"s_nn\": 0, \n        \"B\": 3.0,\n        \"pv\": [0,0,0],\n        \"runstats\": [\"sPHENIX Internal\",        \n\"200 GeV pp\",        \n\"2024-06-25, Run " << runid << "\",        \n\"EMCAL-TPC Matched Tracks\"]  \n    },\n" << endl;

    jason_file_header << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 6.0,\n              \"color\": 16711680\n              } \n          },\n" << endl;
    jason_file_header << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2.0,\n              \"transparent\": 0.6,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    jason_file_header << "          \"CEMC\": {\n              \"type\": \"PROJECTIVE\",\n              \"options\": {\n                  \"rmin\": 90,\n                  \"rmax\": 136.1,\n                  \"deta\": 0.025,\n                  \"dphi\": 0.025,\n                  \"color\": 16766464,\n                  \"transparent\": 0.6,\n                  \"scaleminmax\": true\n              }\n          },\n" << endl;
    jason_file_header << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << endl;
   
}
