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

#include "utilities.h"

void event_file_start(std::ofstream &jason_file_header, std::string date, int runid, int evtid);

void EvtDisplay()
{
    int runid = 46730;
    int evtid = 6747;

    //TFile * fin = new TFile("/sphenix/u/xyu3/hftg01/PhotonConv/macro/Reconstructed/46730/final_46730_ana.root");
    //fin->cd();
    //TTree * tree = (TTree *) fin->Get("tree");

    TChain* chain = new TChain("tree");
    chain->Add(Form("../%d/TrackCalo_*_ana.root",runid));
    setBranch(chain);

    chain->GetEntry(evtid);

    std::ostringstream oss;
    oss << "TrackCalo_run" << runid << "_event" << _eventNumber << ".json";

    std::string outfile = oss.str();
    std::cout << outfile << std::endl;

    std::ofstream outdata;
    outdata.open(outfile);

    std::string date = "2024-06-25";
    event_file_start(outdata, date, runid, _eventNumber);

    bool displayHCal = true;

    outdata << "     \"HITS\": {\n   \n     \"CEMC\":  [";

    int ClusSize =  _cluster_x->size();
    int EMCALSize = _emcal_tower_e->size();
    int HCALSize = _hcal_tower_e->size();

    bool firstOHCALHits = true;
    bool firstIHCALHits = true;
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

    if (!displayHCal)
    {
      outdata <<  "\"HCALIN\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n            ],\n        \"HCALOUT\": [{\"eta\": 0, \"phi\": 0, \"e\": 0}\n \n            ],\n\n" << endl;
    }
    else
    {
      outdata << "\"HCALIN\":  [";
      for (int i = 0; i < HCALSize; i++)
      {
          if (_hcal_tower_io->at(i)!=1)
          {
            continue;
          }
          if (isnan(_hcal_tower_eta->at(i))) continue;
          if (isnan(_hcal_tower_phi->at(i))) continue;
          if (isnan(_hcal_tower_e->at(i))) continue;
          std::ostringstream spts;
          if (firstIHCALHits)
        	firstIHCALHits = false;
          else
        	spts << ",";

          spts << "{ \"eta\": ";
          spts << _hcal_tower_eta->at(i);
          spts << ", \"phi\": ";
          spts << _hcal_tower_phi->at(i);
          spts << ", \"e\": ";
          spts << _hcal_tower_e->at(i);
          spts << ", \"event\": ";
          spts << 0; // adcMax?
          spts << "}";

          outdata << (boost::format("%1%") % spts.str());
          spts.clear();
          spts.str("");
      }
      outdata << "],\n\n " << endl;

      outdata << "\"HCALOUT\":  [";
      for (int i = 0; i < HCALSize; i++)
      {
          if (_hcal_tower_io->at(i)!=2)
          {
            continue;
          }
          if (isnan(_hcal_tower_eta->at(i))) continue;
          if (isnan(_hcal_tower_phi->at(i))) continue;
          if (isnan(_hcal_tower_e->at(i))) continue;
          std::ostringstream spts;
          if (firstOHCALHits)
        	firstOHCALHits = false;
          else
        	spts << ",";

          spts << "{ \"eta\": ";
          spts << _hcal_tower_eta->at(i);
          spts << ", \"phi\": ";
          spts << _hcal_tower_phi->at(i);
          spts << ", \"e\": ";
          spts << _hcal_tower_e->at(i);
          spts << ", \"event\": ";
          spts << 0; // adcMax?
          spts << "}";

          outdata << (boost::format("%1%") % spts.str());
          spts.clear();
          spts.str("");
      }
      outdata << "],\n\n " << endl;
    }

    outdata << "    \"TRACKHITS\": [\n\n ";

    for (int i = 0; i < ClusSize; i++)
    {
        if (isnan(_cluster_x->at(i))) continue;
        if (isnan(_cluster_y->at(i))) continue;
        if (isnan(_cluster_z->at(i))) continue;
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
               if (isnan(_cluster_x->at(i))) continue;
               if (isnan(_cluster_y->at(i))) continue;
               if (isnan(_cluster_z->at(i))) continue;

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
    jason_file_header << "{\n    \"EVENT\": {\n        \"runid\": " << runid << ", \n        \"evtid\": " << evtid << ", \n        \"time\": 0, \n        \"type\": \"Collision\", \n        \"s_nn\": 0, \n        \"B\": 3.0,\n        \"pv\": [0,0,0],\n        \"runstats\": [\"sPHENIX Internal\",        \n        \"200 GeV pp\",        \n        \"" << date <<", Run " << runid << ", Event #" << evtid << "\"]  \n    },\n" << endl;

    jason_file_header << "    \"META\": {\n       \"HITS\": {\n          \"INNERTRACKER\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 6.0,\n              \"color\": 16711680\n              } \n          },\n" << endl;
    jason_file_header << "          \"TRACKHITS\": {\n              \"type\": \"3D\",\n              \"options\": {\n              \"size\": 2.0,\n              \"transparent\": 0.6,\n              \"color\": 16777215\n              } \n          },\n" << endl;
    jason_file_header << "          \"CEMC\": {\n              \"type\": \"PROJECTIVE\",\n              \"options\": {\n                  \"rmin\": 90,\n                  \"rmax\": 136.1,\n                  \"deta\": 0.025,\n                  \"dphi\": 0.025,\n                  \"color\": 16766464,\n                  \"transparent\": 0.6,\n                  \"scaleminmax\": true\n              }\n          },\n" << endl;
    jason_file_header << "          \"HCALIN\": {\n              \"type\": \"PROJECTIVE\",\n              \"options\": {\n                  \"rmin\": 147.27,\n                  \"rmax\": 175.0,\n                  \"deta\": 0.025,\n                  \"dphi\": 0.025,\n                  \"color\": 4290445312,\n                  \"transparent\": 0.6,\n                  \"scaleminmax\": true\n              }\n          },\n" << endl;
    jason_file_header << "          \"HCALOUT\": {\n              \"type\": \"PROJECTIVE\",\n              \"options\": {\n                  \"rmin\": 183.3,\n                  \"rmax\": 348.634,\n                  \"deta\": 0.025,\n                  \"dphi\": 0.025,\n                  \"color\": 24773,\n                  \"transparent\": 0.6,\n                  \"scaleminmax\": true\n              }\n          },\n" << endl;
    jason_file_header << "    \"JETS\": {\n        \"type\": \"JET\",\n        \"options\": {\n            \"rmin\": 0,\n            \"rmax\": 78,\n            \"emin\": 0,\n            \"emax\": 30,\n            \"color\": 16777215,\n            \"transparent\": 0.5 \n        }\n    }\n        }\n    }\n," << endl;
   
}
