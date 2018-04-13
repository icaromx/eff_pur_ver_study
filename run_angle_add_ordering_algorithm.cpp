#include "TROOT.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TError.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TVector3.h"
#include "TCanvas.h"
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#define PI 3.14159265
//#include "/grid/fermiapp/products/larsoft/eigen/v3_3_3/include/eigen3/Eigen/Dense"

#include "/usr/local/Cellar/eigen/3.3.4/include/eigen3/Eigen/Dense" //Needed on MACOS
using namespace std;

struct Point {
	float x;
	float y;
	float z;
float q;
};
struct PCAResults {
  	TVector3 centroid;
	pair<TVector3,TVector3> endPoints;
	float length;
	TVector3 eVals;
	vector<TVector3> eVecs;
};
struct TrkPoint{
    double c;
    double x;
    double y;
    double z;
    double q;
};
struct by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x > b.x;
        else return a.y > b.y;
    }
};
struct reverse_by_y { 
    bool operator()(TrkPoint const &a, TrkPoint const &b) { 
        if(a.y == b.y) return a.x < b.x;
        else return a.y < b.y;
    }
};
typedef vector<TrkPoint> track_def;
typedef vector<Point> PointCloud;
void LoadPointCloud(PointCloud &points, const track_def &ord_trk);
PCAResults DoPCA(const PointCloud &points);
double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2);
vector<double> Unit_Vec(double x1,double y1,double z1);

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////MAIN PROGRAM STARTS////////////////////////////////////////

int main(int argc, char **argv){
	TFile *f_output;
	
	int track_num = atoi(argv[2]);
	//string track_id = argv[2];
	//int int_track_id = atoi(argv[2]);
	///////////////////
	//Define Parameters
	double alpha = atof(argv[3]);
	double ang_gap = atof(argv[4]);
	double ang_deg = acos(ang_gap) * 180. / PI;
	unsigned ang_points = 5;
	int min_points_trk = ang_points*2;
	
  	///////////////////
  	//Define root output file
  	//f_output = TFile::Open("TEST.root","RECREATE");
  	f_output = TFile::Open(Form("results_alpha_%d_ang_%d.root",(int)alpha, (int)ang_deg),"RECREATE");
  	TNtuple *nt_study = new TNtuple("nt_study","nt_study","run_num:ev_num:cluster_id:alpha:dist:angle:ford_x:ford_y:ford_z:lalpha_x:lalpha_y:lalpha_z");
  	TNtuple *nt_vectors = new TNtuple("nt_vectors","nt_vectors", "run_num:ev_num:cluster_id:angle:pca_x_s:pca_y_s:pca_z_s:pca_x_e:pca_y_e:pca_z_e:lalpha_x_e:lalpha_y_e:lalpha_z_e");
  	TNtuple *nt_pca_points = new TNtuple("nt_pca_points","nt_pca_points","run_num:ev_num:cluster_id:pca_x:pca_y:pca_z");
  	///////////////////
  	//OUTPUT COMMENTS
  	cout << "Using only topological cuts with this selection." << endl;

  	///////////////////
  	//READ IN MICHEL LIST
  	ifstream csv_infile("Michel_candidates_vertex_v3.csv");
  	vector<string> TrackData;
  	std::vector<std::vector<double> > Michel_candidates;
  	std::string mline;
  	while (getline(csv_infile, mline,'\n')){
   		TrackData.push_back(mline); //Get each line of the file as a string
  	}
  	int s = TrackData.size();
  	for (unsigned int i=1; i<s; ++i){
   		std::vector<double> v_michel;
    	std::size_t first_comma = TrackData[i].find(",");      // position of the end of the name of each one in the respective string
    	std::size_t second_comma = TrackData[i].find(",", first_comma + 1);
    	std::size_t third_comma = TrackData[i].find(",", second_comma + 1);
    	std::size_t fourth_comma = TrackData[i].find(",", third_comma + 1);
    	std::size_t fifth_comma = TrackData[i].find(",", fourth_comma + 1);
    	double mrun = std::stod(TrackData[i].substr(0,TrackData[i].size()));
    	double meve = std::stod(TrackData[i].substr(first_comma+1,TrackData[i].size()));
    	double mtrk = std::stod(TrackData[i].substr(second_comma+1,TrackData[i].size()));
    	double mvrX = std::stod(TrackData[i].substr(third_comma+1,TrackData[i].size()));
	    double mvrY = std::stod(TrackData[i].substr(fourth_comma+1,TrackData[i].size()));
	    double mvrZ = std::stod(TrackData[i].substr(fifth_comma+1,TrackData[i].size()));
    	//if(mvrX == -1.0 && mvrY == -1.0 && mvrZ == -1.0){
    	//	continue;
	    //}
    	v_michel.push_back(mrun);
    	v_michel.push_back(meve);
    	v_michel.push_back(mtrk);
    	v_michel.push_back(mvrX);
    	v_michel.push_back(mvrY);
    	v_michel.push_back(mvrZ);
    	Michel_candidates.push_back(v_michel);
   	}
   	int mcand_size = Michel_candidates.size();
  	///////////////////////////////////////////////////////////////////////////////////////////////////
  	std::vector< std::vector<int> > kept_tracks;
  	int counter = 0;
  	int counter_x = 0;
  	//int count = ang_count*evr_count;
	double purity, efficiency, ver_rms = 0., ver_mean = 0.;
	int total_num_tracks = 0;
	int tracks_survived_ord_alg = 0;
	int michels_survived_ord_alg = 0;
	int michel_count = 0;
	int track_selected_as_michel = 0;
	int eigenval_cut = 0;
	int ord_alg_cutouts = 0;
	int ord_alg_cutouts_michels = 0;
	std::string line;
	std::ifstream ifs(argv[1]);
		//  cout << angles[ang]  << ", " << EVRs[evr] << endl;
	while(std::getline(ifs, line)){
	  	//cout << line << endl;
	    gROOT->Reset();
	    gErrorIgnoreLevel = kError;
	    TString filename;
	    filename.Form("%s",line.c_str());    
	    TFile *infile = new TFile(filename);
	    //Extract Event Metadata
	    TTree *Trun = (TTree*)infile->Get("Trun");
	    Int_t run_num;
	    Int_t ev_num;
	    Trun->SetBranchAddress("runNo",&run_num);
	    Trun->SetBranchAddress("eventNo",&ev_num);
	    Trun->GetEntry(0);
	    //cout << "Looking at run " << run_num << " from event " << ev_num << endl;
	    //Extract Coordinate information
	    TTree *T_charge_cluster = (TTree*)infile->Get("T_charge_cluster_nfc"); 
	    Double_t cluster_id;
	    Double_t qx;
	    Double_t qy;
	    Double_t qz;
	    Double_t qc;
	    T_charge_cluster->SetBranchAddress("qx",&qx);
	    T_charge_cluster->SetBranchStatus("qx", kTRUE);
	    T_charge_cluster->SetBranchAddress("qy",&qy);
	    T_charge_cluster->SetBranchStatus("qy", kTRUE);
	    T_charge_cluster->SetBranchAddress("qz",&qz);
	    T_charge_cluster->SetBranchStatus("qz", kTRUE);
	    T_charge_cluster->SetBranchAddress("qc",&qc);
	    T_charge_cluster->SetBranchStatus("qc", kTRUE);
	    T_charge_cluster->SetBranchAddress("cluster_id", &cluster_id);
	    T_charge_cluster->SetBranchStatus("cluster_id", kTRUE);
	    int size = T_charge_cluster->GetEntries();
	    /////////////////////////////////////////////////////////////
	    //Extract Clusters///////////////////////////////////////////
	    std::vector<Int_t> clusters;
	    Int_t prev_cval;
    	for (int i = 0; i < size; ++i){
      		T_charge_cluster -> GetEntry(i);
    		if (i == 0){
        		clusters.push_back(cluster_id);
        		prev_cval = cluster_id;
      		}else if(prev_cval != cluster_id){
        		prev_cval = cluster_id;
        		clusters.push_back(cluster_id);
      		}
      	}
    	//Looking at tracks individually
    	int num_clusters = clusters.size();
   		//cout << "HERE " << num_clusters << endl;
    	//Loop through individual clusters
    	std::vector<int> event_kept_trks;
    	for (int c = 0; c < num_clusters; ++c){
    		total_num_tracks += 1;
      		int cluster = clusters[c];
      		string clus_id = to_string(ev_num) + to_string(cluster);
      		int int_clus_id  = stoi(clus_id);
      		//cout << clus_id << endl;
      		
      		track_def trk;
      		//Load every point (cluster id, x, y, z, charge) of a track into the trk object.
      		for (int i = 0; i < size; ++i){
	        	T_charge_cluster -> GetEntry(i);
	        	//Will only store information for current cluster
	        	if(cluster_id != cluster) continue;

	        	if(track_num != -1){
	        		if(int_clus_id != track_num) continue;
	        	}
		        TrkPoint tempPoint;
		        tempPoint.c = cluster_id;
		        tempPoint.x = qx;
		        tempPoint.y = qy;
		        tempPoint.z = qz;
		        tempPoint.q = qc;
		        trk.push_back(tempPoint);
	      	}

	      	//track size has to be larger than the moving window size
	      	if(trk.size() < min_points_trk + 1) continue; // #CUT
	      	//////////////////////////////////////////////////
	      	////////ORDERING ALGORITHM BEGINS/////////////////
	      	//Sort track in descending y value
	      	std::sort(trk.begin(), trk.end(), by_y());
	      	int trk_size = trk.size();
	      	track_def points_left;
	      	track_def points_gd;
	      	//Store points being tested by ordering algorithm
	      	for (int i = 1; i < trk_size; ++i){
	        	TrkPoint tempPoint;
	        	tempPoint.c = trk[i].c;
	        	tempPoint.x = trk[i].x;
	        	tempPoint.y = trk[i].y;
	        	tempPoint.z = trk[i].z;
	        	tempPoint.q = trk[i].q;
	        	points_left.push_back(tempPoint);
	      	}
	      	int pl_size = points_left.size();

	      	track_def ord_trk;
	      	//Highest y value point is the first point in the oredered track
	      	ord_trk.push_back(trk[0]);
		    double old_dist = 10000000.;
		    int low_dist_at = -1;
		    double dist;

		    int m = 0;
		    int i = 0;

		    double low_ord_y = 2000.;
		    double ang_eigV_x=0., ang_eigV_y=0., ang_eigV_z=0.;
		    double ob_temp_x, ob_temp_y, ob_temp_z;
		    double ftest_point_x, ftest_point_y, ftest_point_z; 
		    double pca_eVec_x, pca_eVec_y, pca_eVec_z; 
          	double ford_point_x, ford_point_y, ford_point_z, last_dist;
          	double low_old_dist = 1000.;
		    double ang_dotP, fangle;
		    std::vector<double> vec_pca_ang;
		    std::vector<double> vec_pca;
		    std::vector<double> vec_ob_temp;
		    std::vector<double> vec_temp;
		    std::vector<double> vec_last_ord;
		    track_def pca_points;
		    double ang_decide;

	    	while(pl_size != 0){
	    		track_def ang_chunk;
	    		vec_pca_ang.clear();
	    		if(ord_trk.size() > ang_points){
	    			for (unsigned p = ord_trk.size() - ang_points - 1; p < ord_trk.size() - 1; ++p){
		     			TrkPoint ang_tempPoint;
				        ang_tempPoint.c = ord_trk[p].c;
				        ang_tempPoint.x = ord_trk[p].x;
				        ang_tempPoint.y = ord_trk[p].y;
				        ang_tempPoint.z = ord_trk[p].z;
				        ang_tempPoint.q = ord_trk[p].q;
				        ang_chunk.push_back(ang_tempPoint);	
	    			}
		    		PointCloud ang_pointcloud;
		    		PCAResults ang_results;
		    		LoadPointCloud(ang_pointcloud, ang_chunk);
				    ang_results = DoPCA(ang_pointcloud);
					
				    //ang_eigV_x = ang_results.eVecs[0](0);
				    //ang_eigV_y = ang_results.eVecs[0](1); 
				    //ang_eigV_z = ang_results.eVecs[0](2);

				    vec_pca_ang.push_back(ang_results.eVecs[0](0));
					vec_pca_ang.push_back(ang_results.eVecs[0](1));
					vec_pca_ang.push_back(ang_results.eVecs[0](2));

			    }else{
			    	vec_pca_ang.push_back(0);
			    	vec_pca_ang.push_back(0);
			    	vec_pca_ang.push_back(0);
			    }
	    		
	        	for (int j = 0; j < pl_size; ++j){
	        		dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
	          		if (dist < old_dist){
	            		old_dist = dist;
	            		low_dist_at = j;
	          		}
	        	}
	        	TrkPoint tempPoint;
	        	tempPoint.c = points_left[low_dist_at].c;
	        	tempPoint.x = points_left[low_dist_at].x;
	        	tempPoint.y = points_left[low_dist_at].y;
	        	tempPoint.z = points_left[low_dist_at].z;
	        	tempPoint.q = points_left[low_dist_at].q;

	        	ob_temp_x = tempPoint.x - ord_trk.back().x;
	        	ob_temp_y = tempPoint.y - ord_trk.back().y;
	        	ob_temp_z = tempPoint.z - ord_trk.back().z;

	        	vec_ob_temp = Unit_Vec(ob_temp_x,ob_temp_y,ob_temp_z);

	        	double ang_dist_plus, ang_dist_minus;
	        	ang_dist_plus 	= Pythagoras(vec_ob_temp[0],vec_pca_ang[0],vec_ob_temp[1],vec_pca_ang[1],vec_ob_temp[2],vec_pca_ang[2]);
	        	ang_dist_minus 	= Pythagoras(vec_ob_temp[0],-vec_pca_ang[0],vec_ob_temp[1],-vec_pca_ang[1],vec_ob_temp[2],-vec_pca_ang[2]);

	        	if(ang_dist_plus > ang_dist_minus){
	        		ang_decide = 1.0;
	        	}else{
	        		ang_decide = -1.0;
	        	}
	        	ang_dotP = vec_ob_temp[0] * vec_pca_ang[0] + vec_ob_temp[1] * vec_pca_ang[1] + vec_ob_temp[2] * vec_pca_ang[2];
	        	ang_dotP = ang_dotP*ang_decide;

	        	//if(ang_dotP != 0) cout << ang_dotP << endl;
	        	//cout << ang_dotP << ", " << old_dist << endl;
	        	if (old_dist > alpha && ang_dotP < ang_gap){ //(ang_gap , -1]
	        		//Keep the next point
	        		ord_trk.push_back(tempPoint);
	          		if (tempPoint.y < low_ord_y){
	          			low_ord_y = tempPoint.y;
	          		}
	          		old_dist = 10000000;
	          		points_left.erase(points_left.begin() + low_dist_at);
	          		pl_size = points_left.size();
	          		i = 0;
			    }else if(old_dist > alpha){
			    	if(old_dist < low_old_dist){

			    		vec_last_ord.push_back(ord_trk.back().x);
			    		vec_last_ord.push_back(ord_trk.back().y);
			    		vec_last_ord.push_back(ord_trk.back().z);

			    		//cout << old_dist << endl;
		                last_dist = old_dist;
		                ftest_point_x = tempPoint.x;
		                ftest_point_y = tempPoint.y;
		                ftest_point_z = tempPoint.z;

		                vec_temp = vec_ob_temp;

		                //cout << vec_temp[0] << endl;

		                pca_points = ang_chunk;

		                vec_pca.push_back(ang_decide*vec_pca_ang[0]);
		                vec_pca.push_back(ang_decide*vec_pca_ang[1]);
		                vec_pca.push_back(ang_decide*vec_pca_ang[2]);

		                //pca_eVec_x = ang_eigV_x;
		                //pca_eVec_y = ang_eigV_y;
		                //pca_eVec_z = ang_eigV_z;

		                //cout << acos(vec_temp[0]*ang_decide*vec_pca_ang[0] + vec_temp[1]*ang_decide*vec_pca_ang[1]) *180./PI << endl;

		                fangle = ang_dotP;
		                low_old_dist = old_dist;
              		}
			    	points_gd.push_back(tempPoint);
			        old_dist = 10000000;
			        points_left.erase (points_left.begin() + low_dist_at);
			        pl_size = points_left.size();
			        i++;
	        	}else{
	          		ord_trk.push_back(tempPoint);
	          		if (tempPoint.y < low_ord_y){
	          			low_ord_y = tempPoint.y;
	          		}
	          		old_dist = 10000000;
	          		points_left.erase(points_left.begin() + low_dist_at);
	          		pl_size = points_left.size();
	          		i = 0;
	        	}
	        	
	        	if (pl_size == 0) break;
	    	}
			double bottom_dist;
			bottom_dist = abs(trk.back().y - low_ord_y);
			// If distance between lowest y value of unordered track and 
			// lowest y value of ordered track is greater than 10 cm.
			if(bottom_dist > 10.){
					for (int cand = 0; cand < mcand_size; ++cand){
	    				if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
	    					//cout << line << endl;
	    					ord_alg_cutouts_michels += 1;
	    					fangle = acos(fangle) * 180. / PI;
	    					//double norm_ord_back = Pythagoras(ord_trk.back().x,0.,ord_trk.back().y,0., ord_trk.back().z,0.);
	    					//double norm_eVec = Pythagoras(ord_trk.back().x + pca_eVec_x,0., ord_trk.back().y + pca_eVec_y,0., ord_trk.back().z + pca_eVec_x,0.);
	    					//double norm_test_point = Pythagoras(ftest_point_x,0., ftest_point_y,0., ftest_point_z,0.);
	    					//nt_study->Fill(run_num,ev_num,cluster, alpha, last_dist, fangle, ord_trk.back().x, ord_trk.back().y, ord_trk.back().z, ftest_point_x, ftest_point_y, ftest_point_z);
	    					nt_study->Fill(run_num,ev_num,cluster, alpha, last_dist, fangle, vec_last_ord[0], vec_last_ord[1], vec_last_ord[2], ftest_point_x, ftest_point_y, ftest_point_z);
	    					//nt_vectors->Fill(run_num, ev_num, cluster, fangle, ord_trk.back().x, ord_trk.back().y, ord_trk.back().z, vec_pca[0], vec_pca[1], vec_pca[2], vec_temp[0], vec_temp[1], vec_temp[2]);
	    					nt_vectors->Fill(run_num, ev_num, cluster, fangle, vec_last_ord[0], vec_last_ord[1], vec_last_ord[2], vec_pca[0], vec_pca[1], vec_pca[2], vec_temp[0], vec_temp[1], vec_temp[2]);
	    					//cout << run_num << ", " << ev_num << ", " << cluster << "\t Cut by Ordering Algorithm" << endl;
	    					//cout << "/Users/ivan/Work/Cosmics_Michel_Studies/Track_Plotter/Michel_Plots/Run_" << run_num << "_Event_" << ev_num << "_Cluster_" << cluster << ".pdf" << endl;
	        				for (int i = 0; i < pca_points.size(); ++i){
	        					nt_pca_points->Fill(run_num,ev_num,cluster,pca_points[i].x,pca_points[i].y,pca_points[i].z);
	        				}
	    					break;
	    				}
	    			}
				continue; //#CUT
			}else{
				for (int cand = 0; cand < mcand_size; ++cand){
	    			if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster){
	    				//cout << run_num << ", " << ev_num << ", " << cluster << "\t Cut by Ordering Algorithm" << endl;
	    				michels_survived_ord_alg += 1;
	    				break;
	    			}
	    		}
			}
			tracks_survived_ord_alg += 1;
			//Finished Ordering Points
		    //////////////////////////
		    for (int cand = 0; cand < mcand_size; ++cand){
	    		if(Michel_candidates[cand][0] == run_num && Michel_candidates[cand][1] == ev_num && Michel_candidates[cand][2] == cluster && Michel_candidates[cand][3] == -1 && Michel_candidates[cand][4] == -1 && Michel_candidates[cand][5] == -1){
	    			cout << run_num << ", " <<  ev_num << ", " << cluster << endl;
	    			cout << "/Users/ivan/Work/Cosmics_Michel_Studies/Track_Plotter/Michel_csv_files_no_ord_trk/Not_ordered_Run_" << run_num << "_Event_" << ev_num << "_Track_" << cluster << ".csv"<< endl; 
	    		}
	    	}
   		}
   		infile->Close();
	}
    //purity = ((float)michel_count)/((float)track_selected_as_michel);
    //efficiency = ((float)michel_count)/((float)(s-1));
    //ver_rms = sqrt((1./((float)michel_count)) * ver_rms);
    //ver_mean = (1./((float)michel_count)) * ver_mean;
    cout << "###################################################################" << endl;
    cout << "Total number of tracks = " << total_num_tracks << "; Number of Michels in sample = " << s - 1 << endl;
    cout << "Tracks after ordering algorithm = " << tracks_survived_ord_alg << " with alpha = " << alpha << endl;
    cout << "Michel clusters that survived the ordering algorithm = " << michels_survived_ord_alg << endl;
    cout << "Michel clusters that were cut out by ordering algorithm = " << ord_alg_cutouts_michels << endl;
    cout << "##########################################################" << endl;
    //nt_study -> Fill(alpha,angles[ang],EVRs[evr],purity,efficiency,ver_rms,ver_mean,tracks_survived_ord_alg,track_selected_as_michel,michel_count);  
  	
  	f_output->Write();
  	f_output->Close();
  	cout << "DONE" << endl;
  	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////END OF MAIN PROGRAM/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////START OF FUNCTIONS//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2){
	double dist;
	dist = sqrt(pow(x2-x1,2.) + pow(y2-y1,2.) + pow(z2-z1,2.));
	return dist;
}

vector<double> Unit_Vec(double x1,double y1,double z1){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,0.0,y1,0.0,z1,0.0);
	v.push_back(x1/norm);
	v.push_back(y1/norm);
	v.push_back(z1/norm);
	return v;
}

void LoadPointCloud(PointCloud &points, const track_def &ord_trk) {
  for (int i = 0; i < ord_trk.size(); ++i){
    Point tempPoint;
    tempPoint.x = ord_trk.at(i).x;
    tempPoint.y = ord_trk.at(i).y;
    tempPoint.z = ord_trk.at(i).z;
    tempPoint.q = ord_trk.at(i).q;
    points.push_back(tempPoint);

  }
  return;
}

PCAResults DoPCA(const PointCloud &points) {
  TVector3 outputCentroid;
  pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  vector<TVector3> outputEigenVecs;
  float meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;
  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }
  if (nThreeDHits == 0) {
    PCAResults results;
    return results; 
  }
  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);
  float xi2 = 0.0;
  float xiyi = 0.0;
  float xizi = 0.0;
  float yi2 = 0.0;
  float yizi = 0.0;
  float zi2 = 0.0;
  float weightSum = 0.0;
  for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);
      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  typedef std::pair<float,size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat(eigenMat.eigenvalues());
  eigenValColVector.emplace_back(resultEigenMat(0), 0);
  eigenValColVector.emplace_back(resultEigenMat(1), 1);
  eigenValColVector.emplace_back(resultEigenMat(2), 2);

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

  outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  for (const EigenValColPair &pair : eigenValColVector) {
     outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

  PCAResults results;

  Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

  Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
  Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

  Eigen::Vector3f testPoint;
  Eigen::Vector3f projTestPoint;
  float maxDist1 = -1.0;
  float maxDist2 = -1.0;
  float dist;
  float dotP;
  for (unsigned int i = 0; i < points.size(); i++) {
    testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
    projTestPoint = priAxis.projection(testPoint);
    dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
    dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);


    if ((dotP < 0.0) && (dist > maxDist1)) {
      endPoint1 = projTestPoint;
      maxDist1 = dist;
    }
    else if ((dotP > 0.0) && (dist > maxDist2)) {
      endPoint2 = projTestPoint;
      maxDist2 = dist;
    }
  }
  outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
  outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));
  outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));
  results.centroid = outputCentroid;
  results.endPoints = outputEndPoints;
  results.length = outputLength;
  results.eVals = outputEigenValues;
  results.eVecs = outputEigenVecs;
  return results;
}

