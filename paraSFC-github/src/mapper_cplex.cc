#include "datastructure_rob.h"
#include "util_rob.h"
#include "io_rob.h"
#include "mapper.h"

#include "provs_cplex.h"

#include <chrono>
#include <map>
#include <utility>
#include <memory>
#include <stdio.h>
#include <string>
#include <string.h>

double bwscale = 1; //scale up the bw capacity in the network
long testbw = -1;
string algorithm;
string topology_file, traffic_file, mboxSpec_file;  //added by rob
string parallel = "false";

//parse the parameters, and load data from files
void loadData(int argc, char *argv[]) {
	auto arg_maps = ParseArgs(argc, argv);
	for (auto argument : *arg_maps) {
		if (argument.first == "--per_core_cost") {
			per_core_cost = atof(argument.second.c_str());
		} else if (argument.first == "--per_bit_transit_cost") {
			per_bit_transit_cost = atof(argument.second.c_str());
		} else if (argument.first == "--topology_file") {
			topology_file = argument.second;
		} else if (argument.first == "--middlebox_spec_file") {
			mboxSpec_file = argument.second;
		} else if (argument.first == "--traffic_request_file") {
			traffic_file = argument.second; // added by rob
		} else if (argument.first == "--algorithm") {
			algorithm = argument.second; //viterbi or clpex
		} else if (argument.first == "--max_time") {
			max_time = atoi(argument.second.c_str()); //maximum of time line
		} else if (argument.first == "--parallel") {
			parallel = argument.second; //"true" or "false"
		} else if (argument.first == "--max_iter") {
			maxIterNum = atoi(argument.second.c_str()); //"true" or "false"
		} else if (argument.first == "--outPath") {
			outPath = argument.second; //path and filename prefix of output file
		} else if (argument.first == "--shareVNF") {
			shareVNF = argument.second; //"true" or "false"
		} else if (argument.first == "--TopoScale") {
			topoScale = atof(argument.second.c_str());
		} else if (argument.first == "--bwexpan") {
			bwfactor = atof(argument.second.c_str());
		} else if (argument.first == "--testbw") {
			testbw = atol(argument.second.c_str());
		} else if (argument.first == "--bwscale") {
			bwscale = atof(argument.second.c_str());
		}
	}


	if(parallel == "true"){
		cout<<"initialize parallelized SFCs"<<endl;
		InitializeSFCMap();
	}
	else{
		cout<<"initialize sequential SFCs"<<endl;
		InitializeSeqSFCMap();
	}

	cout << "load topology data: " << topology_file << endl;
	InitializeTopology(topology_file.c_str(), bwscale);

	cout << "load the middlebox information: " << mboxSpec_file << endl;
	InitializeMiddleboxes(mboxSpec_file.c_str());
	//PrintMiddleboxes();

	cout << "load the traffic request info.: " << traffic_file << endl;
	InitializeTrafficRequests_rob(traffic_file.c_str(), testbw);
	//InitializeTrafficRequests(argument.second.c_str());
}

int main(int argc, char *argv[]) {


	if (argc < 6) {
		puts(kUsage.c_str());
		return 1;
	}

	loadData(argc, argv);
	//return 1;

	std::vector<traffic_request> current_traffic_requests;
	std::map<int, std::vector<traffic_request> > current_traffic_requests_SC;
	int current_time = traffic_requests[0].arrival_time;
	double opex = 0.0, running_time = 0.0;
	int processed_traffic = 0;

	// files to write output
	FILE *path_file 		= fopen((outPath + ".paths").c_str(), "a");
	FILE *crtpath_file 		= fopen((outPath + ".crtpath").c_str(), "a");
	FILE *util_file 		= fopen((outPath + ".util.ts").c_str(), "a");
	FILE *time_file 		= fopen((outPath + ".time.ts").c_str(), "a");
	FILE *delay_file 		= fopen((outPath + ".delay.ts").c_str(), "a");
	FILE *sequence_file 	= fopen((outPath + ".sequences").c_str(), "a");
	FILE *sfclen_file 	    = fopen((outPath + ".sfclen.ts").c_str(), "a");
	FILE *pathlen_file 	    = fopen((outPath + ".pathlen.ts").c_str(), "a");
	FILE *copynum_file      = fopen((outPath + ".copynum.ts").c_str(), "a");
	FILE *avg_delay_file 	= fopen((outPath + ".avg.delay.ts").c_str(), "a");
	FILE *bandwidth_file 	= fopen((outPath + ".bandwidth.ts").c_str(), "a");
	FILE *avg_sfclen_file   = fopen((outPath + ".avg.sfclen.ts").c_str(), "a");
	FILE *avg_pathlen_file  = fopen((outPath + ".avg.pathlen.ts").c_str(), "a");
	FILE *avg_copynum_file  = fopen((outPath + ".avg.copynum.ts").c_str(), "a");

	string sSFCres_file = "";
	if(parallel == "true"){
		sSFCres_file = topology_file + ".cplex.bandwidth.ts";
		cout<<sSFCres_file<<endl;
	}

	// print the node and edge count at the begining of the sequence file
	// fprintf(util_file, "%d %d\n", GetNodeCount(graph), GetEdgeCount(graph));

	int time_slot = 0;
	cout << "start to cplex" << endl;
	for (int mcid = 0, cursfc = 0; mcid < traffic_requests.size() && time_slot < max_time; time_slot++) {  //control time slot
	//for (int mcid = 0, cursfc = 0; mcid < traffic_requests.size() && time_slot < max_time;) {

		//cout<<"current time: "<<current_time<<" mcid: "<<mcid<<endl;
		//cout<<"time_slot: "<<time_slot<<" mcid: "<<mcid<<endl;
		for (; mcid < traffic_requests.size() && current_time == traffic_requests[mcid].arrival_time; mcid++) {
		//for (; mcid < traffic_requests.size() && time_slot < max_time; mcid++) {
			// traffic_requests[i].duration = 6000; // 300;
			int curmc = current_traffic_requests.size();
			current_traffic_requests.push_back(traffic_requests[mcid]);
			int sid = traffic_requests[mcid].sid;

			//current_traffic_requests_SC[curmc] = pSFC[sid];

			//cout<<traffic_requests[mcid].middlebox_sequence.size()<<" ";
			for (int scid = 0; scid < traffic_requests_subchains[mcid].size(); scid++) {
				current_traffic_requests_SC[curmc].push_back(traffic_requests_subchains[mcid][scid]);
				//cout<<traffic_requests_subchains[mcid][scid].middlebox_sequence.size()<<" ";
			}
			//cout<<endl;

		}

		if(mcid < traffic_requests.size())
			current_time = traffic_requests[mcid].arrival_time;
		else
			current_time = -1;

		//in order to start from the middle
//		if(time_slot <= 99)
//		{
//			processed_traffic += current_traffic_requests.size();
//			current_traffic_requests.clear();
//			current_traffic_requests_SC.clear();
//			continue;
//		}

		std::map<int, std::vector<int> > scSequence[current_traffic_requests_SC.size()];
		std::vector<std::pair<int, int> > mcEdges[current_traffic_requests.size()];
		std::map<int, std::vector<std::pair<int, int> > > scEdges[current_traffic_requests_SC.size()];
		int delays[current_traffic_requests.size()];
		map<int, std::vector<double> > delay_breakdown;  //delay of each request
		std::vector<int> utilization;  //used cores of each node
		map<int, int> bandwidths;     //used bandwidth of each link
		map<int, int> sfclens;       //length of each sfc
		map<int, int> pathlens;     //length of each provisioned path

		cout << "call cplex" << endl;
		auto solution_start_time = std::chrono::high_resolution_clock::now();
		run_cplex(current_traffic_requests, current_traffic_requests_SC, opex,
				delay_breakdown, running_time, scSequence,
				scEdges, delays, utilization, bandwidths, topology_file, sSFCres_file, bwscale);
		auto solution_end_time = std::chrono::high_resolution_clock::now();
		unsigned long long solution_time = std::chrono::duration_cast
											< std::chrono::nanoseconds
											> (solution_end_time - solution_start_time).count();

		fprintf(time_file, "%llu.%llu\n", solution_time / ONE_GIG, solution_time % ONE_GIG);
		cout<<14<<endl;
		printf("opex: %f, solution time: %llu.%llus\n", opex, solution_time / ONE_GIG, solution_time % ONE_GIG);

		cout<<15<<endl;
		processed_traffic += current_traffic_requests.size();
		cout << processed_traffic * 100.0 / traffic_requests.size() << "% Traffic processed." << endl;

		cout << "save sequence & path of SC res to file" << endl;
		int pathlensum = 0;
		int sfclensum = 0;
		for (int rr = 0; rr < current_traffic_requests_SC.size(); rr++) {
			int scNum = current_traffic_requests_SC[rr].size();
			for(int cc = 0; cc < scNum; cc++){

				// sequence
				std::vector<int> seq = scSequence[rr][cc];
				//fprintf(sequence_file, "%d,", current_traffic_requests_SC[ii].mcID);
				for (int j = 0; j < seq.size(); ++j) {
					fprintf(sequence_file, "%d,", seq[j]);
				}
				fprintf(sequence_file, "\n");

				// path
				std::vector<std::pair<int, int>> scEdge_list = scEdges[rr][cc];
				fprintf(path_file, "%d-%d: ",rr,cc);
				for (std::pair<int, int> edge : scEdge_list) {

					fprintf(path_file, "(%d-%d),", edge.first, edge.second);
				}
				fprintf(path_file, "\n");

				if(cc == 0){
					pathlensum += scEdge_list.size();
					fprintf(pathlen_file, "%lu ", scEdge_list.size());

					sfclensum += seq.size()-2;
					fprintf(sfclen_file,  "%lu ", seq.size()-2);
				}
			}
		}
		fprintf(pathlen_file, "\n");
		fprintf(sfclen_file,  "\n");

		fprintf(avg_pathlen_file, "%f\n", 1.0*pathlensum/current_traffic_requests_SC.size());
		fprintf(avg_sfclen_file,  "%f\n", 1.0*sfclensum /current_traffic_requests_SC.size());


		cout << "delay log" << endl;
		double delaysum = 0.0;
		int copynumsum = 0;
		for (int rr = 0; rr < current_traffic_requests_SC.size(); rr++) {
			double maxdelay = -1.0;
			int maxccid = -1;
			int ccid = 0;
			for (double _delay : delay_breakdown[rr]) {
				//cout<<_delay<<" ";
				if(maxdelay < _delay)
				{
					maxdelay = _delay;
					maxccid = ccid;
				}
				ccid++;
			}
			//cout<<endl;
			fprintf(delay_file,   "%lf ", maxdelay);
			fprintf(crtpath_file, "%d ",  maxccid);

			//cout<<maxdelay<<" "<<maxccid<<endl;

			fprintf(copynum_file,   "%d ", copyNum[rr]);

			delaysum += maxdelay;
			copynumsum += copyNum[rr];
		}
		fprintf(delay_file,   "\n");
		fprintf(crtpath_file, "\n");
		fprintf(copynum_file, "\n");

		fprintf(avg_delay_file,   "%f\n", 1.0*delaysum/current_traffic_requests_SC.size());
		fprintf(avg_copynum_file, "%f\n", 1.0*copynumsum/current_traffic_requests_SC.size());

		cout << "utilization log" << endl;
		for (int cores : utilization) {
			fprintf(util_file, "%d\n", cores);
		}

		cout << "bandwidth log " << bandwidths.size()<< endl;
		for (int r = 0; r < bandwidths.size(); r++) {
			fprintf(bandwidth_file, "%d\n", bandwidths[r]);
		}
		if(parallel=="false"){
			FILE *each_seqbandwidth = fopen((topology_file + ".cplex.bandwidth.ts").c_str(), "w");
			for (int r = 0; r < bandwidths.size(); r++) {
				fprintf(each_seqbandwidth, "%d\n", bandwidths[r]);
			}
			fclose(each_seqbandwidth);
		}

//		fflush(sequence_file);
//		fflush(path_file);
//		fflush(delay_file);
//		fflush(util_file);

		current_traffic_requests.clear();
		current_traffic_requests_SC.clear();
		//exit(0);

		//break;// if only one time slot
	}
	//end of for each time_slot

	// close all the output files
	fclose(sequence_file);
	fclose(path_file);
	fclose(delay_file);
	fclose(util_file);
	fclose(time_file);
	fclose(avg_delay_file);
	fclose(bandwidth_file);
	fclose(pathlen_file);
	fclose(avg_pathlen_file);
	fclose(sfclen_file);
	fclose(avg_sfclen_file);
	fclose(copynum_file);
	fclose(avg_copynum_file);

	return 0;
}
