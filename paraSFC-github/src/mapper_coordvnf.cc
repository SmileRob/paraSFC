#include "datastructure_rob.h"
#include "util_rob.h"
#include "io_rob.h"
#include "provs_coordvnf.h"
#include "mapper.h"

#include <chrono>
#include <map>
#include <utility>
#include <memory>
#include <stdio.h>
#include <string>
#include <string.h>
#include <iostream>

using namespace std;

map<int, unsigned long long> solutiontimes; //var[iter] = solution time
string algorithm;
map<int, int> nodereserve;  //var[nodeid] = amount of reservation

void toProvision() {
	int current_time = traffic_requests[0].arrival_time;
	unsigned long long elapsed_time = 0;
	unsigned long long current_solution_time = 0;
	unsigned long long BWconsumption = 0;
	double totaldelay = 0.0;
	stats.num_accepted = stats.num_rejected = 0;
	const int kNumTrafficRequests = static_cast<int>(traffic_requests.size());
	const int kNumSubchains = static_cast<int>(traffic_requests_subchains.size());

	int time_slot = 0;
	//cout << "start to provision" << endl;
	for (int mcid = 0; mcid < traffic_requests.size() && time_slot < max_time; time_slot++) {    //!!!!!!!!!!!!!!!control time slot
	//for (int mcid = 0; mcid < traffic_requests.size() && mcid < max_time;) {

		unsigned long initBW = GetTotalResidualBandwidth_Rob();

		//cout<<"current time: "<<current_time<<" mcid: "<<mcid<<endl;
		//cout<<"time_slot: "<<time_slot<<" mcid: "<<mcid<<endl;
		for (; mcid < traffic_requests.size() && current_time == traffic_requests[mcid].arrival_time; mcid++) {    //!!!!!!!!!!!!!!!control time slot
		//for (; mcid < traffic_requests.size() && mcid < max_time; mcid++) {
			int currmcid = mcid;

			vector<int> ssfc = traffic_requests[currmcid].middlebox_sequence;
			int sid = traffic_requests[currmcid].sid;
			all_scResults[currmcid].clear();

			int scNum = traffic_requests_subchains[currmcid].size();
			int crtID = NIL;
			double crtdelay = -INF;
			bool hasprov[100] = {false};
			int provnode[100] = {-1};
			reservedEdges.clear();///
			reservedCores.clear();///

			for(int cc = 0; cc < scNum; cc++){

				rejectflag = 0;
				current_solution_time = 0;
				double sfcdelay=INF;

				vector<int> tmpvec;
				all_scResults[currmcid].push_back(tmpvec);
				// Get solution for one traffic by calling the viterbi algorithm.

				auto solution_start_time = std::chrono::high_resolution_clock::now();
				// call ******viterbi algorithm****** to workout the solution
				std::unique_ptr<std::vector<int>> result =
						coordProvision(traffic_requests[currmcid], cc, sfcdelay, hasprov, provnode);
				auto solution_end_time = std::chrono::high_resolution_clock::now();

				if (rejectflag == 1) {
					traffic_requests[currmcid].accepted = false;
					//cout << "[!!!]  reject one request: " << currmcid <<"-"<< cc << endl;
				}
				else{
					traffic_requests_subchains[currmcid][cc].SFCdelay = sfcdelay;
					//cout<<"sfcdelay: "<<sfcdelay<<endl;
					if(sfcdelay > crtdelay){
						crtID = cc;
						crtdelay = sfcdelay;
					}
				}

				//cout << "calculate the running time for this req" << endl;
				unsigned long long solution_time = std::chrono::duration_cast
						< std::chrono::nanoseconds > (solution_end_time - solution_start_time).count();

				//cout << "update time info." << endl;
				current_solution_time += solution_time;	// sum up the running time
				elapsed_time += solution_time;// update the elapsed time, which means time is always going

				//cout << "Save this request's provision result." << endl;
				//cout<<"all results size: "<<all_scResults[mcid].size()<<endl;
				vector<int> tmpsc = traffic_requests_subchains[currmcid][cc].middlebox_sequence;
				for (int tmpi = 0; tmpi < result->size(); tmpi++) {
					all_scResults[currmcid][cc].push_back(result->at(tmpi));

					if(tmpi<result->size()-1){
						if(reservedEdges.find(make_pair(result->at(tmpi),result->at(tmpi+1))) == reservedEdges.end())
						{
							reservedEdges[make_pair(result->at(tmpi),result->at(tmpi+1))] = traffic_requests_subchains[currmcid][cc].min_bandwidth;
						}
					}

					if(!hasprov[tmpsc[tmpi]]){
						reservedCores[result->at(tmpi)] += middleboxes[ssfc[tmpsc[tmpi]]].cpu_requirement;
					}

					if(tmpsc[tmpi] != 0 && tmpsc[tmpi] != ssfc.size()-1){
						hasprov[tmpsc[tmpi]] = true;
						provnode[tmpsc[tmpi]] = result->at(tmpi);
					}
					//cout<<"mcid & cc: "<<mcid<<" "<<cc<<" : "<<result->at(tmpi)<<" ";
					//cout<<result->at(tmpi)<<" ";
				}
				//cout<<endl;
				result->clear();

			}

			if (traffic_requests[currmcid].accepted == true) {
				stats.num_accepted++;
				crtpaths[currmcid] = make_pair(crtID, crtdelay);
				totaldelay += crtdelay;
				traffic_requests[currmcid].SFCdelay = crtdelay;
				//cout<<"crt delay: "<<crtdelay<<endl;
				bool updatedVNF[50]={false};
				bool updatedEdges[50][50]={false};
				for(int cc = 0; cc < scNum; cc++){
					//UpdateResources_Rob(&all_scResults[currmcid][cc], traffic_requests_subchains[currmcid][cc], ssfc);
					//continue;

					//cout<<all_scResults[currmcid][cc].size()<<" "<<traffic_requests_subchains[currmcid][cc].middlebox_sequence.size()<<endl;
					vector<int> tmpresult = all_scResults[currmcid][cc];
					vector<int> tmpseq = traffic_requests_subchains[currmcid][cc].middlebox_sequence;
					for (int tmpi = 0; tmpi < tmpresult.size(); tmpi++) {
						//if(tmpi != tmpresult.size()-1 && !updatedEdges[tmpseq[tmpi]][tmpseq[tmpi+1]]){
						if(tmpi != tmpresult.size()-1 && !updatedEdges[tmpseq[tmpi]][tmpseq[tmpi+1]]){
							ReducePathResidualBandwidth(tmpresult[tmpi], tmpresult[tmpi+1], traffic_requests_subchains[currmcid][cc].min_bandwidth);
							//updatedEdges[tmpseq[tmpi]][tmpseq[tmpi+1]]=true;
							updatedEdges[tmpresult[tmpi]][tmpresult[tmpi+1]]=true;

						}

						if(tmpi !=0 && !updatedVNF[tmpseq[tmpi]]){
							int vnfID = ssfc[tmpseq[tmpi]];
							UpdateMiddleboxInstances(tmpresult[tmpi], &middleboxes[vnfID], traffic_requests_subchains[currmcid][cc]);
							updatedVNF[tmpseq[tmpi]]=true;
						}
					}
				}

				//for(int cc = 0; cc < scNum; cc++){
				//	UpdateResources(&all_scResults[currmcid][cc], traffic_requests_subchains[currmcid][cc]);
				//}
			} else if (traffic_requests[currmcid].accepted == false) {
				stats.num_rejected++;
				traffic_requests[currmcid].SFCdelay = INF;
				//cout << "[!!!] Reject request " << currmcid <<"(sid: "<<sid<<")" << endl;
			}

			//break;   //for testing!!!!!!!!!!!!!!!!!!

		}

		//cout<<"timeslot & elapsed time & totaldelay: "<<time_slot<<" "<<elapsed_time<<" "<<totaldelay<<endl;
		solutiontimes[time_slot] = elapsed_time;

		unsigned long deltaBW = initBW - GetTotalResidualBandwidth_Rob();
		BWconsumption += deltaBW;
		//cout<<"Used BW: "<<deltaBW<<endl;
		GetTotalNodeResidualCores_Rob();
		//cout<<"The average totaldelay: "<<totaldelay*1.0/double(stats.num_accepted)<<endl;

		ReleaseAllResources();

		if(mcid < traffic_requests.size())
			current_time = traffic_requests[mcid].arrival_time;
		else
			current_time = -1;

	}// end for a time slot processing

	// cout<<"Print the solution time for the final round."<<endl;
	// printf("Current time = %d, Solution time = %llu.%llu\n", current_time,
	//		current_solution_time / ONE_GIG, current_solution_time % ONE_GIG);

	// Print the solution time.
	//printf("Solution time: %llu.%llus\n", elapsed_time / ONE_GIG, elapsed_time % ONE_GIG);

	double acceptance = double(stats.num_accepted) / double(stats.num_accepted + stats.num_rejected);
	//printf("Acceptance Ratio: %.8lf%%\n", 100*acceptance);
	FILE *accept_file = fopen((outPath+".acceptance").c_str(), "a");
	fprintf(accept_file, "%f\n", acceptance);
	fclose(accept_file);

	FILE *avrg_latency_file = fopen((outPath+".average_latency").c_str(), "a");
	fprintf(avrg_latency_file, "%f\n", totaldelay*1.0/double(stats.num_accepted));
	fclose(avrg_latency_file);

	//cout << "Write all the SC provision results in a file." << endl;
	FILE *sequence_file = fopen((outPath+".sequences").c_str(), "w");
	FILE *latency_file = fopen((outPath+".latency").c_str(), "w");
	for (int mcid = 0; mcid < all_scResults.size(); mcid++) {
		fprintf(latency_file, "%f\n", traffic_requests[mcid].SFCdelay);

		vector<vector<int> > scrow = all_scResults[mcid];
		for (int tmpi = 0; tmpi < scrow.size(); tmpi++) {
			fprintf(sequence_file, "%d,", mcid);
			for (int tmpjj = 0; tmpjj < scrow[tmpi].size(); tmpjj++) {
				fprintf(sequence_file, "%d,", scrow[tmpi][tmpjj]);
			}
			fprintf(sequence_file, "\n");
		}
		if (scrow.size() == 0)
			fprintf(sequence_file, "%d,\n", mcid);
	}
	fclose(sequence_file);
	fclose(latency_file);

	//cout << "Write solution time to file." << endl;
	//cout<<"solutiontimes.size(): "<<solutiontimes.size()<<endl;
	double lastVal = 0.0;
	FILE *soltionTime_file = fopen((outPath+".solutiontime").c_str(), "w");
	for (int tSlot = 0; tSlot < solutiontimes.size(); tSlot++) {
		unsigned long long tmp = solutiontimes[tSlot]-lastVal;
		lastVal = solutiontimes[tSlot];
		fprintf(soltionTime_file, "%llu.%llu\n", tmp / ONE_GIG, tmp % ONE_GIG);
	}
	fclose(soltionTime_file);

	FILE *total_solutiontime = fopen((outPath+".total_solutiontime").c_str(), "a");
	fprintf(total_solutiontime, "%llu.%llu\n", elapsed_time / ONE_GIG, elapsed_time % ONE_GIG);
	fclose(total_solutiontime);

	FILE *total_bwconsumption = fopen((outPath+".total_bwconsumption").c_str(), "a");
	fprintf(total_bwconsumption, "%llu\n", BWconsumption);
	fclose(total_bwconsumption);
}

//parse the parameters, and load data from files
void loadData(int argc, char *argv[]) {
	auto arg_maps = ParseArgs(argc, argv);
	string topology_file, traffic_file, mboxSpec_file;  //added by rob
	string parallel = "true";
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
		}
	}
	if(parallel == "true"){
		//cout<<"initialize parallelized SFCs"<<endl;
		InitializeSFCMap();
	}
	else{
		//cout<<"initialize sequential SFCs"<<endl;
		InitializeSeqSFCMap();
	}

	//cout<<"load topology data: "<<topology_file<<endl;
	InitializeTopology(topology_file.c_str());

	//cout<<"load the middlebox information: "<<mboxSpec_file<<endl;
	InitializeMiddleboxes(mboxSpec_file.c_str());
	// PrintMiddleboxes();

	//cout<<"load the traffic request info.: "<<traffic_file<<endl;
	InitializeTrafficRequests_rob(traffic_file.c_str());

}

int main(int argc, char *argv[]) {

	if (argc < 6) {
		puts(kUsage.c_str());
		return 1;
	}

	loadData(argc, argv); //parse the parameters, and load data from files

	shareVNF = "false";

	toProvision();

	//int a;  std::cin>>a;

	return 0;
}
