#include "datastructure_rob.h"
#include "io_rob.h"
#include "util_rob.h"
#include "mapper.h"

#include <string>
#include <vector>
using namespace std;


string out_file_prefix = "";
string topology_file, traffic_file, mboxSpec_file, deploy_res_file; //added by rob
string algorithm;
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
		} else if (argument.first == "--deploy_res_file") {
			deploy_res_file = argument.second;
		} else if (argument.first == "--out_file_prefix") {
			out_file_prefix = argument.second;
		}else if (argument.first == "--outPath") {
			outPath = argument.second; //path and filename prefix of output file
		}else if (argument.first == "--shareVNF") {
			shareVNF = argument.second; //"true" or "false"
		}
	}

	//cout<<"load topology data: "<<topology_file<<endl;
	InitializeTopology(topology_file.c_str());

	//cout<<"load the middlebox information: "<<mboxSpec_file<<endl;
	InitializeMiddleboxes(mboxSpec_file.c_str());
	// PrintMiddleboxes();

	//cout<<"load the traffic request info.: "<<traffic_file<<endl;
	InitializeTrafficRequests_rob(traffic_file.c_str());
	//InitializeTrafficRequests(argument.second.c_str());
	if (parallel == "true")
		InitializeTrafficRequestsSubChains_rob(traffic_file.c_str());

	cout<<"load all results"<<endl;

	if(algorithm == "cplex")
		InitializeAllResults_CPLEX_oldversion((deploy_res_file + ".sequences").c_str());
	else
		InitializeAllResults((deploy_res_file + ".sequences").c_str());


	if(algorithm == "cplex")
		InitializeCplexPaths((deploy_res_file+".paths").c_str());
}

int main(int argc, char *argv[]) {

	for (int i = 0; i < traffic_requests.size(); ++i) {
		//    traffic_requests[i].duration = 6000;
	}

	loadData(argc, argv); //parse the parameters, and load data from files
	if(out_file_prefix == "")
		out_file_prefix = deploy_res_file;

	if(algorithm == "cplex"){
		ComputeSolutionCosts_CPLEX_rob(all_mcResults, all_scResults, deploy_res_file);
		ComputeNetworkUtilization_CPLEX_rob(all_mcResults, all_scResults);
		ComputeAllStretches_CPLEX_rob(all_mcResults);
		ComputeKHops_CPLEX_rob(all_mcResults);
	}else{
		ComputeSolutionCosts_rob(all_mcResults, all_scResults);
		ComputeNetworkUtilization_rob(all_mcResults, all_scResults);
		ComputeAllStretches_rob(all_mcResults);
		ComputeKHops_rob(all_mcResults);
	}
	//ComputeServicePoints_rob(all_mcResults);
	//ComputeCloseness_rob(all_mcResults);

	ProcessCostLogs(out_file_prefix);
	ProcessStretchLogs(out_file_prefix);
	ProcessNetUtilizationLogs(out_file_prefix);
	ProcessServerUtilizationLogs(out_file_prefix);
	ProcessKHopsLogs(out_file_prefix);
	ProcessActiveServerLogs(out_file_prefix);
	//ProcessMboxRatio(out_file_prefix);
	//ProcessServicePointLogs(out_file_prefix);
	//ProcessClosenessLogs(out_file_prefix);

	cout<<"Finished!"<<endl;
	return 0;
}
