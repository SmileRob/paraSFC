#include "datastructure_rob.h"
#include "util_rob.h"
#include "io_rob.h"
#include "provs_parasfc.h"
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

double bwscale=1.0;
long testbw = -1;
string parallel = "true";
string topology_file, traffic_file, mboxSpec_file;  //added by rob
string algorithm;
map<int, unsigned long long> solutiontimes; //var[iter] = solution time
map<int, int> nodereserve;  //var[nodeid] = amount of reservation

void preProvision(){
	int current_time = traffic_requests[0].arrival_time;
	int time_slot = 0;
	//cout << "start to pre-provision" << endl;
	for (int mcid = 0; mcid < traffic_requests.size() && time_slot < max_time; time_slot++) { //control time slot
	//for (int mcid = 0; mcid < traffic_requests.size() && mcid < max_time; ) {  //control mcid

		//c out<<"time_slot: "<<time_slot<<" mcid: "<<mcid<<endl;

		for (; mcid < traffic_requests.size() && current_time == traffic_requests[mcid].arrival_time; mcid++) {    //!!!!!!!!!!!!!!!control time slot
		//for (; mcid < traffic_requests.size() && mcid < max_time; mcid++) {

			int scNum = traffic_requests_subchains[mcid].size();
			int sid = traffic_requests[mcid].sid;
			all_scResults[mcid].clear();

			for(int cc = 0; cc < scNum; cc++){

				rejectflag = 0;
				double sfcdelay=INF;
				vector<int> tmpvec;
				all_scResults[mcid].push_back(tmpvec);
				std::unique_ptr<std::vector<int>> result = ViterbiEstimate(traffic_requests[mcid], cc, sfcdelay);

				if (rejectflag == 1) {
					traffic_requests[mcid].accepted = false;
					//cout << "[!!!]  reject one request: " << mcid <<"-"<< cc << endl;
				}

				//cout << "Save this request's provision result." << endl;
				//c out<<"all results size: "<<all_scResults[mcid].size()<<endl;
				for (int tmpi = 0; tmpi < result->size(); tmpi++) {
					all_scResults[mcid][cc].push_back(result->at(tmpi));
					//c out<<"mcid & cc: "<<mcid<<" "<<cc<<" : "<<result->at(tmpi)<<" ";
					//c out<<result->at(tmpi)<<" ";
				}
				//c out<<endl;
				result->clear();
			}

			if (traffic_requests[mcid].accepted == true) {
				for (int cc = 0; cc < scNum; cc++){
					for (int tmpi = 0; tmpi < all_scResults[mcid][cc].size(); tmpi++) {
						nodereserve[all_scResults[mcid][cc][tmpi]]++;
					}
				}
			} else if (traffic_requests[mcid].accepted == false) {
				//cout << "[!!!] Reject request " << mcid <<"(sid: "<<sid<<")" << endl;
			}
			//break;   //for testing!!!!!!!!!!!!!!!!!!
		}

		//c out<<mcid<<": "<<traffic_requests[mcid].min_bandwidth<<" "<<traffic_requests[mcid].sid<<endl;

		if(mcid < traffic_requests.size())
			current_time = traffic_requests[mcid].arrival_time;
		else
			current_time = -1;

	}// end for a time slot processing
}

bool sortFun(const pair<int, int> &a, const pair<int, int> &b)
{
	return a.second < b.second;// ascending sorting
	//return a.second > b.second;// descending sorting
}

void adjustPriority(){
	vector<pair<int, int> > abusiveness;
	int current_time = traffic_requests[0].arrival_time;
	int time_slot = 0;
	int mcid = 0;
	//cout << "start to adjust priority" << endl;
	for (time_slot = 0; time_slot < max_time; time_slot++) {    //!!!!!!!!!!!!!!!control time slot
	//for (int mcid = 0; mcid < traffic_requests.size() && mcid < max_time; ) {

		//c out<<"time_slot: "<<time_slot<<" mcid: "<<mcid<<endl;
		int lowID = mcid, highID = mcid;
		abusiveness.clear();
		for (; mcid < traffic_requests.size() && current_time == traffic_requests[mcid].arrival_time; mcid++,highID++) {    //!!!!!!!!!!!!!!!control time slot
		//for (; mcid < traffic_requests.size() && mcid < max_time; mcid++,highID++) {

			pair<int, int> tmpabs;
			tmpabs.first = mcid;
			tmpabs.second = 0;
			int scNum = traffic_requests_subchains[mcid].size();
			int sid = traffic_requests[mcid].sid;
			if (traffic_requests[mcid].accepted == true) {
				for(int cc = 0; cc < scNum; cc++){
					for (int tmpi = 0; tmpi < all_scResults[mcid][cc].size(); tmpi++) {
						tmpabs.second += nodereserve[all_scResults[mcid][cc][tmpi]];
					}
				}
			} else if (traffic_requests[mcid].accepted == false) {
				tmpabs.second = INF;
				traffic_requests[mcid].accepted = true;
				//cout << "[!!!] Reject request " << mcid <<"(sid: "<<sid<<")" << endl;
			}
			abusiveness.push_back(tmpabs);
			//break;   //for testing!!!!!!!!!!!!!!!!!!
		}

		sort(abusiveness.begin(), abusiveness.end(), sortFun);
		for(int tmpi = 0; tmpi < abusiveness.size(); tmpi++){
			mcOrder[lowID+tmpi] = abusiveness[tmpi].first;
			//c out<<lowID+tmpi<<" "<<abusiveness[tmpi].first<<endl;
		}


		if(mcid < traffic_requests.size())
			current_time = traffic_requests[mcid].arrival_time;
		else
			current_time = -1;

	}// end for a time slot processing
}

void toProvision(map<int, unsigned long> ssfcbwmap, double bw_factor) {
	int current_time = traffic_requests[0].arrival_time;
	unsigned long long elapsed_time = 0;
	unsigned long long current_solution_time = 0;
	unsigned long long BWconsumption = 0;
	double totaldelay = 0.0;
	stats.num_accepted = stats.num_rejected = 0;
	const int kNumTrafficRequests = static_cast<int>(traffic_requests.size());
	const int kNumSubchains = static_cast<int>(traffic_requests_subchains.size());
	map<int, unsigned long> bwmap;
	int  notexpansum = 0;

	int time_slot = 0;
	int mcid = 0;
	//cout << "start to provision" << endl;
	for (time_slot = 0; time_slot < max_time; time_slot++) {     //!!!!!!!!!!!!!!!control time slot
	//for (int mcid = 0; mcid < traffic_requests.size() && mcid < max_time; ) {

		unsigned long initBW = GetTotalResidualBandwidth_Rob();
		//c out<<"current time: "<<current_time<<" mcid: "<<mcid<<endl;
		//c out<<"time_slot: "<<time_slot<<" mcid: "<<mcid<<endl;
		for (; mcid < traffic_requests.size() && current_time == traffic_requests[mcid].arrival_time; mcid++) {    //!!!!!!!!!!!!!!!control time slot
		//for (; mcid < traffic_requests.size() && mcid < max_time; mcid++) {

			int currmcid = mcOrder[mcid];
			//cout<<"currmcid of mcid: "<<currmcid<<" "<<mcid<<endl;
			//c out<<currmcid<<": "<<traffic_requests[currmcid].sid<<" "<<traffic_requests[currmcid].min_bandwidth<<endl;
			//c out<<"+++++++++++++++++++++++++++++++++++++++"<<endl;
			vector<int> ssfc = traffic_requests[currmcid].middlebox_sequence;
			int sid = traffic_requests[currmcid].sid;
			all_scResults[currmcid].clear();

			int scNum = traffic_requests_subchains[currmcid].size();
			int crtID = NIL;
			double crtdelay = -INF;
			bool hasprov[100] = {false};
			int provnode[100] = {-1};
			reservedEdges.clear();
			reservedCores.clear();

			bool countEdges[100][100]={false};
			int  bwsum_r=0;
			int  tmpbw=traffic_requests[currmcid].min_bandwidth;

			for(int cc = 0; cc < scNum; cc++){
				//c out<<"scNum and cc: "<<scNum<<" "<<cc<<endl;

				rejectflag = 0;
				current_solution_time = 0;
				double sfcdelay=INF;

				vector<int> tmpvec;
				all_scResults[currmcid].push_back(tmpvec);
				// Get solution for one traffic by calling the viterbi algorithm.

				auto solution_start_time = std::chrono::high_resolution_clock::now();
				// call ******viterbi algorithm****** to workout the solution
				std::unique_ptr<std::vector<int>> result =
						ViterbiProvision(traffic_requests[currmcid], cc, sfcdelay, hasprov, provnode, 0);
				auto solution_end_time = std::chrono::high_resolution_clock::now();

				if (rejectflag == 1) {
					traffic_requests[currmcid].accepted = false;
					cout << "[!!!]  reject one request: " << currmcid <<"-"<< cc << endl;
				}
				else{
					traffic_requests_subchains[currmcid][cc].SFCdelay = sfcdelay;
					//if(currmcid==6 || currmcid==63)cout<<"cc and sfcdelay: "<<cc<<" "<<sfcdelay<<endl;
					if(sfcdelay > crtdelay){
						crtID = cc;
						crtdelay = sfcdelay;
					}
					//if(currmcid==0)cout<<"sfcdealy & crtdelay: "<<sfcdelay<<" "<<crtdelay<<endl;
				}

				//cout << "calculate the running time for this req" << endl;
				unsigned long long solution_time = std::chrono::duration_cast
						< std::chrono::nanoseconds > (solution_end_time - solution_start_time).count();

				//cout << "update time info." << endl;
				current_solution_time += solution_time;	// sum up the running time
				elapsed_time += solution_time;// update the elapsed time, which means time is always going

				//cout << "Save this request's provision result." << endl;
				//c out<<"all results size: "<<all_scResults[mcid].size()<<endl;
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
						//hasprov[tmpi] = true;
						//provnode[tmpi] = result->at(tmpi);
						hasprov[tmpsc[tmpi]] = true;
						provnode[tmpsc[tmpi]] = result->at(tmpi);////?????????
					}

					//c out<<"mcid & cc: "<<mcid<<" "<<cc<<" : "<<result->at(tmpi)<<" ";
					//c out<<result->at(tmpi)<<" ";

					double resbw = 0.088;
					bool addbw = false;
					if(parallel == "true" && mpmap[sid].find( make_pair(tmpsc[tmpi], tmpsc[tmpi+1])) != mpmap[sid].end() ){
						addbw = true;
					}

					//c out<<"tmpi and result size: "<<tmpi<<" "<<result->size()<<endl;
					if(tmpi < result->size()-1 &&
							!countEdges[result->at(tmpi)][result->at(tmpi+1)] &&
							result->at(tmpi)!=result->at(tmpi+1)){ //link merged
					//if(tmpi < result->size()-1 && !countEdges[tmpsc[tmpi]][tmpsc[tmpi+1]]){ // trunk merge
						//c out<<result->at(tmpi)<<" "<<result->at(tmpi+1)<<endl;;
						int tmplen=ComputeShortestPath(result->at(tmpi),result->at(tmpi+1))->size()-1;
						bwsum_r += tmplen*tmpbw;
						//if(currmcid == 53)
						//cout<<"tmplen and bwsumAAA: "<<result->at(tmpi)<<" -- "<<result->at(tmpi+1)<<" "<<tmplen<<" "<<tmpbw<<endl;
						countEdges[result->at(tmpi)][result->at(tmpi+1)]=true;  //link merged
						//countEdges[tmpsc[tmpi]][tmpsc[tmpi+1]]=true;  //trunk merged
					} else if(tmpi < result->size()-1 &&
							addbw  && result->at(tmpi)!=result->at(tmpi+1)){
						int tmplen=ComputeShortestPath(result->at(tmpi),result->at(tmpi+1))->size()-1;
						bwsum_r += tmplen*tmpbw*resbw;
						//cout<<"tmplen and bwsumBBB: "<<result->at(tmpi)<<" -- "<<result->at(tmpi+1)<<" "<<tmplen<<" "<<tmpbw<<endl;
					}
					//c out<<"bwsum_r: "<<bwsum_r<<endl;

				}
				//c out<<endl;
				result->clear();

			}



			//cout<<"for p-SFCs: check if the bw resources are bounded."<<endl;
			if (parallel=="true" && traffic_requests[currmcid].accepted == true) {
				int bw_req = ssfcbwmap[currmcid];
				//c out<<currmcid<<" "<<traffic_requests[currmcid].sid<<"  "<<traffic_requests[currmcid].min_bandwidth<<" "<<bwsum_r<<" "<<bw_req<<endl;

				//cout<<"+++++++++++++++++++"<<bwsum_r<<" "<<bw_req<<" "<<bw_factor<<endl;
				if(bwsum_r-bw_req > bw_req*bw_factor){

					notexpansum++;
					//c out<<"mcid and sid and ssfc.size(): "<<currmcid<<" "<<sid<<" "<<ssfc.size()<<endl;
					//c out<<"bw expansion too much: "<<bwsum_r-bw_req<<" "<<bw_req*bw_factor<<" ";
					//c out<<" << bwsum, bwreq and bwexpan: "<<bwsum_r<<" "<<bw_req<<" "<<bw_factor<<endl;
					int cc=0;

					traffic_request tmpreq = traffic_requests_subchains[currmcid][cc];
					traffic_requests_subchains[currmcid].clear();
					traffic_requests_subchains[currmcid].push_back(tmpreq);
					//c out<<"traffic_requests_subchains[currmcid].size(): "<<traffic_requests_subchains[currmcid].size()<<endl;
					vector<int> tmpseq;
					for(int tmpi=0; tmpi< ssfc.size(); tmpi++){
						tmpseq.push_back(tmpi);
					}
					traffic_requests_subchains[currmcid][cc].middlebox_sequence = tmpseq;

					//c out<<"crtdelay before: "<<crtdelay<<endl;
					crtID = NIL;
					crtdelay = -INF;
					bool seqhasprov[100] = {false};
					int seqprovnode[100] = {-1};
					reservedEdges.clear();
					reservedCores.clear();

					rejectflag = 0;
					current_solution_time = 0;
					double sfcdelay=INF;
					bwsum_r=0;

					all_scResults[currmcid].clear();
					vector<int> tmpvec;
					all_scResults[currmcid].push_back(tmpvec);
					// Get solution for one traffic by calling the viterbi algorithm.

					auto solution_start_time = std::chrono::high_resolution_clock::now();
					// call ******viterbi algorithm****** to workout the solution
					std::unique_ptr<std::vector<int>> result =
							ViterbiProvision(traffic_requests[currmcid], cc, sfcdelay, seqhasprov, seqprovnode, ssfc.size());
					auto solution_end_time = std::chrono::high_resolution_clock::now();

					if (rejectflag == 1) {
						traffic_requests[currmcid].accepted = false;
						//cout << "[!!!]  reject one request: " << currmcid <<"-"<< cc << endl;
					}
					else{
						traffic_requests_subchains[currmcid][cc].SFCdelay = sfcdelay;
						crtID = 0;
						crtdelay = sfcdelay;
					}


					//c out<<"crtdelay after: "<<crtdelay<<endl;

					//cout << "calculate the running time for this req" << endl;
					unsigned long long solution_time = std::chrono::duration_cast
							< std::chrono::nanoseconds > (solution_end_time - solution_start_time).count();

					//cout << "update time info." << endl;
					current_solution_time += solution_time;	// sum up the running time
					elapsed_time += solution_time;// update the elapsed time, which means time is always going

					bool seqcountEdges[100][100]={false};

					//cout << "Save this request's provision result." << endl;
					//c out<<"all results size: "<<all_scResults[mcid].size()<<endl;
					vector<int> tmpsc = traffic_requests_subchains[currmcid][cc].middlebox_sequence;
					for (int tmpi = 0; tmpi < result->size(); tmpi++) {
						all_scResults[currmcid][cc].push_back(result->at(tmpi));

						if(tmpi<result->size()-1){
							if(reservedEdges.find(make_pair(result->at(tmpi),result->at(tmpi+1))) == reservedEdges.end())
							{
								reservedEdges[make_pair(result->at(tmpi),result->at(tmpi+1))] = traffic_requests_subchains[currmcid][cc].min_bandwidth;
							}
						}

						if(!seqhasprov[tmpsc[tmpi]]){
							reservedCores[result->at(tmpi)] += middleboxes[ssfc[tmpsc[tmpi]]].cpu_requirement;
						}

						if(tmpsc[tmpi] != 0 && tmpsc[tmpi] != ssfc.size()-1){
							//hasprov[tmpi] = true;
							//provnode[tmpi] = result->at(tmpi);
							seqhasprov[tmpsc[tmpi]] = true;
							seqprovnode[tmpsc[tmpi]] = result->at(tmpi);////?????????
						}

						bool addbw = false;
						double resbw = 0.088;
						if(parallel == "true" && mpmap[sid].find( make_pair(tmpsc[tmpi], tmpsc[tmpi+1])) != mpmap[sid].end() ){
							addbw = true;
						}

						//c out<<"tmpi and result size: "<<tmpi<<" "<<result->size()<<endl;
						if(tmpi < result->size()-1 &&
								!seqcountEdges[result->at(tmpi)][result->at(tmpi+1)] &&
								result->at(tmpi)!=result->at(tmpi+1)){ //link merged
						//if(tmpi < result->size()-1 && !seqcountEdges[tmpsc[tmpi]][tmpsc[tmpi+1]]){ // trunk merge
							//c out<<"mcid & cc: "<<mcid<<" "<<cc<<" : "<<result->at(tmpi)<<" ";
							//c out<<result->at(tmpi)<<" ";
							int tmplen=ComputeShortestPath(result->at(tmpi),result->at(tmpi+1))->size();
							bwsum_r += tmplen*tmpbw;
							seqcountEdges[result->at(tmpi)][result->at(tmpi+1)]=true; //link merged
							//seqcountEdges[tmpsc[tmpi]][tmpsc[tmpi+1]]=true; //trunk merged
						} else if(tmpi < result->size()-1 && addbw && result->at(tmpi)!=result->at(tmpi+1)){
							int tmplen=ComputeShortestPath(result->at(tmpi),result->at(tmpi+1))->size();
							bwsum_r += tmplen*tmpbw*resbw;
						}

					}
					//c out<<endl;
					result->clear();
				}
			}


			//cout<<"update resources: ReducePathResidualBandwidth etc."<<endl;
			if (traffic_requests[currmcid].accepted == true) {
				stats.num_accepted++;
				crtpaths[currmcid] = make_pair(crtID, crtdelay);
				totaldelay += crtdelay;
				traffic_requests[currmcid].SFCdelay = crtdelay;
				//c out<<"crtdelay: "<<crtdelay<<endl;
				//c out<<currmcid<<": "<<traffic_requests[currmcid].SFCdelay<<endl;
				bool updatedVNF[50]={false};
				bool updatedEdges[100][100]={false};
				scNum = traffic_requests_subchains[currmcid].size();
				for(int cc = 0; cc < scNum; cc++){
					//UpdateResources_Rob(&all_scResults[currmcid][cc], traffic_requests_subchains[currmcid][cc], ssfc);
					//continue;

					//c out<<"scNum and cc: "<<scNum<<" "<<cc<<endl;
					//c out<<all_scResults[currmcid][cc].size()<<" "<<traffic_requests_subchains[currmcid][cc].middlebox_sequence.size()<<endl;
					vector<int> tmpresult = all_scResults[currmcid][cc];
					vector<int> tmpseq = traffic_requests_subchains[currmcid][cc].middlebox_sequence;
					for (int tmpi = 0; tmpi < tmpresult.size(); tmpi++) {

						bool addbw = false;
						double resbw = 0.088;
						if(parallel == "true" && mpmap[sid].find( make_pair(tmpseq[tmpi], tmpseq[tmpi+1])) != mpmap[sid].end() ){
							addbw = true;
						}

						//if(tmpi != tmpresult.size()-1 && !updatedEdges[tmpseq[tmpi]][tmpseq[tmpi+1]]){
						if(tmpi != tmpresult.size()-1 && !updatedEdges[tmpresult[tmpi]][tmpresult[tmpi+1]]){
							ReducePathResidualBandwidth(tmpresult[tmpi], tmpresult[tmpi+1], traffic_requests_subchains[currmcid][cc].min_bandwidth);
							//updatedEdges[tmpseq[tmpi]][tmpseq[tmpi+1]]=true;
							updatedEdges[tmpresult[tmpi]][tmpresult[tmpi+1]]=true;
						}
						else if(tmpi < tmpresult.size()-1 && addbw && tmpresult[tmpi] != tmpresult[tmpi+1]){
							ReducePathResidualBandwidth(tmpresult[tmpi], tmpresult[tmpi+1],
									resbw*traffic_requests_subchains[currmcid][cc].min_bandwidth);
						}

						if(tmpi !=0 && !updatedVNF[tmpseq[tmpi]]){
							int vnfID = ssfc[tmpseq[tmpi]];
							UpdateMiddleboxInstances(tmpresult[tmpi], &middleboxes[vnfID], traffic_requests_subchains[currmcid][cc]);
							updatedVNF[tmpseq[tmpi]]=true;
						}
					}
				}
			} else if (traffic_requests[currmcid].accepted == false) {
				stats.num_rejected++;
				traffic_requests[currmcid].SFCdelay = INF;
				//cout << "[!!!] Reject request " << currmcid <<"(sid: "<<sid<<")" << endl;
			}
			BWconsumption += bwsum_r;
			bwmap[currmcid] = bwsum_r;

			//c out<<currmcid<<": "<<traffic_requests[currmcid].sid<<" "<<traffic_requests[currmcid].min_bandwidth<<" "<<bwsum_r<<endl;

			//cout<<"going to next mcid"<<endl;
			//break;   //for testing!!!!!!!!!!!!!!!!!!

		}

		//cout<<"timeslot & elapsed time & totaldelay: "<<time_slot<<" "<<elapsed_time<<" "<<totaldelay<<endl;
		solutiontimes[time_slot] = elapsed_time;


		//unsigned long deltaBW = initBW - GetTotalResidualBandwidth_Rob();
		//BWconsumption += deltaBW;

		//c out<<"Used BW: "<<deltaBW<<endl;
		GetTotalNodeResidualCores_Rob();
		//c out<<"The average totaldelay: "<<totaldelay*1.0/double(stats.num_accepted)<<endl;


		//release all resource to start the next round
		ReleaseAllResources();

		//break;

		//time_slot++;
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
	printf("Acceptance Ratio: %.8lf%%\n", 100*acceptance);
	FILE *accept_file = fopen((outPath+".acceptance").c_str(), "a");
	fprintf(accept_file, "%f\n", acceptance);
	fclose(accept_file);

	FILE *avrg_latency_file = fopen((outPath+".average_latency").c_str(), "a");
	fprintf(avrg_latency_file, "%f\n", totaldelay*1.0/double(stats.num_accepted));
	cout<<"average latency: "<<totaldelay*1.0/double(stats.num_accepted)<<endl;
	fclose(avrg_latency_file);

	//cout << "Write all the SC provision results in a file." << endl;
	FILE *sequence_file = fopen((outPath+".sequences").c_str(), "a");
	FILE *latency_file = fopen((outPath+".latency").c_str(), "a");
	for (int mcid = 0; mcid < all_scResults.size(); mcid++) {
		fprintf(latency_file, "%f ", traffic_requests[mcid].SFCdelay);
		//c out<<mcid<<": "<<traffic_requests[mcid].SFCdelay<<endl;
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
	fprintf(latency_file, "\n");
	fclose(sequence_file);
	fclose(latency_file);

	//cout << "Write solution time to file." << endl;
	//c out<<"solutiontimes.size(): "<<solutiontimes.size()<<endl;
	double lastVal = 0.0;
	FILE *soltionTime_file = fopen((outPath+".solutiontime").c_str(), "a");
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

	cout<<"notexpansum: "<<notexpansum<<" "<<bwmap.size()<<endl;
	if(parallel=="false"){
		FILE *each_bwconsumption = fopen((topology_file+".parasfc.each_bwconsumption").c_str(), "w");
		//c out<<topology_file+".each_bwconsumption"<<endl;
		for(int tmpi=0; tmpi<bwmap.size(); tmpi++){
			fprintf(each_bwconsumption, "%lu\n", bwmap[tmpi]);
			//c out<<bwmap[tmpi]<<endl;
		}
		fclose(each_bwconsumption);
	}
}

//parse the parameters, and load data from files
void loadData(int argc, char *argv[]) {
	auto arg_maps = ParseArgs(argc, argv);
	//string topology_file, traffic_file, mboxSpec_file;  //added by rob
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
			maxIterNum = atoi(argument.second.c_str()); //max iteration number for MA procedure
		} else if (argument.first == "--outPath") {
			outPath = argument.second; //path and filename prefix of output file
		} else if (argument.first == "--shareVNF") {
			shareVNF = argument.second; //"true" or "false"
		} else if (argument.first == "--TopoScale") {
			topoScale = atof(argument.second.c_str());
		} else if (argument.first == "--bwfactor") {
			bwfactor = atof(argument.second.c_str());
		} else if (argument.first == "--testbw") {
			testbw = atol(argument.second.c_str());
		} else if (argument.first == "--bwscale") {
			bwscale = atof(argument.second.c_str());
		}
	}
	if(parallel == "true"){
		//c out<<"initialize parallelized SFCs"<<endl;
		InitializeSFCMap();
	}
	else{
		//c out<<"initialize sequential SFCs"<<endl;
		InitializeSeqSFCMap();
	}

	//c out<<"load topology data: "<<topology_file<<endl;
	InitializeTopology(topology_file.c_str(), bwscale);

	//c out<<"load the middlebox information: "<<mboxSpec_file<<endl;
	InitializeMiddleboxes(mboxSpec_file.c_str());
	// PrintMiddleboxes();

	//c out<<"load the traffic request info.: "<<traffic_file<<endl;
	InitializeTrafficRequests_rob(traffic_file.c_str(),testbw);

	//if (parallel == "true")
	//	InitializeTrafficRequestsSubChains_rob(traffic_file.c_str());
}

map<int, unsigned long> getsSFCbw(string filename){
	filename = filename + ".parasfc.each_bwconsumption";
	cout<<filename<<endl;
	map<int, unsigned long> sSFCbw;
	FILE *sSFCres = fopen(filename.c_str(), "r");
	unsigned long tmpbw;
	while(fscanf(sSFCres, "%lu", &tmpbw)>0){
		sSFCbw[sSFCbw.size()] = tmpbw;
		//c out<<tmpbw<<" ";
	}
	//c out<<endl;
	fclose(sSFCres);

	//c out<<"sSFC bw size: "<<sSFCbw.size()<<endl;

	return sSFCbw;
}

int main(int argc, char *argv[]) {

	if (argc < 6) {
		puts(kUsage.c_str());
		return 1;
	}


	loadData(argc, argv); //parse the parameters, and load data from files

	shareVNF = "false";

	preProvision();
	adjustPriority();

	//c out<<"to getsSFCbw"<<endl;
	map<int, unsigned long>seqbwmap;
	if(parallel=="true"){
		seqbwmap = getsSFCbw(topology_file);
	}

	//c out<<"to toProvision"<<endl;
	toProvision(seqbwmap, bwfactor);

	return 0;
}
