#ifndef MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_
#define MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_

#include "datastructure_rob.h"
#include "util_rob.h"
#include <algorithm>

// t_request: a request's mainchain
std::unique_ptr<std::vector<int> > greedProvision(
		const traffic_request &t_request, const int cc, double &delay, bool* hasprov, int* provnode) {

	bool done = false;
	int loopNum = 0;
	bool found = false;
	loopNum++;

	delay = INF;
	bool isgood = true;
	int stage = 0;
	int cost = 0;
	int pre_node, cur_node;
	vector<int> path;
	const static int kNumNodes = graph.size();
	vector<int> ssfc = t_request.middlebox_sequence;
	int sid = t_request.sid;
	vector<int> subchain = pSFC[sid][cc];
	const int kNumStages = subchain.size(); //the length of sub chain
	//cout << "kNumStages: " << kNumStages << endl;


	resource resource_vector;
	for (int i = 0; i < kNumNodes; ++i) {
		resource_vector.cpu_cores.push_back(nodes[i].residual_cores);
		//cout<<nodes[i].residual_cores<<" ";
	}
	//cout<<endl;

	pre_node = t_request.source;
	path.push_back(pre_node);
	int endstage;
	if(subchain[kNumStages-1] == ssfc.size()-1)
		endstage = kNumStages-1;
	else
		endstage = kNumStages;

	for (stage = 1; stage < endstage; stage++) {

		//cout << "stage " << stage << ": "
		//		<< middleboxes[t_request.middlebox_sequence[stage]].middlebox_name
		//		<< endl;
		// get the mbox in current stage
		int fid = subchain[stage];
		int vnf = ssfc[fid];
		const middlebox &m_box = middleboxes[vnf];

		int min_node = NIL, min_cost = INF;

		if(hasprov[fid]){
			cur_node = provnode[fid];
			min_node = provnode[fid];

			bool BWprovisioned = false;
			unsigned long residualBW = GetPathResidualBandwidth(pre_node, cur_node);
			pair<int, int> tmpkey = make_pair(pre_node, cur_node);
			if(reservedEdges.find(tmpkey) != reservedEdges.end()){
				BWprovisioned = true;
				//residualBW -= reservedEdges[tmpkey];
			}

			//if ((GetPathResidualBandwidth(prev_node, current_node) >= t_request.min_bandwidth)) {
			if (residualBW >= t_request.min_bandwidth || BWprovisioned) {
				min_cost = cost + GetCost_rob(pre_node, cur_node, t_request)
						+ GetCost_rob(cur_node, t_request.destination, t_request);
			}
			else
			{
				isgood = false;
				break;
			}
		}
		else {
			for(cur_node = 0; cur_node < kNumNodes;cur_node++){
				// ......check whether the resource is enough for vNF deployment
				//cout<<"Check if current node has enough capacity for deploying the mbox..."<<endl;
				//cout << "pre node and cur node: " << pre_node << " " << cur_node << endl;
				if (IsResourceAvailable_Rob(pre_node, cur_node, resource_vector, m_box, t_request)) {

					//cout << "compute the cost from the prenode to current node" << endl;
					double tmpcost = cost + GetCost_rob(pre_node, cur_node, t_request);
							+ GetCost_rob(cur_node, t_request.destination, t_request);

					if (tmpcost < min_cost) {
						min_cost = tmpcost;
						min_node = cur_node;
					}
				}
			}
		}


		if(min_node != NIL){
			bool new_middlebox_deployed = true;

			// cout << "check whether mbox has existed" << endl;
			// to check whether the instance has been deployed and whether the capacity is enough
			// if do for both, reused the existing instances
			for (middlebox_instance &mbox_instance : deployed_mboxes[min_node]) {

				// vNF instance ***can not*** be shared by multi-user
				if (shareVNF == "false") {
					break;
				}

				// vNF instance can be shared by multi-user
				if (mbox_instance.m_box->middlebox_name == m_box.middlebox_name
						&& mbox_instance.residual_capacity >= t_request.min_bandwidth) {
					new_middlebox_deployed = false;
					break;
				}
			}

			//cout << "id a new instance should be deployed..." << endl;
			if (new_middlebox_deployed) {
				// node resource update
				resource_vector.cpu_cores[min_node] -= m_box.cpu_requirement;
			}


			path.push_back(min_node);
			cost = cost + GetCost_rob(pre_node, min_node, t_request);
			pre_node = min_node;
		}
		else{
			isgood = false;
			break;
		}

		//cout<<cost<<" ";
	}

	if(isgood && subchain[kNumStages-1] == ssfc.size()-1){
		path.push_back(t_request.destination);
		cost = cost + GetCost_rob(pre_node, t_request.destination, t_request);
		//cout<<"("<<pre_node<<" "<<t_request.destination<<" "<<GetCost_rob(pre_node, t_request.destination, t_request)<<") ";
	}
	//cout<<cost<<endl;
	delay = cost;


	if (!isgood) {
		rejectflag = 1;
		delay = INF;
		return std::unique_ptr < std::vector<int> > (new std::vector<int>());
	}
	else{
		// accept the request, and return the path
		std::unique_ptr<std::vector<int> > return_vector(new std::vector<int>());
		for (int pid = 0; pid < path.size(); pid++) {
			return_vector->push_back(path[pid]);
			int fid = subchain[pid];
			if(fid > 0 && fid < ssfc.size()-1){
				int vnf = ssfc[fid];
				const middlebox &m_box = middleboxes[vnf];
				delay += m_box.processing_delay;
			}
		}
		return std::move(return_vector);
	}

}

#endif  //  MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_
