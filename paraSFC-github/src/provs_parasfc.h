#ifndef MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_
#define MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_

#include "datastructure_rob.h"
#include "util_rob.h"
#include <algorithm>

void ViterbiInit() {
	for (int i = 0; i < 50; ++i) {
		for (int j = 0; j < 50; ++j) {
			cost[i][j] = INF;
			pre[i][j] = INF;
		}
	}

}

// t_request: a request
std::unique_ptr<std::vector<int> > ViterbiEstimate(
		const traffic_request &t_request, const int cc, double &delay) {

	//c out<<"viterbi computing..."<<endl;

	ViterbiInit(); // initialization

	int min_index;
	int stage;//, node = NIL;
	const static int kNumNodes = graph.size();
	vector<int> ssfc = t_request.middlebox_sequence;
	int sid = t_request.sid;
	vector<int> subchain = pSFC[sid][cc];
	const int kNumStages = subchain.size(); //the length of sub chain
	std::vector<resource> current_vector, previous_vector;
	current_vector.resize(kNumNodes);
	previous_vector.resize(kNumNodes);

	for (int i = 0; i < kNumNodes; ++i) {
		for (int j = 0; j < kNumNodes; ++j) {
			current_vector[i].cpu_cores.push_back(nodes[j].residual_cores);
			//if(i==0)cout<<nodes[j].residual_cores<< " ";
		}
	}
	//cout << endl;


	// for each vNF in the SFC
	stage = 0;
	cost[stage][t_request.source] = 0;
	int validstage, lastnode;

	if(subchain[kNumStages-1] == ssfc.size()-1){
		validstage = kNumStages-1;
		lastnode = t_request.destination;
	}
	else{
		validstage = kNumStages;
		lastnode = subchain[validstage-1];
	}


	for (stage = 1; stage < validstage; stage++) {

		//c out<<"state "<<stage<<": "<<middleboxes[t_request.middlebox_sequence[subchain[stage]]].middlebox_name<<endl;
		// get the mbox in current stage
		int fid = subchain[stage];
		int vnf = ssfc[fid];
		const middlebox &m_box = middleboxes[vnf];
		previous_vector = current_vector;

		// for each node which can deploy current vNF
		for (int current_node = 0; current_node < kNumNodes; current_node++) {
			if(nodes[current_node].num_cores==0)
				continue;

			// find out the best prev_node for current_node
			min_index = NIL;
			for (int prev_node = 0; prev_node < kNumNodes; prev_node++) {
				if(nodes[prev_node].num_cores==0)
					continue;

				if(stage == 1)
					prev_node = t_request.source;

				// ......the IsResourceAvailable() should consider the bandwidth resource
				// ......check whether the resource is enough
				//c out<<"m box capa: "<<m_box.processing_capacity<<" bw needed: "<<t_request.min_bandwidth<<endl;
				if (IsResourceAvailable(prev_node, current_node,
						previous_vector[prev_node], m_box, t_request)) {

					// ......compute the so-far cost from nodes in previous stage to current node
					double E2Elatency = cost[stage - 1][prev_node] + GetCost_rob(prev_node, current_node, t_request);

					//c out<<"stage: "<<stage<<" vnf: "<<m_box.middlebox_name.c_str()<<" prev node: "<<prev_node;
					//c out<<" cur node: "<<current_node;
					//c out<<" costs: "<<cost[stage - 1][prev_node] <<" "<< GetCost_rob(prev_node, current_node, t_request);
					//c out<<" "<<E2Elatency<<endl;

					// update the viterbi variables
					if (cost[stage][current_node] > E2Elatency) {
						cost[stage][current_node] = E2Elatency;
						pre[stage][current_node] = prev_node;
						min_index = prev_node;
					}
				}
				else
				{
					//c out<<"resoruce is insufficient!"<<endl;
				}

				if(stage == 1)
					break;
			}


			//c out<<cc<<" "<<current_node<<" : "<<min_index<<endl;

			// update the resource in the deployment node
			if (min_index != NIL) {

				// if current_node is chosen to deploy, then the ***resource status*** would become like this:
				current_vector[current_node].cpu_cores = previous_vector[min_index].cpu_cores;
				bool new_middlebox_deployed = true;

				// to check whether the instance has been deployed and whether the capacity is enough
				// if do for both, reused the existing instances
				for (middlebox_instance &mbox_instance : deployed_mboxes[current_node]) {

					// vNF instance ***can not*** be shared by multi-user
					if (shareVNF == "false") {
						break;
					}

					//c out<<"vNF instance can be shared by multi-user"<<endl;
					if (mbox_instance.m_box->middlebox_name
							== m_box.middlebox_name
							&& mbox_instance.residual_capacity
									>= t_request.min_bandwidth) {
						new_middlebox_deployed = false;
						break;
					}
				}

				//c out<<"id a new instance should be deployed..."<<endl;
				if (new_middlebox_deployed) {
					// node resource update
					current_vector[current_node].cpu_cores[current_node] -= m_box.cpu_requirement;

				}

			} else {
				// ...... is it means that the resource is run out??????????????
				//current_vector[current_node].cpu_cores.clear();
			}

		}//end of each node
	}//end of each stage


	//wrapping processing
	if(subchain[kNumStages-1] == ssfc.size()-1){

		// compute the last stage which transfers from the last vNF to egress node
		double min_cost = INF;
		int min_index = NIL;
		for (int prenode = 0; prenode < kNumNodes; prenode++) {
			if(nodes[prenode].num_cores==0)
				continue;

			if(stage == 1)
				prenode = t_request.source;

			// ......why is different from the former ones? because it does not have to consider
			// compute the cost from last node of requested SFC to the egress node
			double transition_cost = cost[kNumStages - 2][prenode]
					+ GetCost_rob(prenode, t_request.destination, t_request);

			//GetTransitCost(cur_node, t_request.destination, t_request) +
			//GetSLAViolationCost(cur_node, t_request.destination, t_request, fake_mbox);

			//c out<<prenode<<" "<<cost[kNumStages - 1][prenode]<<" ";
			//c out<<GetCost_rob(prenode, t_request.destination, t_request)<<" ";
			//c out<<min_cost<<" "<<transition_cost<<endl;

			// find out the node yielding the minimal cost
			if (min_cost > transition_cost) {
				min_cost = transition_cost;
				min_index = prenode;
				pre[stage][t_request.destination] = prenode;
				cost[stage][t_request.destination] = transition_cost;
			}

			if(stage == 1)
				break;
		}

		delay = min_cost;

	}

	//c out<<"backtracking..."<<endl;
	stage = kNumStages-1;
	if(subchain[stage] == ssfc.size()-1){
		lastnode = t_request.destination;
	}
	else{
		int min_index = INF;
		int min_cost  = INF;
		for (int node = 0; node < kNumNodes; node++) {
			if(cost[stage][node] < min_cost)
			{
				min_cost = cost[stage][node];
				lastnode = node;
			}
		}
		delay = min_cost;
	}

	//c out<<"if such node does not exist, reject the request."<<endl;
	if (pre[stage][lastnode] == INF) {
		rejectflag = 1;
		//c out<<"reject cc: "<<cc<<endl;
		//++stats.num_rejected;
		return std::unique_ptr < std::vector<int> > (new std::vector<int>());
	}
	// otherwise, accept the request
	//++stats.num_accepted;

	//c out<<"find the solution sequence"<<endl;
	std::unique_ptr<std::vector<int> > return_vector(new std::vector<int>());
	int current_node = lastnode;
	return_vector->push_back(current_node);
	//c out<<current_node<<" ";
	// backtrack the best path
	for (stage = kNumStages - 1; stage > 0; stage--) {
		current_node = pre[stage][current_node];
		return_vector->push_back(current_node);
		//c out<<current_node<<" ";
	}
	//c out<<endl;

	//c out<<"return results"<<endl;
	std::reverse(return_vector->begin(), return_vector->end());
	return std::move(return_vector);
}

/*********************************************************************************************
 *
 *parameters:
 *  t_request: a request
 *********************************************************************************************/
std::unique_ptr<std::vector<int> > ViterbiProvision(
		const traffic_request &t_request, const int cc, double &delay, bool* hasprov, int* provnode, int ssfclen) {

	//c out<<"viterbi computing..."<<endl;

	//ViterbiInit(); // initialization
	for (int i = 0; i < 50; ++i) {
		for (int j = 0; j < 50; ++j) {
			cost[i][j] = INF;
			pre[i][j] = INF;
		}
	}

	int min_index;
	int stage;//, node = NIL;
	const static int kNumNodes = graph.size();
	vector<int> ssfc = t_request.middlebox_sequence;
	int sid = t_request.sid;
	vector<int> subchain = pSFC[sid][cc];  //<0,1,2,3,4,5,11>
	if(ssfclen > 0){
		subchain.clear();
		for(int tmpi=0; tmpi< ssfclen; tmpi++){
			subchain.push_back(tmpi);
		}
	}

	const int kNumStages = subchain.size(); //the length of sub chain
	std::vector<resource> current_vector, previous_vector;
	current_vector.resize(kNumNodes);
	previous_vector.resize(kNumNodes);

	for (int i = 0; i < kNumNodes; ++i) {
		for (int j = 0; j < kNumNodes; ++j) {
			current_vector[i].cpu_cores.push_back(nodes[j].residual_cores);
			//if(i==0)cout<<nodes[j].residual_cores<< " ";
		}
	}
	//cout << endl;


	// for each vNF in the SFC
	stage = 0;
	cost[stage][t_request.source] = 0;
	int validstage, lastnode;

	if(subchain[kNumStages-1] == ssfc.size()-1){
		validstage = kNumStages-1;
		lastnode = t_request.destination;
	}
	else{
		validstage = kNumStages;
		lastnode = subchain[validstage-1];
	}


	for (stage = 1; stage < validstage; stage++) {

		//c out<<"state "<<stage<<": "<<middleboxes[t_request.middlebox_sequence[subchain[stage]]].middlebox_name<<endl;
		// get the mbox in current stage
		int fid = subchain[stage];
		int vnf = ssfc[fid];
		const middlebox &m_box = middleboxes[vnf];
		previous_vector = current_vector;

		if(hasprov[fid]){
			// find out the best prev_node for current_node
			min_index = NIL;
			int current_node = provnode[fid];
			for (int prev_node = 0; prev_node < kNumNodes; prev_node++) {
				if(nodes[prev_node].num_cores==0)
					continue;
				if(stage==1)
					prev_node = t_request.source;

				// ......the IsResourceAvailable() should consider the bandwidth resource
				// ......check whether the resource is enough
				//c out<<"m box capa: "<<m_box.processing_capacity<<" bw needed: "<<t_request.min_bandwidth<<endl;

				bool BWprovisioned = false;
				long residualBW;
				if(prev_node == current_node){
					residualBW=1999999999;
				}
				else{
					residualBW=GetPathResidualBandwidth(prev_node, current_node);
				}

				pair<int, int> tmpkey = make_pair(prev_node, current_node);
				if(reservedEdges.find(tmpkey) != reservedEdges.end()){
					BWprovisioned = true;
					residualBW -= reservedEdges[tmpkey];
				}


				//if ((GetPathResidualBandwidth(prev_node, current_node) >= t_request.min_bandwidth)) {
				if (residualBW >= t_request.min_bandwidth || BWprovisioned) {

					// ......compute the so-far cost from nodes in previous stage to current node
					double E2Elatency = cost[stage - 1][prev_node] + GetCost_rob(prev_node, current_node, t_request);

					//c out<<"stage: "<<stage<<" vnf: "<<m_box.middlebox_name.c_str()<<" prev node: "<<prev_node;
					//c out<<" cur node: "<<current_node;
					//c out<<" costs: "<<cost[stage - 1][prev_node] <<" "<< GetCost_rob(prev_node, current_node, t_request);
					//c out<<" "<<E2Elatency<<endl;

					// update the viterbi variables
					if (cost[stage][current_node] > E2Elatency) {
						cost[stage][current_node] = E2Elatency;
						pre[stage][current_node] = prev_node;
						min_index = prev_node;
					}
				}
				else
				{
					//cout<<"bandwidth is insufficient!"<<residualBW<<" for "<<prev_node
					//		<<" "<<current_node<<" "<<GetPathResidualBandwidth(prev_node, current_node)<<endl;
				}

				if(stage==1)
					break;
			}
		}
		else{
			// for each node which can deploy current vNF
			for (int current_node = 0; current_node < kNumNodes; current_node++) {

				if(nodes[current_node].num_cores==0)
					continue;

				// find out the best prev_node for current_node
				min_index = NIL;
				for (int prev_node = 0; prev_node < kNumNodes; prev_node++) {

					if(nodes[prev_node].num_cores<=0)
						continue;

					if(stage==1)
						prev_node = t_request.source;

					// ......the IsResourceAvailable() should consider the bandwidth resource
					// ......check whether the resource is enough
					//c out<<"m box capa: "<<m_box.processing_capacity<<" bw needed: "<<t_request.min_bandwidth<<endl;
					if (IsResourceAvailable_Rob(prev_node, current_node, previous_vector[prev_node], m_box, t_request)) {

						// ......compute the so-far cost from nodes in previous stage to current node
						double E2Elatency = cost[stage - 1][prev_node] + GetCost_rob(prev_node, current_node, t_request);

						//c out<<"stage: "<<stage<<" vnf: "<<m_box.middlebox_name.c_str()<<" prev node: "<<prev_node;
						//c out<<" cur node: "<<current_node;
						//c out<<" costs: "<<cost[stage - 1][prev_node] <<" "<< GetCost_rob(prev_node, current_node, t_request);
						//c out<<" "<<E2Elatency<<endl;

						// update the viterbi variables
						if (cost[stage][current_node] > E2Elatency) {
							cost[stage][current_node] = E2Elatency;
							pre[stage][current_node] = prev_node;
							min_index = prev_node;
						}
					}
					else
					{
						//cout<<"resoruce is insufficient!"<<endl;
					}


					if(stage==1)
						break;
				}


				//c out<<cc<<" "<<current_node<<" : "<<min_index<<endl;

				// update the resource in the deployment node
				if (min_index != NIL) {

					// if current_node is chosen to deploy, then the ***resource status*** would become like this:
					current_vector[current_node].cpu_cores = previous_vector[min_index].cpu_cores;
					bool new_middlebox_deployed = true;

					// to check whether the instance has been deployed and whether the capacity is enough
					// if do for both, reused the existing instances
					for (middlebox_instance &mbox_instance : deployed_mboxes[current_node]) {

						// vNF instance ***can not*** be shared by multi-user
						if (shareVNF == "false") {
							break;
						}

						//c out<<"vNF instance can be shared by multi-user"<<endl;
						if (mbox_instance.m_box->middlebox_name == m_box.middlebox_name
								&& mbox_instance.residual_capacity >= t_request.min_bandwidth) {
							new_middlebox_deployed = false;
							break;
						}
					}

					//c out<<" a new instance should be deployed..."<<endl;
					if (new_middlebox_deployed) {
						// node resource update
						current_vector[current_node].cpu_cores[current_node] -= m_box.cpu_requirement;

					}

				} else {

					// ...... is it means that the resource is run out??????????????
					current_vector[current_node].cpu_cores.clear();
				}

			}//end of each node

		}//

	}//end of each stage


	//wrapping processing
	if(subchain[kNumStages-1] == ssfc.size()-1){

		// compute the last stage which transfers from the last vNF to egress node
		double min_cost = INF;
		int min_index = NIL;
		for (int prenode = 0; prenode < kNumNodes; prenode++) {

			if(nodes[prenode].num_cores<=0)
				continue;

			// ......why is different from the former ones? because it does not have to consider
			// compute the cost from last node of requested SFC to the egress node
			double transition_cost = cost[kNumStages - 2][prenode]
					+ GetCost_rob(prenode, t_request.destination, t_request);

			//GetTransitCost(cur_node, t_request.destination, t_request) +
			//GetSLAViolationCost(cur_node, t_request.destination, t_request, fake_mbox);

			//c out<<prenode<<" "<<cost[kNumStages - 1][prenode]<<" ";
			//c out<<GetCost_rob(prenode, t_request.destination, t_request)<<" ";
			//c out<<min_cost<<" "<<transition_cost<<endl;

			// find out the node yielding the minimal cost
			if (min_cost > transition_cost) {
				min_cost = transition_cost;
				min_index = prenode;
				pre[stage][t_request.destination] = prenode;
				cost[stage][t_request.destination] = transition_cost;
			}
		}

		delay = min_cost;
		//c out<<"("<<delay<<") ";

	}

	//c out<<"backtracking..."<<endl;
	stage = kNumStages-1;
	if(subchain[stage] == ssfc.size()-1){
		lastnode = t_request.destination;
	}
	else{
		int min_index = NIL;
		int min_cost  = INF;
		for (int node = 0; node < kNumNodes; node++) {
			if(cost[stage][node] < min_cost)
			{
				min_cost = cost[stage][node];
				lastnode = node;
			}
		}
		delay = min_cost;
		//c out<<"["<<delay<<"] ";
	}

	//c out<<"if such node does not exist, reject the request."<<endl;
	if (pre[stage][lastnode] == INF) {
		rejectflag = 1;
		//c out<<"reject cc: "<<cc<<endl;
		//++stats.num_rejected;
		return std::unique_ptr < std::vector<int> > (new std::vector<int>());
	}
	// otherwise, accept the request
	//++stats.num_accepted;

	//c out<<"find the solution sequence"<<endl;
	std::unique_ptr<std::vector<int> > return_vector(new std::vector<int>());
	int current_node = lastnode;
	return_vector->push_back(current_node);
	//c out<<current_node<<" ";
	// backtrack the best path
	for (stage = kNumStages - 1; stage > 0; stage--) {
		current_node = pre[stage][current_node];
		return_vector->push_back(current_node);
		int fid = subchain[stage];
		if(fid < ssfc.size()-1){
			int vnf = ssfc[fid];
			const middlebox &m_box = middleboxes[vnf];
			delay += m_box.processing_delay;
		}
	}

	//c out<<"return results"<<endl;
	std::reverse(return_vector->begin(), return_vector->end());
	return std::move(return_vector);
}

#endif  //  MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_
