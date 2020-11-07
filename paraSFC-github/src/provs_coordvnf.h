#ifndef MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_
#define MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_

#include "datastructure_rob.h"
#include "util_rob.h"
#include <algorithm>

// t_request: a request's mainchain
std::unique_ptr<std::vector<int> > coordProvision(
		const traffic_request &t_request, const int cc, double &delay, bool* hasprov, int* provnode) {

	//cout<<"Provisioning "<<cc<<"-th CC in sid-"<<t_request.sid<<" of "<<t_request.originalID<<"-th request using coordVNF."<<endl;
	int stage = 0;
	int cost = 0;
	int pre_node, cur_node;
	bool found = false;
	bool done = false;
	int markstage = -1;
	int goforth = 0, backtrack = 0;
	vector<int> path;
	map<int, map<int, bool> > checklist;

	vector<int> ssfc = t_request.middlebox_sequence;
	int sid = t_request.sid;
	vector<int> subchain = pSFC[sid][cc];
	const static int kNumNodes = graph.size();
	const int kNumStages = subchain.size(); //the length of sub chain

	pre_node = t_request.source;
	path.push_back(pre_node);
	cur_node = pre_node;

	//if(subchain[kNumStages-1] == ssfc.size()-1)

	resource resource_vector;
	for (int i = 0; i < kNumNodes; i++) {
		resource_vector.cpu_cores.push_back(nodes[i].residual_cores);
		//cout<<nodes[i].residual_cores<<" ";
	}
	//cout<<endl;
	map<int, resource> resource_map;
	resource_map[-1] = resource_vector;

	for (stage = 1; stage < kNumStages; stage++) {
		if((subchain[kNumStages-1] == ssfc.size()-1) && stage == kNumStages-1){
			//cout << "deploy the last segment of path" << endl;
			//compute the last stage which transfers from the last vNF to egress node

			double tmpcost = cost + GetCost_rob(cur_node, t_request.destination, t_request);
			cost = tmpcost;
			path.push_back(t_request.destination);
			done = true;

			// end of stage == kNumStages-1
		}
		else{

			//get the mbox in current stage
			int fid = subchain[stage];
			int vnf = ssfc[fid];
			const middlebox &m_box = middleboxes[vnf];

			if(hasprov[fid]){
				cur_node = provnode[fid];

				if((checklist.find(stage)) != checklist.end()){
					if((checklist[stage].find(cur_node)) != checklist[stage].end()){
						if(checklist[stage][cur_node])
							break;
					}
				}

				map<int, bool> tmpmap;
				tmpmap[cur_node] = true;
				checklist[stage] = tmpmap;
				found = true;
				path.push_back(cur_node);
				pre_node = cur_node;
				cost = cost + GetCost_rob(pre_node, cur_node, t_request);
				resource_map[stage] = resource_vector;

			}else{

				//cout<<"find out all available nodes within 2 edge to pre_node in graph"<<endl;
				if((checklist.find(stage)) == checklist.end()){
					map<int, bool> tmpmap;
					for(int tmpi  = 0; tmpi < graph[pre_node].size(); tmpi++)
					{
						int tmpid = graph[pre_node][tmpi].u->node_id;
						tmpmap[tmpid]  = false;
						for(int tmpj   = 0; tmpj < graph[tmpid].size(); tmpj++){
							int tmpid2 = graph[tmpid][tmpj].u->node_id;
							tmpmap[tmpid2] = false;
						}
					}
					tmpmap[pre_node] = false;
					checklist[stage] = tmpmap;
				}


		//		map<int, bool>::iterator tmpiit;
		//		for(tmpiit = checklist[stage].begin(); tmpiit!= checklist[stage].end(); tmpiit++)
		//			if(!tmpiit->second)
		//				cout<<tmpiit->first<<" ";
		//		cout<<endl;

				found = false;
				for(int nodeid = -1; nodeid < kNumNodes; nodeid++)
				{
					if(nodeid == -1){  // Firstly, try to deploy in pre node
						cur_node = pre_node;
						//cout<<"Firstly, try to deploy in pre node"<<endl;
					}
					else if(nodeid == pre_node) // has been considered
						continue;
					else if((checklist[stage].find(nodeid)) == checklist[stage].end())  // the node is not with 2-order
						continue;
					else
						cur_node = nodeid;

					 // The node has been considered
					if(checklist[stage][cur_node])
						continue;

					checklist[stage][cur_node] = true;

					// ......check whether the resource is enough for vNF deployment
					//cout<<"Check if current node has enough capacity for deploying the mbox..."<<endl;
					//cout << "pre node & cur node: " << pre_node << " " << cur_node << endl;
					if (IsResourceAvailable_Rob(pre_node, cur_node, resource_vector, m_box, t_request)) {

						//cout << "compute the cost from the prenode to current node" << endl;
						double tmpcost = cost + GetCost_rob(pre_node, cur_node, t_request);

						found = true;
						pre_node = cur_node;
						path.push_back(cur_node);
						cost = tmpcost;
						bool new_middlebox_deployed = true;

						// cout << "check whether mbox has existed" << endl;
						// to check whether the instance has been deployed and whether the capacity is enough
						// if do for both, reused the existing instances
						for (middlebox_instance &mbox_instance : deployed_mboxes[cur_node]) {

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

						//cout << "a new instance should be deployed..." << endl;
						if (new_middlebox_deployed) {
							// node resource update
							resource_vector.cpu_cores[cur_node] -= m_box.cpu_requirement;
							resource_map[stage] = resource_vector;
						}

						break;   //cur stage is deployed
					}// end of checking whether the resource is enough for vNF deployment
				}// end for each node
			}


			if(found && stage == kNumStages-1)//((subchain[kNumStages-1] != ssfc.size()-1) && stage == kNumStages-1)
				done = true;

		}
		//cout<<"found, done: "<<found<<" "<<done<<endl;
		if (!found) {// backtrack
			if(path.size()>1){
				path.pop_back();
				resource_vector = resource_map[stage];
			}
			pre_node = path.back();
			checklist.erase(stage);
			stage -= 2;
		}
		if(stage<=0)
			break;


		//cout<<"end of stage "<<stage<<endl;
	}// end for each stage




	//cout<<"finished deployment :) :) :) :)"<<endl;
	if (!done) {
		rejectflag = 1;
		return std::unique_ptr < std::vector<int> > (new std::vector<int>());
	}

	delay = cost;

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

#endif  //  MIDDLEBOX_PLACEMENT_SRC_VITERBI_H_
