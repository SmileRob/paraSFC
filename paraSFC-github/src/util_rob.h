#ifndef MIDDLEBOX_PLACEMENT_SRC_UTIL_H_
#define MIDDLEBOX_PLACEMENT_SRC_UTIL_H_
#include "datastructure_rob.h"

#include <algorithm>
#include <assert.h>
#include <set>
#include <stack>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#define ONE_GIG 1000000000ULL
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__) " "

#ifdef DBG
#define DEBUG(...) PrintDebugMessage(AT, __VA_ARGS__)
#else
#define DEBUG(...)
#endif

void PrintDebugMessage(const char *location, const char *fmt_string, ...) {
	va_list args;
	va_start(args, fmt_string);
	std::string str = location;
	str += fmt_string;
	vprintf(str.c_str(), args);
	fflush(stdout);
	va_end(args);
}

inline unsigned long CurrentTimeNanos() {
	timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return static_cast<unsigned long>(ts.tv_sec)
			+ static_cast<unsigned long>(ts.tv_nsec);
}

template<class T>
double GetMean(const std::vector<T> &data) {
	T sum = T(0);
	const size_t kNumElements = data.size();
	for (auto &element : data)
		sum += element;
	return sum / static_cast<T>(kNumElements);
}

template<class T>
T GetNthPercentile(const std::vector<T> &data, int n) {

	std::vector<T> temp_data_buffer = data;
	sort(temp_data_buffer.begin(), temp_data_buffer.end());
	const size_t kNumElements = data.size();
	int rank = n * kNumElements;

	// ??? what for?
	if (rank % 100) {
		rank = (rank / 100) + 1;
	} else
		rank /= 100;
	--rank;

	return temp_data_buffer[rank];
}

template<class T>
std::vector<std::pair<T, double> > GetCDF(const std::vector<T> &data) {
	int precision = 1;
	std::vector<T> temp_data_buffer = data;
	if (typeid(temp_data_buffer[0]) == typeid(double)
			|| typeid(temp_data_buffer[0]) == typeid(float)) {
		precision = 1000;
	}
	std::map<int, int> cdf;
	for (int i = 0; i < temp_data_buffer.size(); ++i) {
		int bucket_index = temp_data_buffer[i] * precision;
		if (cdf[bucket_index])
			cdf[bucket_index]++;
		else
			cdf[bucket_index] = 1;
	}
	std::map<int, int>::iterator prev = cdf.begin(), current = cdf.begin();
	current++;
	for (; current != cdf.end(); current++, prev++) {
		current->second += prev->second;
	}
	int total = temp_data_buffer.size();
	std::vector<std::pair<T, double> > ret;
	for (current = cdf.begin(); current != cdf.end(); ++current) {
		T first = static_cast<T>(current->first) / static_cast<T>(precision);
		double second = static_cast<double>(current->second)
				/ static_cast<double>(total);
		ret.push_back(std::make_pair(first, second));
	}
	return ret;
}

inline std::unique_ptr<std::vector<int> > ComputeShortestPath(int source, int destination) {
	//cout<<"src/dst: "<<source<<" "<<destination<<": ";
	std::unique_ptr<std::vector<int> > path(new std::vector<int>());
	while (destination != NIL) {
		//cout<<destination<<"<-";
		path->push_back(destination);
		destination = sp_pre[source][destination];
	}
	//cout<<endl;
	std::reverse(path->begin(), path->end());
	return std::move(path);
}

inline void RefreshServerStats(int timestamp) {
	for (auto &node : nodes) {
		if (node.num_cores > 0) {
			double utilization = static_cast<double>(node.num_cores - node.residual_cores)
					/ static_cast<double>(node.num_cores);
			stats.server_stats.emplace_back(timestamp, node.node_id, utilization);
		}
	}
}


inline unsigned long GetEdgeResidualBandwidth(int source, int destination) {
	return bw[source][destination];
}

inline unsigned long GetPathResidualBandwidth(int source, int destination) {
	std::vector<int> *path_ptr = nullptr;
	std::pair<int, int> cache_index(source, destination);
	if (path_cache[cache_index]) {
		path_ptr = path_cache[cache_index].get();
	} else {
		path_cache[cache_index] = ComputeShortestPath(source, destination);
		path_ptr = path_cache[cache_index].get();
	}

	unsigned long residual_bandwidth = 100000000000000L;
	//unsigned long residual_bandwidth = 0;
	for (int i = 0; i < static_cast<int>(path_ptr->size()) - 1; ++i) {
		DEBUG("edge[%d][%d] = %d\n", path_ptr->at(i), path_ptr->at(i + 1),
				GetEdgeResidualBandwidth(path_ptr->at(i), path_ptr->at(i + 1)));

		residual_bandwidth = std::min(residual_bandwidth,
				GetEdgeResidualBandwidth(path_ptr->at(i), path_ptr->at(i + 1)));
	}

	if(residual_bandwidth==100000000000000L){   //added by Rob 20200907
		residual_bandwidth=0;
	}

	return residual_bandwidth;
}

inline void ReduceEdgeResidualBandwidth(int source, int destination, unsigned long bandwidth) {
	bw[source][destination] -= bandwidth;
	bw[destination][source] -= bandwidth;
}

void DecommissionAllMiddleboxes() {
	for (auto &mboxes : deployed_mboxes)
		mboxes.clear();
}

void ReleaseBandwidth() {
	for (int i = 0; i < graph.size(); ++i) {
		auto &adj_list = graph[i];
		for (auto &endpoint : adj_list) {
			endpoint.residual_bandwidth = endpoint.bandwidth;
			bw[i][endpoint.u->node_id] = bw[endpoint.u->node_id][i] = endpoint.bandwidth;
		}
	}
}

void ReleaseCPU() {
	for (auto &n : nodes) {
		n.residual_cores = n.num_cores;
	}
}

inline void ReleaseAllResources() {
	DecommissionAllMiddleboxes();
	ReleaseCPU();
	ReleaseBandwidth();
}

inline void ReducePathResidualBandwidth(int source, int destination, unsigned long bandwidth) {
	std::pair<int, int> cache_index(source, destination);
	std::vector<int> *path_ptr = nullptr;
	if (path_cache[cache_index]) {
		path_ptr = path_cache[cache_index].get();
	} else {
		auto path = ComputeShortestPath(source, destination);
		path_cache[cache_index] = std::move(path);
		path_ptr = path_cache[cache_index].get();
	}

	for (int i = 0; i < static_cast<int>(path_ptr->size()) - 1; ++i) {
		ReduceEdgeResidualBandwidth(path_ptr->at(i), path_ptr->at(i + 1), bandwidth);
	}
}

inline void ReduceNodeCapacity(int node, const middlebox &m_box) {
	nodes[node].residual_cores -= m_box.cpu_requirement;
}

int UsedMiddleboxIndex(int current_node, const middlebox &m_box, const traffic_request &t_request) {
	if(shareVNF == "false")
		return NIL;
	for (int i = 0; i < deployed_mboxes[current_node].size(); ++i) {
		if (deployed_mboxes[current_node][i].m_box->middlebox_name
				== m_box.middlebox_name) {
			if (deployed_mboxes[current_node][i].residual_capacity
					>= t_request.min_bandwidth) {
				return i;
			}
		}
	}
	return NIL;
}

void UpdateMiddleboxInstances(int current_node, const middlebox *m_box,
		const traffic_request &t_request) {

	//added by rob
	int used_middlebox_index;
	if (shareVNF == "true")
		used_middlebox_index = UsedMiddleboxIndex(current_node, *m_box, t_request);
	else if (shareVNF == "false")
		used_middlebox_index = NIL;

	//int used_middlebox_index = UsedMiddleboxIndex(current_node, *m_box, t_request);
	if (used_middlebox_index != NIL) {
		deployed_mboxes[current_node][used_middlebox_index].residual_capacity -= t_request.min_bandwidth;
	} else {
		deployed_mboxes[current_node].emplace_back(m_box, m_box->processing_capacity - t_request.min_bandwidth);
		ReduceNodeCapacity(current_node, *m_box);
	}
}

void UpdateResources(const std::vector<int> *traffic_sequence, const traffic_request &t_request) {

	for (int i = 0; i < static_cast<int>(traffic_sequence->size()) - 1; ++i) {

		ReducePathResidualBandwidth(traffic_sequence->at(i), traffic_sequence->at(i + 1), t_request.min_bandwidth);
	}
	for (int i = 1; i < static_cast<int>(traffic_sequence->size()) - 1; ++i) {
		const middlebox &m_box = middleboxes[t_request.middlebox_sequence[i - 1]];
		UpdateMiddleboxInstances(traffic_sequence->at(i), &middleboxes[t_request.middlebox_sequence[i - 1]], t_request);
	}
}

void UpdateResources_Rob(const std::vector<int> *traffic_sequence, const traffic_request &t_request, std::vector<int> s_SFC) {

	for (int i = 0; i < static_cast<int>(traffic_sequence->size()) - 1; ++i) {
		ReducePathResidualBandwidth(traffic_sequence->at(i), traffic_sequence->at(i + 1), t_request.min_bandwidth);
	}
	for (int i = 1; i < static_cast<int>(traffic_sequence->size()) - 1; ++i) {

		int vnfID = s_SFC[t_request.middlebox_sequence[i]];
		const middlebox &m_box = middleboxes[vnfID];
		UpdateMiddleboxInstances(traffic_sequence->at(i), &middleboxes[vnfID], t_request);
	}
}


// update the resource for subchains
void UpdateResourcesForSubchain(const std::vector<int> *traffic_sequence,
		const traffic_request &t_request) {

	if (traffic_sequence->size() == 2) {
		//cout << "reduce the bandwidth overhead from links in the path";
		//cout << ", traffic_sequence size: " << traffic_sequence->size() << endl;
		ReducePathResidualBandwidth(traffic_sequence->at(0), traffic_sequence->at(1), t_request.min_bandwidth);

		//cout << "reduce the node corese from node's capacity" << endl;
		const middlebox &m_box = middleboxes[t_request.middlebox_sequence[1]];
		UpdateMiddleboxInstances(traffic_sequence->at(1), &middleboxes[t_request.middlebox_sequence[1]], t_request);
	}

}

inline int IsResourceAvailable(int prev_node, int current_node,
		const resource &resource_vector, const middlebox &m_box,
		const traffic_request &t_request) {

	//cout<<"firstly check whether the bandwidth is enough for the request"<<endl;
	if ((GetPathResidualBandwidth(prev_node, current_node) >= t_request.min_bandwidth)
			|| prev_node == current_node) {

		//cout<<"Check if we can use existing middlebox of the same type."<<endl;
		if (UsedMiddleboxIndex(current_node, m_box, t_request) != NIL) {
			return 1;
		}

		//cout<<"If we cannot use existing ones, then we need to instantiate new one"<<endl;
		if (m_box.processing_capacity >= t_request.min_bandwidth
				&& resource_vector.cpu_cores[current_node] >= m_box.cpu_requirement) {
			//cout<<"we can instantiate a new mbox."<<endl;
			return 1;
		}

		if (0) {
			if (m_box.processing_capacity < t_request.min_bandwidth)
				cout << "mbox's processing capacity is not enough!" << endl;
			if (resource_vector.cpu_cores[current_node] < m_box.cpu_requirement)
				cout << "node " << current_node << "'s cpu cores is not enough!" << endl;
		}

	}
	else{
		//cout<<"bandwidth: "<<GetPathResidualBandwidth(prev_node, current_node)
		//		<<" is not enough for "<<prev_node<<" to "<<current_node<<" for SFC "<<t_request.mcID
		//		<<" needed "<<t_request.min_bandwidth<<endl;
	}

	//cout<<"not available resource"<<endl;
	return 0;
}

inline int IsResourceAvailable_Rob(int prev_node, int current_node,
		const resource &resource_vector, const middlebox &m_box,
		const traffic_request &t_request) {

	//cout<<"firstly check whether the bandwidth is enough for the request"<<endl;

	bool BWprovisioned = false;
	long residualBW = GetPathResidualBandwidth(prev_node, current_node);
	if(prev_node == current_node){
		residualBW = 1999999999;
	}
	//cout<<"residual BW: "<<residualBW<<endl;
	pair<int, int> tmpkey = make_pair(prev_node, current_node);
	if(reservedEdges.find(tmpkey) != reservedEdges.end()){
		BWprovisioned = true;
		residualBW -= reservedEdges[tmpkey];
		//cout<<"GetPathResidualBandwidth(prev_node, current_node) ";
		//cout<<GetPathResidualBandwidth(prev_node, current_node)<<" "
		//		<<reservedEdges[tmpkey]<<" "<<residualBW<<endl;

	}

	//if ((GetPathResidualBandwidth(prev_node, current_node) >= t_request.min_bandwidth)) {
	if (residualBW >= t_request.min_bandwidth || BWprovisioned) {
	//if (residualBW >= t_request.min_bandwidth) {

		//cout<<"Check if we can use existing middlebox of the same type."<<endl;
		if (UsedMiddleboxIndex(current_node, m_box, t_request) != NIL) {
			//cout<<"UsedMiddleboxIndex "<<endl;
			return 1;
		}

		int residualCores = resource_vector.cpu_cores[current_node];
		if(reservedCores.find(current_node)!=reservedCores.end()){
			residualCores -= reservedCores[current_node];
		}
		//cout<<current_node<<": "<<resource_vector.cpu_cores[current_node]<<" "
		//		<<reservedCores[current_node]<<" "<<residualCores<<" "<<m_box.cpu_requirement<<endl;

		//cout<<"If we cannot use existing ones, then we need to instantiate new one"<<endl;
		//if (m_box.processing_capacity >= t_request.min_bandwidth
		//		&& resource_vector.cpu_cores[current_node] >= m_box.cpu_requirement) {
		if (m_box.processing_capacity >= t_request.min_bandwidth
				&& residualCores >= m_box.cpu_requirement) {
			//cout<<"we can instantiate a new mbox."<<endl;
			return 1;
		}

		if (0) {
			if (m_box.processing_capacity < t_request.min_bandwidth)
				cout << "mbox's processing capacity is not enough!" << endl;
			if (resource_vector.cpu_cores[current_node] < m_box.cpu_requirement)
				cout << "node " << current_node << "'s cpu cores is not enough!" << endl;
		}

	}

	//cout<<"not available resource"<<endl;
	return 0;
}

inline double GetSLAViolationCost(int prev_node, int current_node,
		const traffic_request &t_request, const middlebox &m_box) {
	const int kNumSegments = t_request.middlebox_sequence.size() + 1;
	const double kPerSegmentLatencyBound = (1.0 * t_request.max_delay)
			/ kNumSegments;

	if (shortest_path[prev_node][current_node] + m_box.processing_delay
			> kPerSegmentLatencyBound)
		return (shortest_path[prev_node][current_node] + m_box.processing_delay
				- kPerSegmentLatencyBound) * t_request.delay_penalty;
	return 0.0;
}

double GetSLAViolationCost(int source, int destination, double max_delay,
		double penalty) {
	if (shortest_path[source][destination] > max_delay) {
		return penalty * (shortest_path[source][destination] - max_delay);
	}
	return 0;
}

inline double GetTransitCost(int prev_node, int current_node,
		const traffic_request &t_request) {
	int path_length = shortest_edge_path[prev_node][current_node];
	if (path_length >= INF)
		return INF;
	return (1.0 / 1000.0) * path_length * per_bit_transit_cost
			* t_request.min_bandwidth * t_request.duration;
}

double GetServerEnergyConsumption(int num_cores_used) {
	int full_servers_used = num_cores_used / NUM_CORES_PER_SERVER;
	double energy_consumed = static_cast<double>(full_servers_used
			* SERVER_PEAK_ENERGY);
	int residual_cores = num_cores_used % NUM_CORES_PER_SERVER;
	energy_consumed += POWER_CONSUMPTION_ONE_SERVER(residual_cores);
	return energy_consumed;
}

inline double GetEnergyCost(int current_node, const middlebox &m_box,
		const resource &resource_vector, const traffic_request &t_request) {
	if (UsedMiddleboxIndex(current_node, m_box, t_request) != NIL) {
		return 0;
	}
	int previously_used_cores = nodes[current_node].num_cores
			- resource_vector.cpu_cores[current_node];
	int currently_used_cores = previously_used_cores + m_box.cpu_requirement;
	double duration_hours = static_cast<double>(t_request.duration)
			/ (60.0 * 60.0);
	double previous_cost = GetServerEnergyConsumption(previously_used_cores) *
	duration_hours * PER_UNIT_ENERGY_PRICE;
	double current_cost = GetServerEnergyConsumption(currently_used_cores) *
	duration_hours * PER_UNIT_ENERGY_PRICE;
	double energy_cost = current_cost - previous_cost;
	if (previously_used_cores == 0) {
		energy_cost += (GetServerEnergyConsumption(0) * duration_hours *
		PER_UNIT_ENERGY_PRICE);
	}
	return energy_cost;
}

inline double GetDeploymentCost(int current_node, const middlebox &m_box,
		const traffic_request &t_request) {
	// If we can use existing middlebox then there is no deployment cost.
	if (UsedMiddleboxIndex(current_node, m_box, t_request) != NIL) {
		return 0.0;
	}
	return m_box.deployment_cost;
}

// get the delay when flow traversing from prev_node to current node.
double GetCost_rob(int prev_node, int current_node, const traffic_request &t_request) {
	//delay = processing delay + transmission delay + transporting delay
	double transmitDelay = 0;
	double transportDelay = 0;

	std::vector<int> *path_ptr = nullptr;
	std::pair<int, int> cache_index(prev_node, current_node);
	if (path_cache_rob[cache_index]) {
		path_ptr = path_cache_rob[cache_index].get();
	} else {
		path_cache_rob[cache_index] = ComputeShortestPath(prev_node, current_node);
		path_ptr = path_cache_rob[cache_index].get();
	}

	//transportDelay is summation of link latency over all links comprising the path from prev to current node
	for (int i = 0; i < static_cast<int>(path_ptr->size()) - 1; ++i) {
		transportDelay += shortest_path[path_ptr->at(i)][path_ptr->at(i + 1)];
		//cout<<shortest_path[path_ptr->at(i)][path_ptr->at(i + 1)]<<" ";
	}
	//if shortest_path save the end to end latency for each path, we can have
	//transportDelay = shortest_path[prev_node][current_node];

	//cout<<" transportdelay: "<<transportDelay<<endl;
	//transportDelay = (vol / bandwidth)*  (number of links)
	//transmitDelay = (1.0 * t_request.dataVol / 10000)
	//		* shortest_edge_path[prev_node][current_node]; // in ms, assume bandwidth of link is 10Gbps

	// suppose the user data is divided into 1000 blocks,
	// so the transmit delay is just cause by only one blocks
	//transmitDelay = 1.0 * t_request.dataVol / t_request.min_bandwidth;// in ms, assume bandwidth of link is 10Gbps

	//cout<<" transmitDelay: "<<transmitDelay<<endl;
	DEBUG("transmitDelay = %lf, transportDelay = %lf\n", transmitDelay, transportDelay);
	return transmitDelay + transportDelay;
}

inline int GetLatency(int source, int destination) {
	for (edge_endpoint endpoint : graph[source]) {
		if (endpoint.u->node_id == destination)
			return endpoint.delay;
	}
	return NIL;
}

unsigned long GetBandwidthUsage(const std::vector<int> &traffic_sequence,
		const traffic_request &t_request) {
	unsigned long bandwidth_usage = 0;
	for (int i = 0; i < traffic_sequence.size() - 1; ++i) {
		bandwidth_usage += (t_request.min_bandwidth * (ComputeShortestPath(traffic_sequence[i], traffic_sequence[i + 1])->size() - 1));
	}
	return bandwidth_usage;
}

map<pair<int, int>, long> GetLinkUsage(const std::vector<int> &traffic_sequence, const traffic_request &t_request) {
	map<pair<int, int>, long> edgeutil;
	for (int i = 0; i < traffic_sequence.size() - 1; ++i) {
		std::unique_ptr<std::vector<int> > path = ComputeShortestPath(traffic_sequence[i], traffic_sequence[i + 1]);
		if(path->size()>1){
			int src = path->at(0);
			for(int j = 1; j<path->size(); j++)
			{
				int dst = path->at(j);
				edgeutil[make_pair(src,dst)] += t_request.min_bandwidth;
				src = dst;
			}
		}
	}
	return edgeutil;
}

unsigned long GetTotalNetworkBandwidth() {
	unsigned long total_bandwidth = 0;
	for (int i = 0; i < graph.size(); ++i) {
		for (auto &endpoint : graph[i]) {
			total_bandwidth += endpoint.bandwidth;
		}
	}
	return total_bandwidth;
}

unsigned long GetTotalResidualBandwidth_Rob() {
	unsigned long total_bandwidth = 0;
	for (int i = 0; i < nodes.size(); i++) {
		for (int j = 0; j < nodes.size(); j++) {
			total_bandwidth += bw[i][j];
		}
	}
	return total_bandwidth;
}

unsigned long GetTotalNodeResidualCores_Rob() {
	long total_cores = 0;
	for (int i = 0; i < nodes.size(); ++i) {
		total_cores += nodes[i].residual_cores;
		//cout<<nodes[i].residual_cores<<" ";
	}
	//cout<<endl;
	//cout<<"total residual cores: "<<total_cores<<endl;
	return total_cores;
}


int GetShortestPathLength(const std::vector<int> &result) {
	int embedded_path_length = 0;
	const int kSequenceLength = result.size();
	int kSource = result[0];
	int kDestination = result[kSequenceLength - 1];
	//get shortest path between node A and B
	int shortest_path_length = ComputeShortestPath(kSource, kDestination)->size() - 1;
	return shortest_path_length;
}

int GetSolutionProvisionPathLength(const std::vector<int> &result) {
	int embedded_path_length = 0;
	const int kSequenceLength = result.size();
	int kSource = result[0];
	int kDestination = result[kSequenceLength - 1];
	//get provision path length
	for (int i = 0; i < kSequenceLength - 1; ++i) {
		embedded_path_length +=
				ComputeShortestPath(result[i], result[i + 1])->size() - 1;
	}
	return embedded_path_length;
}

#if 1 // to be used
void ComputeNetworkUtilization_rob(std::map<int, int> mcOrder,
		std::map<int, std::vector<std::vector<int> > > SC_solutions) {

	//cout<<"ComputeNetworkUtilization_rob"<<endl;
	const unsigned long kNetworkCapacity = GetTotalNetworkBandwidth();
	cout<<"kNetworkCapacity: "<<kNetworkCapacity<<endl;
	unsigned long bandwidth_usage = 0;
	for (int mcid = 0; mcid < SC_solutions.size(); mcid++) {

		int currmcid = mcOrder[mcid];

		if (SC_solutions[mcid].size() > 0) {
			//get the bandwidth usage of main-chain
			//unsigned long bandwidth_usage = GetBandwidthUsage(solutions[mcid], traffic_requests[mcid]);

			//get the bw usage in each link
			//linkUsage[mcid] = GetLinkUsage(solutions[mcid], traffic_requests[mcid]);

			//get the bandwidth usage of sub-chains
			for (int scid = 0; scid < SC_solutions[mcid].size(); scid++) {

				bandwidth_usage += GetBandwidthUsage(SC_solutions[mcid][scid], traffic_requests_subchains[currmcid][scid]);

				//get the bw usage in each link
				map<pair<int, int>, long> utils = GetLinkUsage(SC_solutions[mcid][scid], traffic_requests_subchains[currmcid][scid]);
				map<pair<int, int>, long>::iterator piit;
				for(piit=utils.begin(); piit!=utils.end(); piit++)
				{
					if((linkUsage[mcid].find(piit->first))==linkUsage[mcid].end())
						linkUsage[mcid][piit->first] = piit->second;
					else
						linkUsage[mcid][piit->first] += piit->second;
				}
			}

			//net_util.push_back(static_cast<double>(bandwidth_usage) / static_cast<double>(kNetworkCapacity));

		}else{
			//net_util.push_back(-1.0);
		}

	}
}
#endif


#if 0
void ComputeSolutionCosts_rob(std::map<int, std::vector<int> > solutions,
		std::map<int, std::vector<std::vector<int> > > SC_solutions) {

	ofstream outfile("deploypath");

	cout<<"ComputeSolutionCosts_rob......"<<solutions.size()<<endl;
	int current_time = traffic_requests[0].arrival_time;

	int nodeUtilsCopyID = 0;  //for nodeUtils
	int nodeUtilCopyNum = 10; //for nodeUtils
	for (int i = 0; i < solutions.size(); i++) {

		if (current_time != traffic_requests[i].arrival_time) {

			//for nodeUtils
			if(nodeUtilsCopyID < nodeUtilCopyNum)
			{
				vector<int> nodeutil;
				for (auto &n : nodes) {
					if (n.num_cores <= 0)
						nodeutil.push_back(0);
					else{
						int used_cores = n.num_cores - n.residual_cores;
						nodeutil.push_back(used_cores);
					}
				}
				nodeUtils[nodeUtilsCopyID] = nodeutil;
				nodeUtilsCopyID++;
			}

			//cout<<"stats active server of time slot "<<current_time<<endl;
			int active_servers = 0;
			for (auto &n : nodes) {
				if (n.num_cores <= 0)
					continue;
				int used_cores = n.num_cores - n.residual_cores;
				if (used_cores > 0)
					++active_servers;
				//printf("ts = %d, nodeID = %d, Used cores = %d, duration = %d\n",
				//		current_time, n.node_id, used_cores, traffic_requests[i - 1].duration);
			}
			num_active_servers.push_back(std::pair<int, int>(current_time, active_servers));

			//cout<<"deployed vNFs over all nodes"<<endl;
			int n_deployed = 0;
			for (int j = 0; j < deployed_mboxes.size(); ++j)
				n_deployed += deployed_mboxes[j].size();
			mbox_count.push_back(n_deployed);

			//cout<<"record the server stats info."<<endl;
			RefreshServerStats(current_time);

			//cout<<"release all resource for provision in time slot"<<endl;
			ReleaseAllResources();
			current_time = traffic_requests[i].arrival_time;
		}

		vector<int> current_solution = solutions[i];
		const int kLastIndex = static_cast<int>(current_solution.size()) - 1;
		int embedded_path_length = 0;

		double transmitTau = 1.0 * traffic_requests[i].dataVol / traffic_requests[i].min_bandwidth;
		double total_delay = transmitTau;

		resource resource_vector;
		for (auto &node : nodes)
			resource_vector.cpu_cores.push_back(node.residual_cores);

		//cout<<"for each vNFs in main chain"<<endl;
		for (int kk = 1; kk < current_solution.size(); kk++) {

			//cout<<"i: "<<i<<" kk: "<<kk<<" seq.size: "<<traffic_requests[i].middlebox_sequence.size()
			//		<<" kLastIndex: "<<kLastIndex<<endl;
			auto &m_box = middleboxes[traffic_requests[i].middlebox_sequence[kk - 1]];
			int current_node = current_solution[kk];
			int prev_node = current_solution[kk - 1];

			std::unique_ptr<std::vector<int> > path = ComputeShortestPath(prev_node, current_node);
			for(int tmpi=1; tmpi < path->size(); tmpi++){
				int node0 = path->at(tmpi-1);
				int node1 = path->at(tmpi);
				outfile<<"("<<node0<<"-"<<node1<<"),";
				total_delay += shortest_path[node0][node1];
			}

#if 0
			// a different calculation for latency
			//cout<<"Transit latency from prev to current node, one part of total_delay"<<endl;
			total_delay += GetCost_rob(prev_node, current_node, traffic_requests[i]);
#endif
			//cout<<"vNF processing time, one part of total_delay"<<endl;
			//cout<<"middlebox_sequence[kk - 1]: "<<traffic_requests[i].middlebox_sequence[kk - 1]<<endl;
			if (kk != 0 && kk != kLastIndex) {
				total_delay += middleboxes[traffic_requests[i].middlebox_sequence[kk-1]].processing_delay;
				total_delay += transmitTau;
			}

			//cout<<"Update the resource vector with any new middleboxes."<<endl;
			if (kk != kLastIndex && UsedMiddleboxIndex(current_node, m_box, traffic_requests[i]) == NIL) {
				resource_vector.cpu_cores[current_node] -= m_box.cpu_requirement;
			}

		}

		//cout<<"SLA violation cost."<<endl;
		double sla_cost = 0.0;
		sla_cost = total_delay - traffic_requests[i].max_delay;

		sla_costs.push_back(sla_cost);
		total_costs.push_back(total_delay);
		DEBUG("current traffic request = %d\n", i);

		//cout<<"update resource utilization for main-chain"<<endl;
		UpdateResources(&current_solution, traffic_requests[i]);

#if 0
		//cout<<"for sub-chains"<<endl;
		for(int tmpscid = 0; tmpscid < SC_solutions[i].size(); tmpscid++)
		{
			if(mcIDscIDmap[i].size() == 0)
				continue;

			int scid = mcIDscIDmap[i][tmpscid];
			current_solution = SC_solutions[i][tmpscid];

			//cout<<"for each vNFs in its SFC"<<endl;
			for (int kk = 1; kk < current_solution.size(); ++kk) {
				auto &m_box = middleboxes[traffic_requests_subchains[scid].middlebox_sequence[1]];
				int current_node = current_solution[kk];
				int prev_node = current_solution[kk - 1];

				// Update the resource vector with any new middleboxes.
				if (UsedMiddleboxIndex(current_node, m_box, traffic_requests[i]) == NIL) {
					resource_vector.cpu_cores[current_node] -= m_box.cpu_requirement;
				}
			}

			//cout<<"update resource utilization for sub-chain"<<endl;
			UpdateResourcesForSubchain(&current_solution, traffic_requests_subchains[scid]);
		}
#endif

		outfile<<endl;

	}

	//cout<<"processing for the last time slot"<<endl;
	//int K = traffic_requests.size() - 1;

	//cout<<"stats active server"<<endl;
	int active_servers = 0;
	for (auto &n : nodes) {
		if (n.num_cores <= 0)
			continue;
		int used_cores = n.num_cores - n.residual_cores;
		if (used_cores > 0)
			++active_servers;
		//printf("ts = %d, nodeID = %d, Used cores = %d, duration = %d\n",
		//		current_time, n.node_id, used_cores, traffic_requests[K].duration);
	}
	num_active_servers.push_back(std::pair<int, int>(current_time, active_servers));

	//cout<<"deployed vNFs over all nodes"<<endl;
	int n_deployed = 0;
	for (int j = 0; j < deployed_mboxes.size(); ++j)
		n_deployed += deployed_mboxes[j].size();
	mbox_count.push_back(n_deployed);

	//record the server stats info.
	RefreshServerStats(current_time);

	//release all resource
	ReleaseAllResources();

	outfile.close();
}

void ComputeAllStretches_rob(std::map<int, std::vector<int> > solutions) {

	//cout<<"ComputeAllStretches_rob......"<<solutions.size()<<endl;
	for (int mcid = 0; mcid < solutions.size(); mcid++) {

		//cout<<"mcid: "<<mcid<<endl;
		vector<int> current_solution = solutions[mcid];
		if(current_solution.size() > 0){
			//cout << "get shortest path between node A and B" << endl;
			int shortestPathLen = GetShortestPathLength(current_solution);

			if(shortestPathLen == 0){
				cout<<"shortestPathLen == 0, mcid: "<<mcid<<endl;
			}

			shortestPathLength.push_back(shortestPathLen);

			//cout << "get provision path length" << endl;
			int embedded_path_length = GetSolutionProvisionPathLength(current_solution);
			provisionPathLength.push_back(embedded_path_length);// added by rob


			//cout << "compute the strectch of provision path" << endl;
			double stretch = 1.0 * embedded_path_length / shortestPathLen;
			stretches.push_back(stretch);
		} else {
			shortestPathLength.push_back(-1.0);
			provisionPathLength.push_back(-1.0);// added by rob
			stretches.push_back(-1.0);
		}
	}
}

void ComputeNetworkUtilization_rob(std::map<int, std::vector<int> > solutions,
		std::map<int, std::vector<std::vector<int> > > SC_solutions) {

	//cout<<"ComputeNetworkUtilization_rob"<<endl;
	const unsigned long kNetworkCapacity = GetTotalNetworkBandwidth();
	cout<<"kNetworkCapacity: "<<kNetworkCapacity<<endl;
	for (int mcid = 0; mcid < solutions.size(); mcid++) {

		if (solutions[mcid].size() > 0) {
			//get the bandwidth usage of main-chain
			unsigned long bandwidth_usage = GetBandwidthUsage(solutions[mcid], traffic_requests[mcid]);

			//get the bw usage in each link
			linkUsage[mcid] = GetLinkUsage(solutions[mcid], traffic_requests[mcid]);

			//get the bandwidth usage of sub-chains
			for (int tmpi = 0; tmpi < SC_solutions[mcid].size(); tmpi++) {
				if(mcIDscIDmap[mcid].size() == 0)
					continue;

				int scid = mcIDscIDmap[mcid][tmpi];
				bandwidth_usage += GetBandwidthUsage(SC_solutions[mcid][tmpi], traffic_requests_subchains[scid]);

				//get the bw usage in each link
				map<pair<int, int>, long> utils = GetLinkUsage(SC_solutions[mcid][tmpi], traffic_requests_subchains[scid]);
				map<pair<int, int>, long>::iterator piit;
				for(piit=utils.begin(); piit!=utils.end(); piit++)
				{
					if((linkUsage[mcid].find(piit->first))==linkUsage[mcid].end())
						linkUsage[mcid][piit->first] = piit->second;
					else
						linkUsage[mcid][piit->first] += piit->second;
				}
			}

			net_util.push_back(static_cast<double>(bandwidth_usage) / static_cast<double>(kNetworkCapacity));

		}else{
			net_util.push_back(-1.0);
		}

	}
}

void ComputeKHops_rob(std::map<int, std::vector<int> > solutions) {
	//cout<<"ComputeKHops_rob"<<endl;
	ingress_k.resize(solutions.size());
	egress_k.resize(solutions.size());
	for (int i = 0; i < solutions.size(); ++i) {
		if (solutions[i].size() > 0) {
			int ingress = solutions[i].front();
			int egress = solutions[i].back();
			int ihops = 0, ehops = 0;

			for (int j = 1; j < solutions[i].size() - 1; ++j) {
				ihops += ComputeShortestPath(solutions[i][j-1], solutions[i][j])->size() - 1;
				ingress_k[i].push_back(ihops);
			}

			for (int j = solutions[i].size() - 2; j >= 1; --j) {
				ehops += ComputeShortestPath(solutions[i][j+1], solutions[i][j])->size() - 1;
				egress_k[i].push_back(ehops);
			}
		}
	}
}

void ComputeCloseness_rob(std::map<int, std::vector<int> > solutions) {
	//cout<<"ComputeCloseness_rob"<<endl;
	sol_closeness.resize(solutions.size());
	for (int mcid = 0; mcid < solutions.size(); mcid++) {
		for (int nodeid = 0; nodeid < solutions[mcid].size(); nodeid++) {
			sol_closeness[mcid].push_back(closeness[nodeid]);
		}
	}
}

void ComputeServicePoints_rob(std::map<int, std::vector<int> > solutions) {
	//cout<<"ComputeServicePoints_rob"<<endl;
	for (int mcid = 0; mcid < solutions.size(); mcid++) {
		vector<int> current_solution = solutions[mcid];
		std::set<int> S;
		for (int i = 1; i < current_solution.size() - 1; ++i) {
			S.insert(current_solution[i]);
		}
		num_service_points.push_back(S.size());
	}
}
#endif


#if 0
void ComputeSolutionCosts_CPLEX_rob(std::map<int, std::vector<int> > solutions,
		std::map<int, std::vector<std::vector<int> > > SC_solutions, string deploy_res_file) {

	vector<double> costs;
	string filename = deploy_res_file + ".cost.ts";
	cout<<filename<<endl;

	string strline;
	ifstream infile;
	infile.open(filename.c_str());
	while(getline(infile, strline))
	{
		//cout<<strline<<endl;
		istringstream ss(strline);
		double tmpnum = 0;
		while(ss>>tmpnum)
		{
			//cout<<tmpnum<<endl;
			costs.push_back(tmpnum);
		}
	}
	infile.close();

	int reqid = 0;
	filename = deploy_res_file + ".paths-MC";
	infile.clear();
	infile.open(filename.c_str());
	while(getline(infile, strline))
	{
		double delay = costs[reqid];
		int prepos = -1,curpos = -1;
		while((curpos = strline.find(",",prepos+1)) != -1){
			string str = strline.substr(prepos+1,curpos-prepos);
			int pos = str.find("-");
			string tmpstr1 = str.substr(1,pos-1);
			string tmpstr2 = str.substr(pos+1,str.length()-pos-1);
			int n1 = atoi(tmpstr1.c_str());
			int n2 = atoi(tmpstr2.c_str());
			delay += shortest_path[n1][n2];
			//delay -= 1; // correction for the plusing 1 ms in the cplex program
			prepos = curpos;
		}
		costs[reqid] = delay;
		reqid++;
	}
	infile.close();


	//cout<<"ComputeSolutionCosts_rob......"<<solutions.size()<<endl;
	int current_time = traffic_requests[0].arrival_time;


	int nodeUtilsCopyID = 0;
	int nodeUtilCopyNum = 10;
	for (int i = 0; i < solutions.size(); ++i) {

		if (current_time != traffic_requests[i].arrival_time) {

			//for nodeUtils
			if (nodeUtilsCopyID < nodeUtilCopyNum) {
				vector<int> nodeutil;
				for (auto &n : nodes) {
					if (n.num_cores <= 0)
						nodeutil.push_back(0);
					else {
						int used_cores = n.num_cores - n.residual_cores;
						nodeutil.push_back(used_cores);
					}
				}
				nodeUtils[nodeUtilsCopyID] = nodeutil;
				nodeUtilsCopyID++;
			}

			//cout<<"stats active server of time slot "<<current_time<<endl;
			int active_servers = 0;
			for (auto &n : nodes) {
				if (n.num_cores <= 0)
					continue;
				int used_cores = n.num_cores - n.residual_cores;
				if (used_cores > 0)
					++active_servers;
				//printf("ts = %d, nodeID = %d, Used cores = %d, duration = %d\n",
				//		current_time, n.node_id, used_cores, traffic_requests[i - 1].duration);
			}
			num_active_servers.push_back(std::pair<int, int>(current_time, active_servers));

			//cout<<"deployed vNFs over all nodes"<<endl;
			int n_deployed = 0;
			for (int j = 0; j < deployed_mboxes.size(); ++j)
				n_deployed += deployed_mboxes[j].size();
			mbox_count.push_back(n_deployed);

			//cout<<"record the server stats info."<<endl;
			RefreshServerStats(current_time);

			//cout<<"release all resource for provision in time slot"<<endl;
			ReleaseAllResources();
			current_time = traffic_requests[i].arrival_time;
		}

		vector<int> current_solution = solutions[i];
		const int kLastIndex = static_cast<int>(current_solution.size()) - 1;
		double total_delay = costs[i];

		//cout<<"SLA violation cost."<<endl;
		double sla_cost = 0.0;
		sla_cost = total_delay - traffic_requests[i].max_delay;

		sla_costs.push_back(sla_cost);
		total_costs.push_back(total_delay);
		DEBUG("current traffic request = %d\n", i);

		//cout<<"update resource utilization for main-chain"<<endl;
		UpdateResources(&current_solution, traffic_requests[i]);

		//cout<<"for sub-chains"<<endl;
		for(int tmpscid = 0; tmpscid < SC_solutions[i].size(); tmpscid++)
		{
			if(mcIDscIDmap[i].size() == 0)
				continue;

			int scid = mcIDscIDmap[i][tmpscid];
			current_solution = SC_solutions[i][tmpscid];

			//cout<<"update resource utilization for sub-chain"<<endl;
			UpdateResourcesForSubchain(&current_solution, traffic_requests_subchains[scid]);
		}

	}

	//cout<<"*** processing for the last time slot ***"<<endl;
	//int K = traffic_requests.size() - 1;

	//cout<<"stats active server"<<endl;
	int active_servers = 0;
	for (auto &n : nodes) {
		if (n.num_cores <= 0)
			continue;
		int used_cores = n.num_cores - n.residual_cores;
		if (used_cores > 0)
			++active_servers;
		//printf("ts = %d, nodeID = %d, Used cores = %d, duration = %d\n",
		//		current_time, n.node_id, used_cores, traffic_requests[K].duration);
	}
	num_active_servers.push_back(std::pair<int, int>(current_time, active_servers));

	//cout<<"deployed vNFs over all nodes"<<endl;
	int n_deployed = 0;
	for (int j = 0; j < deployed_mboxes.size(); ++j)
		n_deployed += deployed_mboxes[j].size();
	mbox_count.push_back(n_deployed);

	//record the server stats info.
	RefreshServerStats(current_time);

	//release all resource
	ReleaseAllResources();
}

void ComputeAllStretches_CPLEX_rob(std::map<int, std::vector<int> > solutions) {

	//cout<<"ComputeAllStretches_CPLEX_rob......"<<solutions.size()<<endl;
	for (int mcid = 0; mcid < solutions.size(); mcid++) {

		//cout<<"mcid: "<<mcid<<" "<<solutions[mcid].size()<<" "<<all_mcPaths[mcid].size()<<endl;
		vector<int> current_solution = solutions[mcid];
		if(current_solution.size() > 0 && all_mcPaths[mcid].size() > 0){

			//cout << "get shortest path between node A and B" << endl;
			int shortestPathLen = GetShortestPathLength(current_solution);

			if(shortestPathLen == 0){
				cout<<"shortestPathLen == 0, mcid: "<<mcid<<endl;
			}

			shortestPathLength.push_back(shortestPathLen);

			//cout << "get provision path length" << endl;
			int embedded_path_length = all_mcPaths[mcid].size();
			provisionPathLength.push_back(embedded_path_length);// added by rob

			//cout << "compute the strectch of provision path" << endl;
			double stretch = 1.0 * embedded_path_length / shortestPathLen;
			stretches.push_back(stretch);

		} else {
			shortestPathLength.push_back(-1.0);
			provisionPathLength.push_back(-1.0);// added by rob
			stretches.push_back(-1.0);
		}
	}
}

void ComputeNetworkUtilization_CPLEX_rob(std::map<int, std::vector<int> > solutions,
		std::map<int, std::vector<std::vector<int> > > SC_solutions) {
	//cout<<"ComputeNetworkUtilization_rob"<<endl;
	const unsigned long kNetworkCapacity = GetTotalNetworkBandwidth();
	cout<<"kNetworkCapacity: "<<kNetworkCapacity<<endl;
	for (int mcid = 0; mcid < solutions.size(); mcid++) {

		if (solutions[mcid].size() > 0 && all_mcPaths[mcid].size() > 0) {

			//cout<<"get the bandwidth usage of main-chain"<<endl;
			unsigned long bandwidth_usage = traffic_requests[mcid].min_bandwidth * all_mcPaths[mcid].size();
			//get the bw usage in each link
			linkUsage[mcid] = GetLinkUsage(solutions[mcid], traffic_requests[mcid]);

			//cout<<"get the bandwidth usage of sub-chains"<<endl;
			for (int tmpi = 0; tmpi < SC_solutions[mcid].size(); tmpi++) {
				if(mcIDscIDmap[mcid].size() == 0 || all_scPaths[mcid].size() == 0)
					continue;

				int scid = mcIDscIDmap[mcid][tmpi];
				bandwidth_usage += traffic_requests[mcid].min_bandwidth * all_scPaths[mcid][tmpi].size();


				//get the bw usage in each link
				map<pair<int, int>, long> utils = GetLinkUsage(SC_solutions[mcid][tmpi], traffic_requests_subchains[scid]);
				map<pair<int, int>, long>::iterator piit;
				for(piit=utils.begin(); piit!=utils.end(); piit++)
				{
					if((linkUsage[mcid].find(piit->first))==linkUsage[mcid].end())
						linkUsage[mcid][piit->first] = piit->second;
					else
						linkUsage[mcid][piit->first] += piit->second;
				}
			}

			net_util.push_back(static_cast<double>(bandwidth_usage) / static_cast<double>(kNetworkCapacity));

		} else {
			net_util.push_back(-1.0);
		}
	}
	cout<<"finished ComputeNetworkUtilization_CPLEX_rob"<<endl;
}

void ComputeKHops_CPLEX_rob(std::map<int, std::vector<int> > solutions) {
	//cout<<"ComputeKHops_rob"<<endl;
	ingress_k.resize(solutions.size());
	egress_k.resize(solutions.size());
	for (int i = 0; i < solutions.size(); ++i) {
		if (solutions[i].size() > 0) {
			int ingress = solutions[i].front();
			int egress = solutions[i].back();
			int ihops = 0, ehops = 0;

			for (int j = 1; j < solutions[i].size() - 1; ++j) {
				ihops += ComputeShortestPath(solutions[i][j-1], solutions[i][j])->size() - 1;
				ingress_k[i].push_back(ihops);
			}

			for (int j = solutions[i].size() - 2; j >= 1; --j) {
				ehops += ComputeShortestPath(solutions[i][j+1], solutions[i][j])->size() - 1;
				egress_k[i].push_back(ehops);
			}
		}
	}
}
void ProcessCostLogs(const std::string &output_file_prefix) {
	const std::string kCostTsFileName = output_file_prefix + ".cost.ts";
	const std::string reqSummaryFile = output_file_prefix + ".reqSummary.ts";

	FILE *cost_ts_file = fopen(kCostTsFileName.c_str(), "w");
	FILE *summ_ts_file = fopen(reqSummaryFile.c_str(), "w");
	std::vector<double> cost_ts_data;
	//cout << "Log time series data for cost." << endl;
	int current_time = traffic_requests[0].arrival_time;
	double current_cost = 0.0;
	double current_sla_cost = 0.0;
	int t = 0;
	int reqNumTs = 0;
	for (int i = 0; i < all_mcResults.size(); ++i) {
		if (current_time != traffic_requests[i].arrival_time) {
			cost_ts_data.push_back(current_cost / reqNumTs);
			fprintf(cost_ts_file, "%d %lf %lf\n", current_time, current_cost / reqNumTs, current_sla_cost / reqNumTs);
			fprintf(summ_ts_file, "%d\n", reqNumTs);

			cout<<reqNumTs<<endl;

			current_time = traffic_requests[i].arrival_time;
			current_cost = current_sla_cost = 0.0;
			t++;
			reqNumTs = 0;
		}
		if(all_mcResults[i].size()>0){
			reqNumTs++;
			current_cost += total_costs[i];
			current_sla_cost += sla_costs[i];
		}
	}
	cost_ts_data.push_back(current_cost / reqNumTs);
	fprintf(cost_ts_file, "%d %lf %lf\n", current_time, current_cost / reqNumTs, current_sla_cost / reqNumTs);
	fclose(cost_ts_file);
	fclose(summ_ts_file);

	/*
	//cout << "Log mean, 5th, and 95th percentile of the total cost." << endl;
	const std::string kCostSummaryFileName = output_file_prefix + ".cost.summary";
	double mean_cost = GetMean(cost_ts_data);
	double fifth_percentile_cost = GetNthPercentile(cost_ts_data, 5);
	double ninety_fifth_percentile_cost = GetNthPercentile(cost_ts_data, 95);
	FILE *cost_summary_file = fopen(kCostSummaryFileName.c_str(), "w");
	fprintf(cost_summary_file, "%lf %lf %lf\n", mean_cost, fifth_percentile_cost, ninety_fifth_percentile_cost);
	fclose(cost_summary_file);
	*/

}

void ProcessStretchLogs(const std::string &output_file_prefix) {
	string cdfFileName = output_file_prefix + ".stretch.cdf";
	string sumFileName = output_file_prefix + ".stretch.summary";
	string serviceLenFile = output_file_prefix + ".service.lengths";
	FILE *stretch_file = fopen(cdfFileName.c_str(), "w");
	FILE *stretch_summary_file = fopen(sumFileName.c_str(), "w");
	FILE *serviceLen_file = fopen(serviceLenFile.c_str(), "w");

	int current_time = traffic_requests[0].arrival_time;
	vector<double> current_stretches;
	int shortestPathLen = 0;
	int embededPathlen = 0;
	double stretch = 0.0;
	int sfcLen = 0;
	int mcLen = 0;
	int reqNumTs = 0;
	int rejectNum = 0;

	for (int i = 0; i < traffic_requests.size(); ++i) {

		//cout<<"i: "<<i<<endl;
		if (current_time != traffic_requests[i].arrival_time) {

			//cout<<"print the stretch CDF"<<endl;
			std::vector<std::pair<double, double> > cdf = GetCDF(current_stretches);
			for (int i = 0; i < cdf.size(); ++i) {
				fprintf(stretch_file, "%lf %lf\n", cdf[i].first, cdf[i].second);
			}

			//cout<<"print the stretch percentile point"<<endl;
			double mean_stretch = GetMean(current_stretches);
			double first_percentile = GetNthPercentile(current_stretches, 1);
			double ninety_ninth_percentile = GetNthPercentile(current_stretches, 99);
			fprintf(stretch_summary_file, "%lf %lf %lf\n", mean_stretch,
					first_percentile, ninety_ninth_percentile);

			//cout<<"print the stretch, sfc length, ..."<<endl;
			fprintf(serviceLen_file, "%d %lf %lf %lf %lf %lf\n", current_time,
					double(sfcLen) / reqNumTs, double(mcLen) / reqNumTs,
					double(shortestPathLen) / reqNumTs,
					double(embededPathlen) / reqNumTs, stretch / reqNumTs);

			current_stretches.clear();
			current_time = traffic_requests[i].arrival_time;
			sfcLen = 0;
			mcLen = 0;
			shortestPathLen = 0;
			embededPathlen = 0;
			stretch = 0;
			reqNumTs = 0;
		}

		if (stretches[i] > 0) {
			current_stretches.push_back(stretches[i]);
			mcLen += traffic_requests[i].middlebox_sequence.size();
			sfcLen += reqSFCLen[i];
			shortestPathLen += shortestPathLength[i];
			embededPathlen += provisionPathLength[i];
			stretch += stretches[i];
			reqNumTs++;
		}
	}
	cout<<"print the last round results"<<endl;
	// print the stretch CDF
	std::vector<std::pair<double, double> > cdf = GetCDF(current_stretches);
	for (int i = 0; i < cdf.size(); ++i) {
		fprintf(stretch_file, "%lf %lf\n", cdf[i].first, cdf[i].second);
	}

	// print the stretch percentile point
	double mean_stretch = GetMean(current_stretches);
	double first_percentile = GetNthPercentile(current_stretches, 1);
	double ninety_ninth_percentile = GetNthPercentile(current_stretches, 99);
	fprintf(stretch_summary_file, "%lf %lf %lf\n", mean_stretch,
			first_percentile, ninety_ninth_percentile);

	//print the stretch, sfc length, ...
	fprintf(serviceLen_file, "%d %lf %lf %lf %lf %lf\n", current_time,
			double(sfcLen) / reqNumTs, double(mcLen) / reqNumTs,
			double(shortestPathLen) / reqNumTs,
			double(embededPathlen) / reqNumTs, stretch / reqNumTs);

	fclose(stretch_file);
	fclose(stretch_summary_file);
	fclose(serviceLen_file);
}

void ProcessNetUtilizationLogs(const std::string &output_file_prefix) {
	cout<<"ProcessNetUtilizationLogs"<<endl;
	// Write time series data for utilization.
	const std::string netUtilTsFileName = output_file_prefix + ".netutil.ts";
	string linkUtilTsfile = output_file_prefix + ".linkUtils.ts";
	FILE *netutil_ts_file = fopen(netUtilTsFileName.c_str(), "w");
	FILE *linkUtil_ts_file = fopen(linkUtilTsfile.c_str(), "w");

	int current_time = traffic_requests[0].arrival_time;
	double current_util = 0.0;
	std::vector<double> netutil_ts_data;
	map<pair<int, int>, long> tmpLinkUtil;
	int linkUtilCopyNum = 10;
	int linkUtilCopyID = 0;
	for (int i = 0; i < traffic_requests.size(); ++i) {
		if (current_time != traffic_requests[i].arrival_time) {
			netutil_ts_data.push_back(current_util);
			fprintf(netutil_ts_file, "%d %lf\n", current_time, current_util);
			current_time = traffic_requests[i].arrival_time;
			current_util = 0.0;

			if (tmpLinkUtil.size() > 0 && linkUtilCopyID < linkUtilCopyNum) {
				map<pair<int, int>, long>::iterator piit;
				std::vector<long> tmpvec;
				for (piit = tmpLinkUtil.begin(); piit != tmpLinkUtil.end(); piit++) {
					tmpvec.push_back(piit->second);
				}
				sort(tmpvec.begin(), tmpvec.end());
				for(int tmpi = tmpvec.size()-1; tmpi>=0; tmpi--)
				{
					fprintf(linkUtil_ts_file, "%lu ", tmpvec[tmpi]);
				}

				for(int tmpi = 0; tmpi< 200-tmpLinkUtil.size(); tmpi++)
					fprintf(linkUtil_ts_file, "0 ");
				fprintf(linkUtil_ts_file, "\n");
				linkUtilCopyID++;
			}

			tmpLinkUtil.clear();
		}

		if (net_util[i] > 0)
			current_util += net_util[i];

		if(linkUsage[i].size()>0 && linkUtilCopyID < linkUtilCopyNum)
		{
			map<pair<int, int>, long>::iterator piit;
			for(piit=linkUsage[i].begin(); piit!= linkUsage[i].end(); piit++)
			{
				pair<int, int> _link = piit->first;
				pair<int, int> _link_dual = make_pair(_link.second, _link.first);
				if((tmpLinkUtil.find(_link)) == tmpLinkUtil.end()){
					if((tmpLinkUtil.find(_link_dual)) == tmpLinkUtil.end()){
						tmpLinkUtil[piit->first] = piit->second;
					}
					else
						tmpLinkUtil[_link_dual] += piit->second;
				}
				else
					tmpLinkUtil[_link] += piit->second;
			}
		}

	}
	//to print info. for the last round
	netutil_ts_data.push_back(current_util);
	fprintf(netutil_ts_file, "%d %lf\n", current_time, current_util);



	fclose(netutil_ts_file);
	fclose(linkUtil_ts_file);

	/*
	// Write mean, 5th and 95th percentile of this utilization data to file.
	double mean_util = GetMean(netutil_ts_data);
	double fifth_percentile_util = GetNthPercentile(netutil_ts_data, 5);
	double ninety_fifth_percentile_util = GetNthPercentile(netutil_ts_data, 95);
	const std::string summaryFileName = output_file_prefix + ".netutil.summary";
	FILE *netutil_summary_file = fopen(summaryFileName.c_str(), "w");
	fprintf(netutil_summary_file, "%lf %lf %lf\n", mean_util,
			fifth_percentile_util, ninety_fifth_percentile_util);
	fclose(netutil_summary_file);
	*/
}

void ProcessServerUtilizationLogs(const std::string &output_file_prefix) {
	cout<<"ProcessServerUtilizationLogs"<<endl;
	// Process utilization data. Also derive fragmentation data from utilization
	// data: fragmentation = 1 - utilization.
	string utilTsFileName = output_file_prefix + ".serverutil.ts";
	string fragmentationFileName = output_file_prefix + ".serverfrag.ts";
	string utilPerNodeTsFile = output_file_prefix + ".nodeUtils.ts";
	FILE *util_ts_file = fopen(utilTsFileName.c_str(), "w");
	FILE *fragmentation_ts_file = fopen(fragmentationFileName.c_str(), "w");
	FILE *nodeUtils_ts_file = fopen(utilPerNodeTsFile.c_str(), "w");

	for(int tmpi = 0; tmpi < nodeUtils.size(); tmpi++)
	{
		vector<int> tmpvec = nodeUtils[tmpi];
		for(int tmpj = 0; tmpj < tmpvec.size(); tmpj++)
			fprintf(nodeUtils_ts_file, "%d ", tmpvec[tmpj]);
		fprintf(nodeUtils_ts_file, "\n");
	}


	int current_time = traffic_requests[0].arrival_time;
	std::vector<double> util_data;
	//stats.server_stats.emplace_back(INF, NIL, INF);
	std::vector<std::vector<double> > per_server_util;
	per_server_util.resize(graph.size());
	for (auto &server_stat : stats.server_stats) {

		//get per-server-stats
		if (server_stat.server_id != NIL) {
			if (fabs(server_stat.utilization - 0.0) > EPS) {
				per_server_util[server_stat.server_id].push_back(server_stat.utilization);
			}
		}

		if (current_time != server_stat.timestamp) {
			double mean_util = GetMean(util_data);
			double fifth_percentile_util = GetNthPercentile(util_data, 5);
			double ninety_fifth_percentile_util = GetNthPercentile(util_data, 95);
			fprintf(util_ts_file, "%d %lf %lf %lf\n", current_time, mean_util,
					fifth_percentile_util, ninety_fifth_percentile_util);
			double mean_fragmentation = 1 - mean_util;
			double fifth_percentile_fragmentation = 1 - fifth_percentile_util;
			double ninety_fifth_percentile_fragmentation = 1 - ninety_fifth_percentile_util;
			fprintf(fragmentation_ts_file, "%d %lf %lf %lf\n", current_time,
					mean_fragmentation, fifth_percentile_fragmentation,
					ninety_fifth_percentile_fragmentation);
			current_time = server_stat.timestamp;
			util_data.clear();
		}
		if (fabs(server_stat.utilization - 0.0) > EPS) {
			util_data.push_back(server_stat.utilization);
		}
	}

	// print the last round results
	double mean_util = GetMean(util_data);
	double fifth_percentile_util = GetNthPercentile(util_data, 5);
	double ninety_fifth_percentile_util = GetNthPercentile(util_data, 95);
	fprintf(util_ts_file, "%d %lf %lf %lf\n", current_time, mean_util,
			fifth_percentile_util, ninety_fifth_percentile_util);
	double mean_fragmentation = 1 - mean_util;
	double fifth_percentile_fragmentation = 1 - fifth_percentile_util;
	double ninety_fifth_percentile_fragmentation = 1 - ninety_fifth_percentile_util;
	fprintf(fragmentation_ts_file, "%d %lf %lf %lf\n", current_time,
			mean_fragmentation, fifth_percentile_fragmentation,
			ninety_fifth_percentile_fragmentation);

	fclose(util_ts_file);
	fclose(fragmentation_ts_file);
	fclose(nodeUtils_ts_file);

	// Process per server utilization data.
	const std::string perServerUtilFileName = output_file_prefix + ".per_server_util";
	FILE *per_server_util_file = fopen(perServerUtilFileName.c_str(), "w");
	std::vector<double> mean_util_data;
	for (int i = 0; i < per_server_util.size(); ++i) {
		if (per_server_util[i].size() > 0) {
			double mean_util = GetMean(per_server_util[i]);
			double fifth_percentile_util = GetNthPercentile(per_server_util[i], 5);
			double ninety_fifth_percentile_util = GetNthPercentile(per_server_util[i], 95);
			fprintf(per_server_util_file, "%lf %lf %lf\n", mean_util,
					fifth_percentile_util, ninety_fifth_percentile_util);
			mean_util_data.push_back(mean_util);
		}
	}
	fclose(per_server_util_file);

	// Process CDF of mean server utilization.
	const std::string kServerUtilCdfFile = output_file_prefix + ".sutil.cdf";
	FILE *server_util_cdf_file = fopen(kServerUtilCdfFile.c_str(), "w");
	std::vector<std::pair<double, double> > util_cdf = GetCDF(mean_util_data);
	for (auto &cdf_data : util_cdf) {
		fprintf(server_util_cdf_file, "%lf %lf\n", cdf_data.first,cdf_data.second);
	}
	fclose(server_util_cdf_file);
}

void ProcessKHopsLogs(const std::string &output_file_prefix) {
	cout<<"ProcessKHopsLogs"<<endl;
	std::vector<int> ihops, ehops;
	const std::string ingressKHopsFileName = output_file_prefix + ".ingress_k.cdf";
	const std::string egressKHopsFileName = output_file_prefix + ".egress_k.cdf";
	FILE *ingress_k_file = fopen(ingressKHopsFileName.c_str(), "w");
	FILE *egress_k_file = fopen(egressKHopsFileName.c_str(), "w");

	for (int i = 0; i < ingress_k.size(); ++i) {
		for (auto &elem : ingress_k[i])
			ihops.push_back(elem);
		for (auto &elem : egress_k[i])
			ehops.push_back(elem);
	}
	std::vector<std::pair<int, double> > ingress_k_cdf = GetCDF(ihops);
	std::vector<std::pair<int, double> > egress_k_cdf = GetCDF(ehops);
	for (auto &cdf : ingress_k_cdf) {
		fprintf(ingress_k_file, "%d %lf\n", cdf.first, cdf.second);
	}
	for (auto &cdf : egress_k_cdf) {
		fprintf(egress_k_file, "%d %lf\n", cdf.first, cdf.second);
	}

	fclose(ingress_k_file);
	fclose(egress_k_file);
}

void ProcessActiveServerLogs(const std::string &output_file_prefix) {
	const std::string kActiveServerLogFile = output_file_prefix
			+ ".active_server.ts";
	FILE *active_server_log = fopen(kActiveServerLogFile.c_str(), "w");
	for (auto &ts_data : num_active_servers) {
		fprintf(active_server_log, "%d %d\n", ts_data.first, ts_data.second);
	}
}

void ProcessMboxRatio(const std::string &output_file_prefix) {
	int traffic_count = 0;
	int current_time = traffic_requests[0].arrival_time;
	const std::string kMboxRatioFileName = output_file_prefix + ".mbox.ratio";
	FILE *mbox_ratio_file = fopen(kMboxRatioFileName.c_str(), "w");
	const int kMboxSeqSize = traffic_requests[0].middlebox_sequence.size();
	for (int i = 0; i < traffic_requests.size(); ++i) {
		if (current_time != traffic_requests[i].arrival_time) {
			int nmbox = mbox_count.front();
			mbox_count.pop_front();
			fprintf(mbox_ratio_file, "%d %d %lu\n", current_time, nmbox,
					traffic_requests[i].middlebox_sequence.size()
							* traffic_count);
			traffic_count = 0;
			current_time = traffic_requests[i].arrival_time;
		}
		++traffic_count;
	}
	if (!mbox_count.empty()) {
		fprintf(mbox_ratio_file, "%d %d %d\n", current_time, mbox_count.front(),
				kMboxSeqSize * traffic_count);
	}
	fclose(mbox_ratio_file);
}

void ProcessServicePointLogs(const std::string &output_file_prefix) {
	const std::string kServicePointLogFile = output_file_prefix
			+ ".service_points";
	FILE *service_point_log = fopen(kServicePointLogFile.c_str(), "w");
	std::vector<std::pair<int, double> > cdf = GetCDF(num_service_points);
	for (auto &cdf_element : cdf) {
		fprintf(service_point_log, "%d %lf\n", cdf_element.first,
				cdf_element.second);
	}
	fclose(service_point_log);
}

void ProcessClosenessLogs(const std::string &output_file_prefix) {
	const std::string kClosenessLogFile = output_file_prefix + ".closeness.cdf";
	FILE *closeness_log = fopen(kClosenessLogFile.c_str(), "w");
	std::vector<double> data;
	for (int i = 0; i < sol_closeness.size(); ++i) {
		for (int j = 0; j < sol_closeness[i].size(); ++j) {
			data.push_back(sol_closeness[i][j]);
		}
	}
	std::vector<std::pair<double, double> > cdf = GetCDF(data);
	for (int i = 0; i < cdf.size(); ++i) {
		fprintf(closeness_log, "%lf %lf\n", cdf[i].first, cdf[i].second);
	}
	fclose(closeness_log);
}

std::vector<int> CplexComputePath(
		const std::vector<std::pair<int, int> > &edges,
		const std::vector<int> sequence) {
	int source = sequence.front();
	int destination = sequence.back();
	std::vector<std::vector<int> > adj;
	adj.resize(graph.size());
	std::vector<int> indeg(graph.size(), 0);
	std::vector<int> outdeg(graph.size(), 0);
	for (auto &edge : edges) {
		DEBUG("(%d, %d)\n", edge.first, edge.second);
		adj[edge.first].push_back(edge.second);
		++outdeg[edge.first];
		++indeg[edge.second];
	}
	for (int i = 0; i < graph.size(); ++i) {
		if ((indeg[i] + outdeg[i]) != 0) {
			if (i != source && i != destination) {
				assert(indeg[i] == outdeg[i]);
			} else {
				DEBUG("node = %d, indeg = %d, outdeg = %d\n", i, indeg[i], outdeg[i]);
				assert(abs(indeg[i] - outdeg[i]) == 1);
			}
		}
	}
	std::stack<int> s;
	std::vector<int> path;
	int current_node = source;
	while (true) {
		if (adj[current_node].empty()) {
			path.push_back(current_node);
			if (!s.empty()) {
				current_node = s.top();
				s.pop();
			} else
				break;
		} else {
			s.push(current_node);
			int neighbor = adj[current_node].back();
			adj[current_node].pop_back();
			current_node = neighbor;
		}
	}
	std::reverse(path.begin(), path.end());
	return path;
}

#endif

#endif  // MIDDLEBOX_PLACEMENT_SRC_UTIL_H_
