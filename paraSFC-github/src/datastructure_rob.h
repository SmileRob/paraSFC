#ifndef MIDDLEBOX_PLACEMENT_SRC_DATASTRUCTURE_H_
#define MIDDLEBOX_PLACEMENT_SRC_DATASTRUCTURE_H_

#include <list>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <map>
#include <set>
#include <memory>
#include <stdlib.h>

using namespace std;

#define INF 99999999
#define MAXN 1000
#define NIL -1
#define EPS 1e-9

#define NUM_CORES_PER_SERVER 160
// #define SERVER_IDLE_ENERGY 0.0   // Kilo Watt
// #define SERVER_PEAK_ENERGY 0.135 // Kilo Watt
#define SERVER_IDLE_ENERGY 0.0805  // Kilo Watt
#define SERVER_PEAK_ENERGY 2.735   // Kilo Watt
#define POWER_CONSUMPTION_ONE_SERVER(cores)    \
  (SERVER_IDLE_ENERGY +                        \
   (SERVER_PEAK_ENERGY - SERVER_IDLE_ENERGY) * \
       ((1.0 * (cores)) / (1.0 * (NUM_CORES_PER_SERVER))))
#define PER_UNIT_ENERGY_PRICE 0.10
#define HW_MBOX_IDLE_ENERGY 1.1
#define HW_MBOX_PEAK_ENERGY 1.7
//#define HW_MBOX_TRAFFIC_CAPACITY 1000 // 1000Mbps
#define HW_MBOX_TRAFFIC_CAPACITY 50000  // 50Gbps
#define HW_MBOX_POWER_CONSUMPTION(traffic)       \
  (HW_MBOX_IDLE_ENERGY +                         \
   (HW_MBOX_PEAK_ENERGY - HW_MBOX_IDLE_ENERGY) * \
       ((1.0 * traffic) / (1.0 * (HW_MBOX_TRAFFIC_CAPACITY))))

#if 0
struct LINKNODE{
	int vnfID;
	int forkNum;
	struct LINKNODE * pre;
	struct LINKNODE * next[10];
};
#endif

struct sfcnode{
	int vnfID;
	int preNum;
	int forkNum;
	int pre[10];
	int next[10];
};


struct middlebox {
  std::string middlebox_name;
  int cpu_requirement;
  int processing_delay;
  int processing_capacity;
  double deployment_cost;  // Hourly cost of running a middlebox.

  middlebox(const std::string &mb_name, const std::string &mb_cpu,
            const std::string &mb_delay, const std::string &mb_capacity,
            const std::string &mb_cost)
      : middlebox_name(mb_name),
        cpu_requirement(atoi(mb_cpu.c_str())),
        processing_delay(atoi(mb_delay.c_str())),
        processing_capacity(atoi(mb_capacity.c_str())),
        deployment_cost(atof(mb_cost.c_str())) {}

  std::string GetDebugString() {
    return "middlebox_name : " + middlebox_name + ", cpu_requirement : " +
           std::to_string(cpu_requirement) + ", processing_delay : " +
           std::to_string(processing_delay) + ", processing_capacity : " +
           std::to_string(processing_capacity) + ", deployment_cost : " +
           std::to_string(deployment_cost);
  }
};

struct middlebox_instance {
  const middlebox *m_box;
  long residual_capacity;
  middlebox_instance(const middlebox *m_box, long res_cap)
      : m_box(m_box), residual_capacity(res_cap) {}
};

struct traffic_request {
  int arrival_time;
  int duration;
  int source, destination;

  long dataVol; // data volume of this request, added by rob

  // Deprecated
  int sla_specification;

  // Minimum bandwidth demand in Kbps.
  int min_bandwidth;

  // Maximum delay according to the SLA.
  int max_delay;

  // Penalty for per 1ms guarantee violation.
  double delay_penalty;

  // ...... branching node ID for sub-chains, added by rob
  int branchNodeID;
  bool accepted;
  int mcID;
  int originalID;
  int seqSFCLen;
  int sid;
  double SFCdelay;

  std::vector<int> middlebox_sequence;
  std::vector<std::vector<int> > psfcpaths;
  traffic_request(const std::string &tr_arrival_time,
                  const std::string &tr_source,
				  const std::string &tr_dest,
                  const std::string &tr_min_bandwidth,
                  const std::string &tr_max_delay,
                  const std::string &tr_data_volume,
				  const int tr_sid,
                  const std::vector<int> &tr_mbox_seq,
				  const std::vector<std::vector<int> > &tr_paths)
      : arrival_time(atoi(tr_arrival_time.c_str())),
        source(atoi(tr_source.c_str())),
        destination(atoi(tr_dest.c_str())),
        min_bandwidth(atoi(tr_min_bandwidth.c_str())),
        max_delay(atoi(tr_max_delay.c_str())),
		dataVol(atol(tr_data_volume.c_str())),   //added by rob
		sid(tr_sid),
        middlebox_sequence(tr_mbox_seq),
  	  	psfcpaths(tr_paths){}

  std::string GetDebugString() {
    std::string seq_string;
    for (auto value : middlebox_sequence) {
      seq_string += std::to_string(value) + " ";
    }
    return "arrival_time = " + std::to_string(arrival_time) + ", duration = " +
           std::to_string(duration) + ", source : " + std::to_string(source) +
           ", destination : " + std::to_string(destination) +
           ", min_bandwidth : " + std::to_string(min_bandwidth) +
           ", max_delay : " + std::to_string(max_delay) + ", delay_penalty : " +
           std::to_string(delay_penalty) + ", middlebox_sequence_length : " +
           std::to_string(middlebox_sequence.size()) +
           ", middlebox_sequence: " + seq_string;
  }

};

struct node {
  int node_id;
  int num_cores;
  int residual_cores;
  std::string GetDebugString() {
    return "node_id : " + std::to_string(node_id) + ", num_cores : " +
           std::to_string(num_cores);
  }
};

struct edge_endpoint {
  node *u;
  long bandwidth;
  int delay;
  long residual_bandwidth;
  edge_endpoint(node *u_ptr, long bw, int del)
      : u(u_ptr), bandwidth(bw), delay(del), residual_bandwidth(bw) {}
};

struct resource {
  std::vector<int> cpu_cores;
};

// Statistics for each traffic embedding.
struct traffic_statistics {
  // Arrival time of the traffic request.
  int arrival_time;

  // Energy cost of embedding the traffic request.
  double cost;

  // Hop distance of each middlebox from the ingress and egress.
  std::vector<int> ingress_hops, egress_hops;

  // Ratio of embedded path length and shortest path route length
  double stretch;
  traffic_statistics(int a_time, double c) : arrival_time(a_time), cost(c) {}
};

struct server_statistics {
  int timestamp;
  int server_id;
  double utilization;
  server_statistics(int ts, int sid, double util)
      : timestamp(ts), server_id(sid), utilization(util) {}
};

// Statistics for the whole solution.
struct solution_statistics {
  // Start time of the simulation in nano seconds.
  unsigned long long start_time;

  // End time of the simulation in nano seconds.
  unsigned long long end_time;

  // The number of accepted and rejected embeddings.
  int num_accepted, num_rejected;
  std::vector<traffic_statistics> t_stats;
  std::vector<server_statistics> server_stats;
};

inline int GetNodeCount(const std::vector<std::vector<edge_endpoint>> &g) {
  return g.size();
}

inline int GetEdgeCount(const std::vector<std::vector<edge_endpoint>> &g) {
  int edge_count = 0;
  for (int i = 0; i < g.size(); ++i) {
    edge_count += g[i].size();
  }
  edge_count /= 2;
  return edge_count;
}

// added by rob
extern double bwfactor;//bw expansion factor
extern std::string shareVNF;
extern double topoScale;
extern int rejectflag;
//extern std::map<int, int> scID2mcID; //map<subchain ID, its associated mainchain ID>
extern std::map<int, std::vector<int>> all_mcResults; // map<mcID, provision results>
extern std::map<int, std::vector<std::vector<int>>> all_scResults; // map<scID, provision resutls>
extern std::map<int, std::vector<std::pair<int, int> > > all_mcPaths;// map<mcID, var[pathID] = path >
extern std::map<int, std::vector<std::vector<std::pair<int, int> > > > all_scPaths;// map<mcID, var[scID][pathID] = path >
extern std::map<int, std::vector<int>> mcIDscIDmap; //map<mcid, vector<scIDs>>
extern std::map<int, int> scIDmcIDmap; //map<scid, mcid>
extern std::vector<int> provisionPathLength; //var[mcid]: provision path length
extern std::vector<int> shortestPathLength; //var[mcid]: the length of the shortest path between A and B
extern std::vector<int> reqSFCLen; //var[mcid]: the length of SFCs of request mcid
extern std::map<int, std::map<int, int> > mcRank; //map<mcid, map<vnf, rank in MC> >
extern std::map<int, std::map<std::pair<int,int>, long> > linkUsage; //map<mcid, map<link, netUtil> >
extern std::map<int, std::vector<int> > nodeUtils; //map[timeslot] = resrc_util[nodeid]
extern std::map<int, std::vector<int> > pathmap;
extern std::map<int, std::set<int> > vnfnext;
extern std::map<int, std::set<int> > vnfpre;
extern std::map<int, int> copyNum;
extern std::map<int, std::pair<int, double> > crtpaths;

extern std::map<int,map<pair<int, int>, int> > mpmap;
extern std::vector<std::vector<int> > 					sSFC;
extern std::map<int, std::vector<std::vector<int> > > 	pSFC;
extern std::map<std::pair<int, int>, unsigned long>  				reservedEdges;
extern std::map<int, unsigned int> 								reservedCores;


// Global data structures;
extern std::vector<middlebox> middleboxes;
extern std::vector<traffic_request> traffic_requests;
extern std::map<int,std::vector<traffic_request> > traffic_requests_subchains; // for request's subchains
extern std::vector<node> nodes;
extern std::vector<std::vector<edge_endpoint>> graph;
extern std::vector<double> closeness;
extern std::vector<std::vector<middlebox_instance>> deployed_mboxes;
extern std::list<int> mbox_count;
extern std::map<std::pair<int, int>, std::unique_ptr<std::vector<int>>> path_cache;
extern std::map<std::pair<int, int>, std::unique_ptr<std::vector<int>>> path_cache_rob;
extern std::vector<double> deployment_costs, energy_costs, transit_costs, sla_costs, total_costs, stretches;
extern std::vector<double> e_cost_ts;
extern std::vector<std::vector<int>> ingress_k;
extern std::vector<std::vector<int>> egress_k;
extern std::vector<std::pair<int, int>> num_active_servers;
extern std::vector<std::vector<double>> sol_closeness;
extern std::vector<int> num_service_points;
extern std::vector<double> net_util;
extern solution_statistics stats;
extern double per_core_cost, per_bit_transit_cost;
extern double cost[MAXN][MAXN];
extern int pre[MAXN][MAXN];
extern int shortest_path[MAXN][MAXN], sp_pre[MAXN][MAXN]; //former: shortest delay from u to v, latter: through which gain shortest delay
extern int shortest_edge_path[MAXN][MAXN]; //path length (by edge or hop)
extern long bw[MAXN][MAXN];
extern int max_time;
extern middlebox fake_mbox;
extern std::vector<std::vector<int>> results;
extern std::vector<std::vector<int>> paths;
#endif  // MIDDLEBOX_PLACEMENT_SRC_DATASTRUCTURE_H_
