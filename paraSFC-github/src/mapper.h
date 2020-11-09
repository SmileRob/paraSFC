
const std::string kUsage =
		"./middleman "
				"--per_core_cost=<per_core_cost>\n\t--per_bit_transit_cost=<per_bit_transit"
				"_cost>\n\t--topology_file=<topology_file>\n\t"
				"--middlebox_spec_file=<middlebox_spec_file>\n\t--traffic_r"
				"equest_file=<traffic_request_file>\n\t--algorithm=<algorithm>";

map<int,map<pair<int, int>, int> > mpmap;
vector<vector<int> > 					sSFC;
map<int, vector<vector<int> > > 		pSFC;
std::map<std::pair<int, int>, unsigned long>  				reservedEdges;
std::map<int, unsigned int> 									reservedCores;

map<int, vector<int> > pathmap;
std::map<int, std::set<int> > vnfnext;
std::map<int, std::set<int> > vnfpre;
std::map<int, int> copyNum;
std::vector<middlebox> middleboxes; //var[mboxID]
std::vector<traffic_request> traffic_requests; 	//var[req ID]
std::map<int, std::vector<traffic_request> > traffic_requests_subchains; //var[req ID], added by rob
std::vector<node> nodes;
std::vector<std::vector<edge_endpoint>> graph;
std::vector<double> closeness;
std::vector<std::vector<middlebox_instance>> deployed_mboxes; // var[nodeID][mboxID]: deployed mboxes in each node
std::vector<double> deployment_costs, energy_costs, transit_costs, sla_costs, total_costs, stretches;
std::vector<double> e_cost_ts;
std::vector<std::vector<int>> ingress_k, egress_k;
std::vector<std::pair<int, int>> num_active_servers;
std::vector<std::vector<double>> sol_closeness;
std::list<int> mbox_count;
std::vector<int> num_service_points;
std::vector<double> net_util;
double per_core_cost, per_bit_transit_cost;
double cost[MAXN][MAXN];
int pre[MAXN][MAXN];
int shortest_path[MAXN][MAXN], sp_pre[MAXN][MAXN];
int shortest_edge_path[MAXN][MAXN];
long bw[MAXN][MAXN];
int max_time;
std::map<std::pair<int, int>, std::unique_ptr<std::vector<int>>>path_cache;
std::map<std::pair<int, int>, std::unique_ptr<std::vector<int>>>path_cache_rob; //added by rob!!!
solution_statistics stats;
std::vector<vector<int> > CurrentMC_results; //var[reqID][vnfID]
std::map<int, vector<vector<int> > > CurrentSC_results; //var[reqID][sc_index][vnfID], added by rob!!!
std::vector<std::vector<int>> results;
std::vector<std::vector<int>> paths;
std::vector<std::vector<int>> optMC_results;	//var[mcID][vnfID]
map<int, vector<vector<int>>> optSC_results;	//map<mcID, vector[sc_index][vnfID]>
std::map<int, int> opt_mcOrder; //map<processing order, mainchain ID>
std::map<int, vector<int> > opt_scOrder; //map<processing order, subchain ID>
std::map<int, std::map<int, int> > mcRank; //var[mcid][vnf]: vnf's rank in main-chain
std::map<int, std::map<std::pair<int,int>, long> > linkUsage; //map<mcid, map<link, netUtil> >
std::map<int, std::vector<int> > nodeUtils; //map[timeslot] = resrc_util[nodeid]
std::map<int, pair<int, double> > crtpaths;

map<int, vector<int>> all_mcResults; // map<mcID, provision results>
map<int, vector<vector<int>>> all_scResults; // map<scID, provision resutls>
map<int, vector<std::pair<int, int> > > all_mcPaths;// map<mcID, var[pathID] = path >
map<int, vector<vector<std::pair<int, int> > > > all_scPaths;// map<mcID, var[scID][pathID] = path >
std::map<int, int> mcOrder; //map<processing order, abs mainchain ID>
std::map<int, vector<int> > scOrder; //map<processing order, abs subchain ID>
map<int, vector<int>> mcIDscIDmap; //map<mcid, vector<scIDs>>
map<int, int> scIDmcIDmap; //map<scid, mcid>
std::vector<int> provisionPathLength; //var[mcid]: provision path length
std::vector<int> shortestPathLength; //var[mcid]: the length of the shortest path between A and B
std::vector<int> reqSFCLen; //var[mcid]: the length of SFCs of request mcid
double bwfactor = 99;//bw expansion factor
string shareVNF = "false";
double topoScale = 1.0;
int rejectflag;
int opt_delay;
int opt_fitness;
double opt_utilization;
int total_delay;
int total_fitness;
double total_utilization;
int worst_fitness;
int worstMCID;
int maxIterNum = 100;


string outPath;

middlebox fake_mbox("switch", "0", "0", TOSTRING(INF), "0.0");
