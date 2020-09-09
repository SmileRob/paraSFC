#ifndef MIDDLEBOX_PLACEMENT_SRC_CPLEX_H
#define MIDDLEBOX_PLACEMENT_SRC_CPLEX_H

#include "datastructure_rob.h"
#include "util_rob.h"
#include <string>
#include <cmath>
#include <climits>
#include <utility>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

typedef IloArray<IloIntVarArray> IloIntVar2dArray;
typedef IloArray<IloIntVar2dArray> IloIntVar3dArray;
typedef IloArray<IloIntVar3dArray> IloIntVar4dArray;
typedef IloArray<IloIntVar4dArray> IloIntVar5dArray;

typedef IloArray<IloIntArray> IloInt2dArray;
typedef IloArray<IloInt2dArray> IloInt3dArray;

typedef IloArray<IloExprArray> IloExpr2dArray;


void run_cplex(std::vector<traffic_request> _traffic_requests,
		map<int, std::vector<traffic_request> > _traffic_requests_SC,
		double &opex,
		map<int,std::vector<double> > &delay_breakdown,
		double &running_time,
		map<int, std::vector<int> > *scSequence,
		map<int, std::vector<std::pair<int, int> > > *scEdges,
		int *delays,
		std::vector<int> &utilization, map<int, int> &bandwidths,
		string topology_filename,
		string sSFCres_file,
		double bwscale) {

	int maxsfclen = 20;

	IloEnv env;
	try {
		cout << "declare the model and the solver" << endl;
		IloModel model(env);
		IloCplex cplex(model);
		cplex.setParam(IloCplex::DataCheck, 1);

		printf("Modeling Physical Network...\n");
		cout << "CPLEX Variable: counters" << endl;
		int kLinkCount = 0;
		int kServerCount = 0;
		int kVNFCount = middleboxes.size();
		int kTrafficCount = _traffic_requests.size();
		cout << "kTrafficCount: " << kTrafficCount << endl;
		int mcOffset = _traffic_requests[0].mcID;          //???
		map<int, unsigned long> sSFCbw;

		cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<sSFCres_file<<endl;
		if(sSFCres_file != ""){
			FILE *sSFCres = fopen(sSFCres_file.c_str(), "r");
			unsigned long tmpbw;
			for(int tmpi=0; tmpi<kTrafficCount; tmpi++){
				fscanf(sSFCres, "%lu", &tmpbw);
				sSFCbw[sSFCbw.size()] = tmpbw;
				cout<<tmpbw<<" ";
			}
			fclose(sSFCres);
		}
		cout<<endl;

		// read topology file and populate associated variables
		FILE *topology_file = fopen(topology_filename.c_str(), "r");

		// read switch and link count
		fscanf(topology_file, "%d %d", &kServerCount, &kLinkCount);
		cout<<"kServerCount & kLinkCount: "<<kServerCount<<" "<<kLinkCount<<endl;

		cout << "CPLEX Variable: vNF attributes" << endl;
		// R_func[v] = R+, denote the resource requirement of function f
		// tau_F[f] = R+, denote the processing time of function f
		int R_func[kVNFCount];
		double tau_F[kVNFCount];
		R_func[0] = 0;
		R_func[kVNFCount] = 0;
		tau_F[0] = 0;
		tau_F[kVNFCount] = 0;
		for (int f = 1; f < kVNFCount; f++) {
			R_func[f] = middleboxes[f].cpu_requirement;
			tau_F[f] = middleboxes[f].processing_delay;
		}

		// tau_T[r] = R+, denote the transmission latency of SFC request r
		// BW_req[r] = R+, denote the bandwidth requirement of SFC request r
		int tau_T[kTrafficCount];
		int BW_req[kTrafficCount];
		for (int r = 0; r < kTrafficCount; r++) {
			//tau_T[r] = _traffic_requests[r].dataVol / _traffic_requests[r].min_bandwidth;
			BW_req[r] = _traffic_requests[r].min_bandwidth;
		}

		// C_link[l] = R+, denote the bandwidth capacity of link l
		// tau_L[l] = R+, denote the transporting latency of link l
		int C_link[kLinkCount * 2];
		double tau_L[kLinkCount * 2];
		for (int l = 0; l < kLinkCount * 2; l++) {
			C_link[l] = 0;
			tau_L[l] = 0;
		}

		// C_node[v] = R+, denote the resource capacity of node v
		int C_node[kServerCount];
		for (int v = 0; v < kServerCount; v++) {
			C_node[v] = 0;
		}

		cout << "CPLEX Variable: indicators" << endl;
		//   I_l_in[l][v]: whether l is the incoming link of node v
		//  I_l_out[l][v]: whether l is the outgoing link of node v
		// I_src_mc[r][v]: whether v is the ingress node of request r
		// I_dst_mc[r][v]: whether v is the egress node of request r
		int  I_l_in[kLinkCount * 2][kServerCount];
		int I_l_out[kLinkCount * 2][kServerCount];
		int I_src_mc[kTrafficCount][kServerCount];
		int I_dst_mc[kTrafficCount][kServerCount];

		//initialization
		for (int l = 0; l < kLinkCount * 2; l++) {
			for (int v = 0; v < kServerCount; v++) {
				I_l_in[l][v] = 0;
				I_l_out[l][v] = 0;
			}
		}
		for (int r = 0; r < kTrafficCount; r++) {
			for (int v = 0; v < kServerCount; v++) {
				I_src_mc[r][v] = 0;
				I_dst_mc[r][v] = 0;
			}
		}

		cout << "CPLEX Variable: ID mappers" << endl;
		//int nodePair2link[kServerCount][kServerCount];
		map<int, pair<int, int> > link2nodePair;

		cout << "load node info" << endl;
		for (int i = 0; i < kServerCount; i++) {
			int aa, bb;
			fscanf(topology_file, "%d %d", &aa, &bb);

			C_node[i] = nodes[i].num_cores;
			cout<<C_node[i]<<" ";
		}
		//cout<<endl;

		cout << "read link info from file" << endl;
		for (int _l = 0; _l < kLinkCount; _l++) {
			int source, destination;
			unsigned long bandwidth;
			double latency;
			fscanf(topology_file, "%d %d %lu %lf", &source, &destination, &bandwidth, &latency);
			//source = source-1;
			//destination = destination -1;
			bandwidth = bandwidth*bwscale;
			int delay = latency;
			cout<<source<<" "<<destination<<" "<<bandwidth<<" "<<latency<<endl;

			C_link[_l] = bandwidth;
			C_link[_l + kLinkCount] = bandwidth;
			tau_L[_l] = delay;
			tau_L[_l + kLinkCount] = delay;

			//cout<<_l<<" & "<<_l+kLinkCount<<": "<<bandwidth<<endl;

			//nodePair2link[_u][_v] = _l;
			link2nodePair[_l] = make_pair(source, destination);
			//nodePair2link[_v][_u] = _l + kLinkCount;
			link2nodePair[_l + kLinkCount] = make_pair(destination, source);

			I_l_in[_l][source] = 1;
			I_l_out[_l][destination] = 1;
			I_l_in[_l + kLinkCount][destination] = 1;
			I_l_out[_l + kLinkCount][source] = 1;

		}
		//cout<<endl;
		fclose(topology_file);

		cout << "get the ingress and egress node info." << endl;
		for (int r = 0; r < _traffic_requests.size(); r++) {
			I_src_mc[r][_traffic_requests[r].source] = 1;
			I_dst_mc[r][_traffic_requests[r].destination] = 1;
			//cout<<_traffic_requests[r].source<<" -> "<<_traffic_requests[r].destination<<"; ";
		}
		//cout<<endl;

		cout << "CPLEX Variable: decision variables" << endl;
		IloIntVar3dArray alpha(env, kTrafficCount); //alpha[r][vnfID][v]
		for (int r = 0; r < kTrafficCount; r++) {
			int dim = _traffic_requests[r].middlebox_sequence.size();
			//int vnfNum = middleboxes.size();
			alpha[r] = IloIntVar2dArray(env, maxsfclen);
			for (int vnfID = 0; vnfID < maxsfclen; vnfID++) {
				alpha[r][vnfID] = IloIntVarArray(env, kServerCount, 0, 1);
			}
		}



		IloIntVar4dArray beta(env, kTrafficCount);//beta[r][vnfID+2][vnfID+2][l]
		for (int r = 0; r < kTrafficCount; r++) {
			//int vnfNum = middleboxes.size();
			beta[r] = IloIntVar3dArray(env, maxsfclen);
			for(int vnfID = 0; vnfID < maxsfclen; vnfID++){
				beta[r][vnfID] = IloIntVar2dArray(env, maxsfclen);
				for(int vnfID2 = 0; vnfID2 < maxsfclen; vnfID2++){
					beta[r][vnfID][vnfID2] = IloIntVarArray(env, kLinkCount * 2, 0, 1);
				}
			}
		}

		cout << "CPLEX Constraint: node constraints" << endl;
		// each vNF in a SFC deployed once
		for (int r = 0; r < kTrafficCount; r++) {
			std::vector<int> tmpvec = _traffic_requests[r].middlebox_sequence;
			for(int i = 1; i < tmpvec.size()-1; i++) {
				IloExpr sum(env);
				for (int v = 0; v < kServerCount; v++) {
					sum += alpha[r][i][v];
				}
				//model.add(sum == 1);
				char name[100];
				sprintf(name, "deployOnce@r%d-i%d", r, i);
				model.add(IloRange(env, 1, sum, 1, name));
			}

			int tmpsrc=_traffic_requests[r].source;
			int tmpdst=_traffic_requests[r].destination;

			if(1){
				IloExpr sum(env);
				sum += alpha[r][0][tmpsrc];
				char name[100];
				sprintf(name, "deployIngress@r%d", r);
				model.add(IloRange(env, 1, sum, 1, name));

			}

			if(1){
				int tmpi=tmpvec.size()-1;
				IloExpr sum(env);
				sum += alpha[r][tmpi][tmpdst];
				char name[100];
				sprintf(name, "deployEgress@r%d", r);
				model.add(IloRange(env, 1, sum, 1, name));

			}

			if(1){
				IloExpr sum(env);
				for (int v = 0; v < kServerCount; v++) {
					if(v==tmpsrc)
						continue;
					sum += alpha[r][0][v];
				}
				//model.add(sum == 1);
				char name[100];
				sprintf(name, "alpha0_nonIngress@r%d", r);
				model.add(IloRange(env, 0, sum, 0, name));

			}

			if(1){
				int tmpi=tmpvec.size()-1;
				IloExpr sum(env);
				for (int v = 0; v < kServerCount; v++) {
					if(v==tmpdst)
						continue;
					sum += alpha[r][tmpi][v];
				}
				//model.add(sum == 1);
				char name[100];
				sprintf(name, "alpha11_nonEgress@r%d", r);
				model.add(IloRange(env, 0, sum, 0, name));

			}

		}



		cout << "consumed resource can not exceed the node's capacity" << endl;
		for (int v = 0; v < kServerCount; v++) {
			IloExpr sum(env);
			for (int r = 0; r < kTrafficCount; r++) {
				std::vector<int> tmpvec = _traffic_requests[r].middlebox_sequence;
				for(int fid = 1; fid < tmpvec.size()-1; fid++) {
					int vnfID = tmpvec[fid];
					sum += alpha[r][fid][v]* R_func[vnfID];
				}
			}
			sum -= C_node[v];
			//model.add(sum <= 0);
			char name[100];
			sprintf(name, "nodeResourceConstraint@v%d", v);
			model.add(IloRange(env, -IloInfinity, sum, 0, name));
		}

		cout << "CPLEX Constraint: link bandwidth constraints" << endl;
		// consumed bandwidth can not exceed the link's capacity
		for (int l = 0; l < kLinkCount * 2; l++) {
			IloExpr sum(env);
			for (int r = 0; r < kTrafficCount; r++) {
				int scNum = _traffic_requests_SC[r].size();
				//int vnfNum = middleboxes.size();
				int ssfclen = _traffic_requests[r].middlebox_sequence.size();
				bool pairfilter[ssfclen][ssfclen] = {false};  //vnfID=0: ingress, vnfID=vnfNum-1: egress
				for(int cc = 0; cc < scNum; cc++){
					std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
					for (int i = 0; i < tmpvec.size()-1; i++) {
						int currvnf = tmpvec[i];
						int nextvnf = tmpvec[i+1];

						double resbw = 0.088;
						bool addbw = false;
						int sid=_traffic_requests[r].sid;
						if(sSFCres_file != "" && mpmap[sid].find( make_pair(currvnf, nextvnf)) != mpmap[sid].end() ){
							addbw = true;
						}

						int tmpbw = _traffic_requests[r].min_bandwidth;
						if(cc==0){
							sum += beta[r][currvnf][nextvnf][l] * tmpbw;
						} else if(!pairfilter[currvnf][nextvnf]){
							sum += beta[r][currvnf][nextvnf][l] * tmpbw;
							pairfilter[currvnf][nextvnf] = true;
						}else if (addbw){
							//if the link is used to transmit packet header
							//bw = ...;
							//the average packet size in data centers is around 724 bytes.
							//For TCP packets, the header only occupies 8.8% of the total size.
							tmpbw = resbw * tmpbw;
							sum += beta[r][currvnf][nextvnf][l] * tmpbw;
						}
					}
				}
			}
			sum -= C_link[l];
			//if(C_link[l] < 0) cout<<"<0 @ l="<<l<<": "<<C_link[l]<<endl;;
			//model.add(sum <= 0);
			char name[100];
			sprintf(name, "linkBWConstraint@l%d", l);
			model.add(IloRange(env, -IloInfinity, sum, 0, name));
		}
		//cout<<endl;


		if(sSFCres_file!=""){
			cout<<"Extra resource consumption constraints"<<endl;

			for (int r = 0; r < kTrafficCount; r++) {
				IloExpr sum(env);

				int scNum = _traffic_requests_SC[r].size();
				//int vnfNum = middleboxes.size();
				int ssfclen = _traffic_requests[r].middlebox_sequence.size();
				bool pairfilter[ssfclen][ssfclen] = {false};  //vnfID=0: ingress, vnfID=vnfNum-1: egress
				for(int cc = 0; cc < scNum; cc++) {
					std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
					for (int i = 0; i < tmpvec.size()-1; i++) {
						int currvnf = tmpvec[i];
						int nextvnf = tmpvec[i+1];

						double resbw = 0.088;
						bool addbw = false;
						int sid=_traffic_requests[r].sid;
						if(sSFCres_file != "" && mpmap[sid].find( make_pair(currvnf, nextvnf)) != mpmap[sid].end() ){
							addbw = true;
						}

						for (int l = 0; l < kLinkCount * 2; l++) {
							int tmpbw = _traffic_requests[r].min_bandwidth;
							if(cc==0){
								sum += beta[r][currvnf][nextvnf][l] * tmpbw;
							} else if(!pairfilter[currvnf][nextvnf]){
								sum += beta[r][currvnf][nextvnf][l] * tmpbw;
							}else if (addbw){
								//if the link is used to transmit packet header
								//bw = ...;
								//the average packet size in data centers is around 724 bytes.
								//For TCP packets, the header only occupies 8.8% of the total size.
								tmpbw = 0.088*tmpbw;
								tmpbw = resbw * tmpbw;
								sum += beta[r][currvnf][nextvnf][l] * tmpbw;
							}
						}
						pairfilter[currvnf][nextvnf] = true;
					}
				}
				sum -= sSFCbw[r];
				sum /= sSFCbw[r];

				//if(C_link[l] < 0) cout<<"<0 @ l="<<l<<": "<<C_link[l]<<endl;;
				//model.add(sum <= 0);
				char name[100];
				sprintf(name, "extraBWcontrol@r%d", r);
				model.add(IloRange(env, 0, sum, bwfactor, name));
			}
			//cout<<endl;
		}


#if 0
		cout << "connectivity guarantee for each path" << endl;
		for (int r = 0; r < kTrafficCount; r++) {
			int scNum = _traffic_requests_SC[r].size();
			int ssfclen = _traffic_requests[r].middlebox_sequence.size();
			for (int cc = 0; cc < scNum; cc++) {
				//IloExpr sum(env);
				std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;

				for (int v = 0; v < kServerCount; v++) {

					for (int i = 0; i < tmpvec.size()-1; i++) {
						int currvnf = tmpvec[i];
						int nextvnf = tmpvec[i+1];

						IloExpr sum(env);
						for (int l = 0; l < kLinkCount * 2; l++) {
							//sum -= beta[r][currvnf][nextvnf][l] * I_l_in[l][v];
							sum += beta[r][currvnf][nextvnf][l] * I_l_out[l][v];
						}
						if(nextvnf==ssfclen-1)
							sum += I_dst_mc[r][_traffic_requests[r].destination];

						sum -= alpha[r][currvnf][v];

						//model.add(sum == 0);
						char name[100];
						sprintf(name, "connectivity@r%d-vnf%d-%d", r, currvnf, nextvnf);
						model.add(IloRange(env, 0, sum, IloInfinity, name));

					}
				}

				//sum += (I_src_mc[r][v] - I_dst_mc[r][v]);
				//model.add(sum == 0);
				//char name[100];
				//sprintf(name, "in==out@r%d-v%d", r, v);
				//model.add(IloRange(env, 0, sum, 0, name));
			}
		}
#endif


#if 1///
		cout << "connectivity guarantee for each path" << endl;
		for (int r = 0; r < kTrafficCount; r++) {
			int scNum = _traffic_requests_SC[r].size();
			int ssfclen = _traffic_requests[r].middlebox_sequence.size();
			for (int cc = 0; cc < scNum; cc++) {
				//IloExpr sum(env);


				std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
				for (int i = 0; i < tmpvec.size()-1; i++) {
					int currvnf = tmpvec[i];
					int nextvnf = tmpvec[i+1];

				    IloExpr sum(env);
					for (int v = 0; v < kServerCount; v++) {
						for (int l = 0; l < kLinkCount * 2; l++) {
							sum += beta[r][currvnf][nextvnf][l] * I_l_in[l][v] * v;
							sum -= beta[r][currvnf][nextvnf][l] * I_l_out[l][v] * v;
						}
						sum -= alpha[r][currvnf][v]*v;
						sum += alpha[r][nextvnf][v]*v;
					}


					char name[100];
					sprintf(name, "connectivity@r%d-cc%d", r, cc);
					//model.add(IloRange(env, 1, sum, IloInfinity, name));
					model.add(IloRange(env, 0, sum, 0, name));
				}

			}
		}
#endif


#if 1
		cout << "inflow == outflow in each node" << endl;
		for (int v = 0; v < kServerCount; v++) {
			for (int r = 0; r < kTrafficCount; r++) {
				int scNum = _traffic_requests_SC[r].size();
				//int ssfclen = _traffic_requests[r].middlebox_sequence.size();
				for (int cc = 0; cc < scNum; cc++) {
					//IloExpr sum(env);

					int vnfNum = middleboxes.size();
					std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
					for (int i = 0; i < tmpvec.size()-1; i++) {
						int currvnf = tmpvec[i];
						int nextvnf = tmpvec[i+1];

						IloExpr sum(env);
						//sum = alpha[r][currvnf][v]-alpha[r][nextvnf][v];
						//for (int l = 0; l < kLinkCount * 2; l++) {
						//	sum -= beta[r][currvnf][nextvnf][l] * I_l_in[l][v];
						//}


						for (int l = 0; l < kLinkCount * 2; l++) {
							sum -= beta[r][currvnf][nextvnf][l] * I_l_in[l][v];
							sum += beta[r][currvnf][nextvnf][l] * I_l_out[l][v];
						}

						//if(currvnf == 0)
						//	sum += I_src_mc[r][v];
						//else
							sum += alpha[r][currvnf][v];

						//if(nextvnf == middleboxes.size()-1)
						//	sum -= I_dst_mc[r][v];
						//else
							sum -= alpha[r][nextvnf][v];


						//model.add(sum == 0);
						char name[100];
						sprintf(name, "in==out@r%d-v%d-vnf%d-%d", r, v, currvnf, nextvnf);
						model.add(IloRange(env, 0, sum, 0, name));

					}
					//sum += (I_src_mc[r][v] - I_dst_mc[r][v]);
					//model.add(sum == 0);
					//char name[100];
					//sprintf(name, "in==out@r%d-v%d", r, v);
					//model.add(IloRange(env, 0, sum, 0, name));


				}
			}
		}
#endif

		/*
		cout<<"consider ingress/egress nodes or vNF deployed nodes"<<endl;
		// at least one outgoing/incoming link for each ingress/egress node or vNF deployed node in mainchain
		for (int r = 0; r < kTrafficCount; r++) {
			for (int v = 0; v < kServerCount; v++) {

				// for each vNF deployed node
				// it has at least one incoming and one outgoing links
				if (I_src_mc[r][v] == 0 && I_dst_mc[r][v] == 0) {
					int scNum = _traffic_requests_SC[r].size();
					for (int cc = 0; cc < scNum; cc++){
						int vnfNum = middleboxes.size();
						std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
						int prevnf = 0;
						int curvnf = -1;
						for (int i = 0; i < tmpvec.size(); i++) {
							curvnf = tmpvec[i];

							// at least one incoming link in the deployed node
							IloExpr sumIn(env);
							for (int l = 0; l < kLinkCount * 2; l++) {
								sumIn += beta[r][prevnf][curvnf][l] * I_l_out[l][v];
							}

//							int prev = 0;
//							int curv = -1;
//							for(int tmpi = 0; tmpi < tmpvec.size(); tmpi++){
//								for (int l = 0; l < kLinkCount * 2; l++) {
//									sumIn += beta[r][prev][curv][l] * I_l_out[l][v];
//								}
//							}

							// The following is to force the deployed node lie in selected path
							sumIn -= alpha[r][curvnf][v];
							char nameIn[100];
							sprintf(nameIn, "atleastOneIn@r%d-v%d-i%d", r, v, curvnf);
							model.add(IloRange(env, -1, sumIn, 1, nameIn));

							prevnf = curvnf;
						}

						prevnf = tmpvec[0];
						curvnf = -1;
						for (int i = 1; i < tmpvec.size()+1; i++) {
							if(i < tmpvec.size())
								curvnf = tmpvec[i];
							else
								curvnf = vnfNum-1;

							// at least one outgoing link in the deployed node
							IloExpr sumOut(env);
							for (int l = 0; l < kLinkCount * 2; l++) {
								sumOut -= beta[r][prevnf][curvnf][l] * I_l_in[l][v];
							}
							// The following is to force the deployed node lie in selected path
							sumOut += alpha[r][prevnf][v];
							char nameOut[100];
							sprintf(nameOut, "atleastOneOut@r%d-v%d-i%d", r, v, prevnf);
							model.add(IloRange(env, -1, sumOut, 1, nameOut));

							prevnf = curvnf;
						}
					}
				}
			}
		}
		*/


		/*
		for (int l = 0; l < kLinkCount * 2; l++) {
			cout<<tau_L[l]<<" ";
		}
		cout<<endl;
		for (int f = 0; f < kVNFCount; f++) {
			cout<<R_func[f]<<" ";
		}
		cout<<endl;
		for (int f = 0; f < kVNFCount; f++) {
			cout<<tau_F[f]<<" ";
		}
		cout<<endl;
		*/

		cout << "Objective & Solution" << endl;
		// build the objective
		IloExpr objective(env);

		for (int r = 0; r < kTrafficCount; r++) {
			int scNum = _traffic_requests_SC[r].size();
			vector<int> ssfc = _traffic_requests[r].middlebox_sequence;
			for(int cc = 0; cc < scNum; cc++){

				int vnfNum = middleboxes.size();
				std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
				for (int i = 0; i < tmpvec.size()-1; i++) {
					int currvnf = tmpvec[i];
					int nextvnf = tmpvec[i+1];

					for (int l = 0; l < kLinkCount * 2; l++) {
						objective += beta[r][currvnf][nextvnf][l] * tau_L[l];
					}
				}

				for (int v = 0; v < kServerCount; v++) {
					std::vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
					for (int i = 1; i < tmpvec.size(); i++) {
						int fid = tmpvec[i];
						if(fid == ssfc.size()-1)
							continue;
						int vnfID = _traffic_requests[r].middlebox_sequence[fid];
						//int vnfID = tmpvec[i];
						//if(vnfID == middleboxes.size()-1)
						//	continue;
						objective += alpha[r][fid][v] * (tau_F[vnfID]);
					}
				}
			}

		}

		model.add(objective >= 0);
		model.add(IloMinimize(env, objective));

		/*
		 cout << "R_func: " << endl;
		 for (int tmpi = 0; tmpi < kVNFCount; tmpi++) {
		 cout << R_func[tmpi] << " ";
		 }
		 cout << endl;

		 cout << "C_link: " << endl;
		 for (int tmpi = 0; tmpi < kLinkCount; tmpi++) {
		 cout << C_link[tmpi] << " ";
		 }
		 cout << endl;

		 cout << "C_node: " << endl;
		 for (int tmpi = 0; tmpi < kServerCount; tmpi++) {
		 cout << C_node[tmpi] << " ";
		 }
		 cout << endl;

		 cout << "tau_F: " << endl;
		 for (int tmpi = 0; tmpi < kVNFCount; tmpi++) {
		 cout << tau_F[tmpi] << " ";
		 }
		 cout << endl;

		 cout << "tau_T: " << endl;
		 for (int tmpi = 0; tmpi < kTrafficCount; tmpi++) {
		 cout << tau_T[tmpi] << " ";
		 }
		 cout << endl;

		 cout << "tau_L: " << endl;
		 for (int tmpi = 0; tmpi < kLinkCount; tmpi++) {
		 cout << tau_L[tmpi] << " ";
		 }
		 cout << endl;

		 cout << "I_l_in: " << endl;
		 for (int tmpi = 0; tmpi < kLinkCount * 2; tmpi++) {
		 cout << "link " << tmpi << ":" << endl;
		 for (int tmpj = 0; tmpj < kServerCount; tmpj++) {
		 cout << I_l_in[tmpi][tmpj] << " ";
		 }
		 cout << endl;
		 }
		 cout << endl;

		 cout << "I_l_out: " << endl;
		 for (int tmpi = 0; tmpi < kLinkCount * 2; tmpi++) {
		 cout << "link " << tmpi << ":" << endl;
		 for (int tmpj = 0; tmpj < kServerCount; tmpj++) {
		 cout << I_l_out[tmpi][tmpj] << " ";
		 }
		 cout << endl;
		 }
		 cout << endl;

		 cout << "I_src_mc: " << endl;
		 for (int tmpj = 0; tmpj < kServerCount; tmpj++) {
		 cout << "node " << tmpj << ":" << endl;
		 for (int tmpi = 0; tmpi < kTrafficCount; tmpi++) {
		 cout << I_src_mc[tmpi][tmpj] << " ";
		 }
		 cout << endl;
		 }
		 cout << endl;

		 cout << "I_dst_mc: " << endl;
		 for (int tmpj = 0; tmpj < kServerCount; tmpj++) {
		 cout << "node " << tmpj << ":" << endl;
		 for (int tmpi = 0; tmpi < kTrafficCount; tmpi++) {
		 cout << I_dst_mc[tmpi][tmpj] << " ";
		 }
		 cout << endl;
		 }
		 cout << endl;
		 */

		cout << "Gonna solve the problem" << endl;
		IloTimer timer(env);
		timer.restart();
		cout << "set time limit" << endl;
		const IloInt timeLimit = 24 * 60 * 60; // one day
		const IloNum relativeGap = 0.0001; // find Integer solution within 0.1% of optimal
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::EpGap, relativeGap);
		//cplex.setParam(IloCplex::Threads, 8);
		cplex.setParam(IloCplex::MemoryEmphasis, true);
		cplex.setParam(IloCplex::PreDual, true);

		cplex.exportModel("lpex1.lp");

		cout << "solving the problem" << endl;
		if (!cplex.solve()) {
			timer.stop();
			cout << "!!!!!!!!!!!!!!!!!!! Solution Status = "
					<< cplex.getStatus() << " !!!!!!!!!!!!!!!!!!!" << endl;
			throw(-1);
		}
		timer.stop();
		opex = cplex.getObjValue();

		cout<<"opex: "<<opex<<endl;

		cout << "Solution Status = " << cplex.getStatus() << endl;
		cout << "Solution Value = " << opex << endl;
		for (int v = 0; v < kServerCount; v++) {
			utilization.push_back(0);
		}

		cout << "obtain provisioned nodes for each path" << endl;
		for (int r = 0; r < kTrafficCount; r++) {
			int scNum = _traffic_requests_SC[r].size();
			vector<int> ssfc = _traffic_requests[r].middlebox_sequence;
			//cout<<"r: "<<r<<endl;
			for(int cc = 0; cc < scNum; cc++){
				double cur_delay = 0;
				scSequence[r][cc].push_back(_traffic_requests[r].source);//.............
				vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
				int dim = tmpvec.size();
				for (int i = 1; i < dim; i++) {
					int fid = tmpvec[i];
					int vnfID = ssfc[fid];
					//cout<<vnfID<<" ";
					if(fid == ssfc.size()-1)
						continue;

					//int vnfID = tmpvec[i];
					//if(vnfID==middleboxes.size()-1)
					//	continue;
					for (int v = 0; v < kServerCount; v++) {
						if (fabs(cplex.getValue(alpha[r][fid][v]) - 1) < EPS) {
							scSequence[r][cc].push_back(v);
							utilization[v]++;
							cur_delay += middleboxes[vnfID].processing_delay;
						}
					}
				}
				//cout<<endl;
				if(tmpvec[dim-1] == ssfc.size()-1)
					scSequence[r][cc].push_back(_traffic_requests[r].destination);
				//this is the sum of transmission delay + vnf processing time, so the total delay should plus transport delay additionally
				delay_breakdown[r].push_back(cur_delay);
			}
		}

/*
		cout << "obtain provisioned nodes for each sub-chain" << endl;
		for (int c = 0; c < kSubchainCount; c++) {
			std::vector<int> tmpvec = _traffic_requests_SC[c].middlebox_sequence;
			if (tmpvec.size() > 0) {

				//find out the first provisioned node for SC
				if (tmpvec[0] == -1)
					scSequence[c].push_back(_traffic_requests_SC[c].source);
				else {
					int mcid = _traffic_requests_SC[c].mcID - mcOffset;
					int srcid = mcRank[mcid + mcOffset][tmpvec[0]];
					scSequence[c].push_back(mcSequence[mcid][srcid + 1]);
				}

				//find out the second provisioned ndoe for SC
				for (int v = 0; v < kServerCount; v++) {
					if (fabs(cplex.getValue(alpha_sc[c][v]) - 1) < EPS) {
						scSequence[c].push_back(v);
						utilization[v]++;
					}
				}
			}
		}
*/

		// retrieve provisioned edges for path
		cout << "obtain provisioned edges for each path" << endl;
		int edgesum = 0;
		for (int r = 0; r < kTrafficCount; r++) {
			unsigned long bw = _traffic_requests[r].min_bandwidth;
			int scNum = _traffic_requests_SC[r].size();
			vector<int> ssfc = _traffic_requests[r].middlebox_sequence;
			bool updatedVNF[50]={false};
			bool updatedEdges[50][50]={false};
			updatedVNF[0]=true;
			for(int cc = 0; cc < scNum; cc++){
				double cur_delay = delay_breakdown[r][cc];
				std::vector<std::pair<int, int> > links;
				vector<int> tmpvec = _traffic_requests_SC[r][cc].middlebox_sequence;
				vector<int> tmpres = scSequence[r][cc];
				for(int i = 0; i < tmpvec.size()-1; i++){
					int currvnf = tmpvec[i];
					int nextvnf = tmpvec[i+1];


					for (int l = 0; l < kLinkCount * 2; l++) {
						if (fabs(cplex.getValue(beta[r][currvnf][nextvnf][l]) - 1) < EPS) {
							links.push_back(link2nodePair[l]);
							//cur_delay += 1; // suppose the latency in each link is 1
							cur_delay += shortest_path[link2nodePair[l].first][link2nodePair[l].second];
							//if (r < 5)
							//	cout << "(" << link2nodePair[l].first << ", " << link2nodePair[l].second << ") ";
							//bandwidths[l]++;
							//bandwidths[l] += tmpbw;

							double resbw = 0.088;
							bool addbw = false;
							int sid=_traffic_requests[r].sid;
							if(sSFCres_file != "" && mpmap[sid].find( make_pair(currvnf, nextvnf)) != mpmap[sid].end() ){
								addbw = true;
							}

							unsigned long tmpbw = bw;
							//if the link is used to transmit packet header
							if(cc==0){
								bandwidths[r] += tmpbw;
							} else if(!updatedEdges[currvnf][nextvnf]){  //merge trunk
							//if(!updatedEdges[tmpres[i]][tmpres[i+1]]){   //merge link
								bandwidths[r] += tmpbw;
							} else if (addbw){
								//if the link is used to transmit packet header
								//bw = ...;
								//the average packet size in data centers is around 724 bytes.
								//For TCP packets, the header only occupies 8.8% of the total size.
								bandwidths[r] += tmpbw * resbw;
							}
						}
						edgesum += cplex.getValue(beta[r][currvnf][nextvnf][l]);
					}

					updatedEdges[currvnf][nextvnf] = true;
					//updatedEdges[tmpres[i]][tmpres[i+1]]=true;
					if(!updatedVNF[nextvnf]){
						updatedVNF[nextvnf]=true;
					}

				}
				delay_breakdown[r][cc] = cur_delay;
				//if (r < 5)
				//	cout << ";" << endl;
				scEdges[r][cc] = links;
			}
		}


		running_time = timer.getTime();

		cout << "path edge sum: " << edgesum << endl;
	} catch (IloException &e) {
		cerr << "Concert exception caught: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();
}

#endif  // MIDDLEBOX_PLACEMENT_SRC_CPLEX_H
