#ifndef MIDDLEBOX_PLACEMENT_SRC_IO_H_
#define MIDDLEBOX_PLACEMENT_SRC_IO_H_

#include "datastructure_rob.h"
#include "util_rob.h"
#include <algorithm>
#include <string>
#include <string.h>
#include <iostream>
using namespace std;

std::unique_ptr<std::map<std::string, std::string> > ParseArgs(int argc,
		char *argv[]) {
	std::unique_ptr<std::map<std::string, std::string> > arg_map(
			new std::map<std::string, std::string>());
	for (int i = 1; i < argc; ++i) {
		char *key = strtok(argv[i], "=");
		char *value = strtok(NULL, "=");
		//DEBUG(" [%s] => [%s]\n", key, value);
		arg_map->insert(std::make_pair(key, value));
	}
	return std::move(arg_map);
}

inline int GetMiddleboxIndex(const std::string &middlebox_name) {
	//DEBUG("Finding middlebox: %s\n", middlebox_name.c_str());
	for (int i = 0; i < middleboxes.size(); ++i) {
		//DEBUG("i = %d, name = %s\n", i, middleboxes[i].middlebox_name.c_str());
		if (middleboxes[i].middlebox_name == middlebox_name) {
			//DEBUG("Middlebox %s found at %d\n", middlebox_name.c_str(), i);
			return i;
		}
	}
	//DEBUG("Middlebox %s not found :(\n", middlebox_name.c_str());
	return -1;  // Not found.
}

std::unique_ptr<std::vector<std::vector<std::string> > > ReadCSVFile_rob(
		const char *filename) {
	//DEBUG("[Parsing %s]\n", filename);
	FILE *file_ptr = fopen(filename, "r");
	const static int kBufferSize = 1024;
	char line_buffer[kBufferSize];

	std::unique_ptr<std::vector<std::vector<std::string> > > ret_vector(
			new std::vector<std::vector<std::string> >());
	std::vector<std::string> current_line;
	int row_number = 0;

	while (fgets(line_buffer, kBufferSize, file_ptr)) {
		current_line.clear();
		std::string strline = line_buffer;
		//cout<<"strline: ("<<strline<<")"<<endl;
		if (strline.length() > 1) {
			char *token = strtok(line_buffer, ",\n\r");
			current_line.push_back(token);
			while ((token = strtok(NULL, ",\n"))) {
				current_line.push_back(token);
			}
		}
		ret_vector->push_back(current_line);
	}

	fclose(file_ptr);
	//printf("Parsed %d lines\n", static_cast<int>(ret_vector->size()));
	return std::move(ret_vector);
}

std::unique_ptr<std::vector<std::vector<std::string> > > ReadCSVFile(
		const char *filename) {
	//DEBUG("[Parsing %s]\n", filename);
	FILE *file_ptr = fopen(filename, "r");
	const static int kBufferSize = 1024;
	char line_buffer[kBufferSize];
	std::unique_ptr<std::vector<std::vector<std::string> > > ret_vector(
			new std::vector<std::vector<std::string> >());
	std::vector<std::string> current_line;
	int row_number = 0;
	while (fgets(line_buffer, kBufferSize, file_ptr)) {
		current_line.clear();
		char *token = strtok(line_buffer, ",\n\r");
		current_line.push_back(token);
		while ((token = strtok(NULL, ",\n"))) {
			current_line.push_back(token);
		}
		ret_vector->push_back(current_line);
	}
	fclose(file_ptr);
	//DEBUG("Parsed %d lines\n", static_cast<int>(ret_vector->size()));
	return std::move(ret_vector);
}


void InitializeAllResults_CPLEX_oldversion(const char *filename) {

	string mcFile = filename;
	mcFile += "-MC";
	cout << "loading mc results: " << mcFile << endl;
	all_mcResults.clear();
	all_scResults.clear();

	auto csv_vector = ReadCSVFile_rob(mcFile.c_str());
	for (int i = 0; i < csv_vector->size(); ++i) {
		std::vector<std::string> &row = (*csv_vector)[i];
		std::vector<int> current_result;
		for (auto &element : row) {
			current_result.push_back(atoi(element.c_str()));
			//DEBUG("Pushing %s\n", element.c_str());
		}//DEBUG("\n");
		all_mcResults[i] = current_result;
		vector<vector<int> > tmpvec;
		all_scResults[i] = tmpvec;
	}

	if (traffic_requests_subchains.size() > 0) {
		string scFile = filename;
		scFile += "-SC";
		cout << "loading sc results: " << scFile << endl;
		auto SC_vector = ReadCSVFile_rob(scFile.c_str());
		for (int i = 0; i < SC_vector->size(); ++i) {
			std::vector<std::string> &row = (*SC_vector)[i];
			std::vector<int> current_result;
			for (int tmpi = 0; tmpi < row.size(); tmpi++) {
				current_result.push_back(atoi(row[tmpi].c_str()));
			}
			int mcid = scIDmcIDmap[i];
			all_scResults[mcid].push_back(current_result);
		}
	}

}

void InitializeAllResults(const char *filename) {

	string mcFile = filename;
	mcFile += "-MC";
	cout << "loading mc results: " << mcFile << endl;
	all_mcResults.clear();
	auto csv_vector = ReadCSVFile_rob(mcFile.c_str());
	for (int i = 0; i < csv_vector->size(); ++i) {
		std::vector<std::string> &row = (*csv_vector)[i];
		std::vector<int> current_result;
		for (auto &element : row) {
			current_result.push_back(atoi(element.c_str()));
		}
		all_mcResults[i] = current_result;
	}

	if (traffic_requests_subchains.size() > 0) {
		string scFile = filename;
		scFile += "-SC";
		cout << "loading sc results: " << scFile << endl;
		all_scResults.clear();
		auto SC_vector = ReadCSVFile_rob(scFile.c_str());
		for (int i = 0; i < SC_vector->size(); ++i) {
			std::vector<std::string> &row = (*SC_vector)[i];
			int mcid = atoi(row[0].c_str());
			std::vector<int> current_result;
			for (int tmpi = 1; tmpi < row.size(); tmpi++) {
				current_result.push_back(atoi(row[tmpi].c_str()));
			}
			all_scResults[mcid].push_back(current_result);
		}
	}

	cout << "finished loading all results" << endl;

	//print for debug
//	for (int tmpi = 0; tmpi < all_mcResults.size(); tmpi++) {
//		vector<int> vec = all_mcResults[tmpi];
//		for (int tmpjj = 0; tmpjj < vec.size(); tmpjj++) {
//			cout << vec[tmpjj] << ",";
//		}
//		cout << endl;
//	}
//
//	for (int tmpi = 0; tmpi < all_scResults.size(); tmpi++) {
//		vector<vector<int>> vecvec = all_scResults[tmpi];
//		for (int tmpjj = 0; tmpjj < vecvec.size(); tmpjj++) {
//			cout<<tmpi<<",";
//			for(int tmpkk = 0; tmpkk< vecvec[tmpjj].size(); tmpkk++)
//				cout << vecvec[tmpjj][tmpkk] << ",";
//			cout<<endl;
//		}
//		if(vecvec.size() == 0)
//			cout<<tmpi<<","<<endl;
//	}

}


void InitializeCplexPaths(const char *filename) {

	string mcFile = filename;
	mcFile += "-MC";
	cout << "loading mc results: " << mcFile << endl;

	auto csv_vector = ReadCSVFile_rob(mcFile.c_str());
	for (int i = 0; i < csv_vector->size(); ++i) {
		std::vector<std::string> &row = (*csv_vector)[i];
		std::vector<std::pair<int, int> > current_result;
		for (auto &element : row) {
			string tmpstr = element;
			if(tmpstr.length()<5)
				continue;
			int pos = tmpstr.find("-",1);
			string n1 = tmpstr.substr(1,pos-1);
			string n2 = tmpstr.substr(pos+1,tmpstr.length()-pos-2);
			current_result.push_back(make_pair(atoi(n1.c_str()),atoi(n2.c_str())));
			//cout<<atoi(n1.c_str())<<" "<<atoi(n2.c_str())<<endl;
		}
		all_mcPaths[i] = current_result;

		vector<vector<std::pair<int, int> > > tmpvar;
		all_scPaths[i] = tmpvar;
	}

	if (traffic_requests_subchains.size() > 0) {
		string scFile = filename;
		scFile += "-SC";
		cout << "loading sc results: " << scFile << endl;
		auto SC_vector = ReadCSVFile_rob(scFile.c_str());
		//cout<<"SC_vector size: "<<SC_vector->size()<<endl;
		for (int i = 0; i < SC_vector->size(); i++) {
			std::vector<std::string> &row = (*SC_vector)[i];
			std::vector<std::pair<int, int> > current_result;
			if(row.size() != 0){
				for (auto &element : row) {
					string tmpstr = element;
					if (tmpstr.length() < 5)
						continue;
					int pos = tmpstr.find("-", 1);
					string n1 = tmpstr.substr(1, pos - 1);
					string n2 = tmpstr.substr(pos + 1, tmpstr.length() - pos - 2);
					current_result.push_back(make_pair(atoi(n1.c_str()), atoi(n2.c_str())));
					//if (i < 10) cout << atoi(n1.c_str()) << " -- " << atoi(n2.c_str()) << endl;
				}
			}

			int mcid = scIDmcIDmap[i];
			all_scPaths[mcid].push_back(current_result);
			//if(i<10)cout<<"mcid: "<<mcid<<" scid: "<<i<<" "<<row.size()<<" "<<all_scPaths[mcid].size()<<endl;
		}
	}

	//cout << "finished loading all results" << endl;

}


void InitializeSolutionPaths(const char *filename) {
	paths.clear();
	auto csv_vector = ReadCSVFile(filename);
	for (int i = 0; i < csv_vector->size(); ++i) {
		std::vector<std::string> &row = (*csv_vector)[i];
		std::vector<int> current_path;
		for (auto &element : row) {
			current_path.push_back(atoi(element.c_str()));
		}
		paths.push_back(current_path);
	}
}

void InitializeMiddleboxes(const char *filename) {
	middleboxes.clear();
	middleboxes.emplace_back("in","0","0","0","0");
	auto csv_vector = ReadCSVFile_rob(filename);
	for (int i = 0; i < csv_vector->size(); ++i) {
		std::vector<std::string> &row = (*csv_vector)[i];
		middleboxes.emplace_back(row[0], row[1], row[2], row[3], row[4]);
		//cout<<row[0]<<" "<<row[1]<<" "<< row[2]<<" "<<row[3]<<" "<<row[4]<<endl;
	}
	middleboxes.emplace_back("e","0","0","0","0");
}

void PrintMiddleboxes() {
	printf("[Middleboxes (count = %d)]\n",
			static_cast<int>(middleboxes.size()));
	for (int i = 0; i < middleboxes.size(); ++i) {
		printf("[i = %d] %s\n", i, middleboxes[i].GetDebugString().c_str());
	}
}

void PrintTrafficRequests() {
	printf("[Traffic Requests (count = %d)\n",
			static_cast<int>(traffic_requests.size()));
	for (int i = 0; i < traffic_requests.size(); ++i) {
		printf("[i = %d] %s\n", i,
				traffic_requests[i].GetDebugString().c_str());
	}
}

void InitializeTrafficRequestsSubChains_rob(const char *filename) {

	string reqfile = filename;
	reqfile += "-SCs";
	//cout<<"filename: "<<reqfile<<endl;

	int time_slot = 0, pre_time = -1;

	int sfcOffset = 7;  // the offset to locate the SFC
	traffic_requests_subchains.clear();
	auto csv_vector = ReadCSVFile(reqfile.c_str());

	//cout<<"for each line the file"<<endl;
	for (int i = 0; i < csv_vector->size(); ++i) {
		std::vector<int> mbox_sequence;
		std::vector<std::string> &row = (*csv_vector)[i];

		//the second column is time slot
		int cur_time = atoi(row[1].c_str());
		if (cur_time != pre_time) {
			pre_time = cur_time;
			time_slot++;
			if (time_slot > max_time)
				break;
		}

		//cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<" "<<row[3]<<" "<<row[4]<<" "<<row[5]<<" "<<row[6]<<": ";
		for (int mbox_seq_index = sfcOffset; mbox_seq_index < row.size();
				++mbox_seq_index) {
			mbox_sequence.push_back(atoi(row[mbox_seq_index].c_str()));
			if(atoi(row[mbox_seq_index].c_str()) == 0)
				break;
			//mbox_sequence.push_back(GetMiddleboxIndex(row[mbox_seq_index]));
			//cout<<GetMiddleboxIndex(row[mbox_seq_index])<<" ("<<row[mbox_seq_index]<<") ";
		}
		//cout<<endl;

		int mcid = atoi(row[0].c_str());

		//traffic_requests_subchains[mcid].emplace_back(row[1], row[2], row[3], row[4], row[5], row[6], mbox_sequence);
		//traffic_requests_subchains[mcid].emplace_back(row[0], row[1], row[2], "1", mbox_sequence);


	}


	//cout<<"Finished InitializeTrafficRequestsSubChains_rob"<<endl;
}

//Tree
void InitializeSFCMap_tree(){

	//cout<<"Initialize the SFC maps"<<endl;
	vector<int> webs = {0,6,5,1,2,11};
	vector<int> voip = {0,6,9,3,1,10,6,11};
	vector<int> vdio = {0,1,2,7,8,3,4,6,11};
	vector<int> olgm = {0,6,3,5,6,11};
	sSFC = {webs,voip,vdio,olgm};

	vector<vector<int> > webpsfc = {{0,1,2,5},{0,1,5},{0,1,3},{0,1,4}};
	//vector<vector<int> > webpsfc = {{0,6,5,11},{0,6,11},{0,6,1},{0,6,2}};
	pSFC[0] = webpsfc;

	vector<vector<int> > voipsfc = {{0,1,2,5,6,7},{0,1,2,3,7},{0,1,2,4}};
	//vector<vector<int> > voipsfc = {{0,6,9,10,6,11},{0,6,9,3,11},{0,6,9,1}};
	pSFC[1] = voipsfc;

	vector<vector<int> > vdiosfc = {{0,3,4,5,7,8},{0,3,4,6,7,8},{0,3,4,7,8},{0,1},{0,2}};
	//vector<vector<int> > vdiosfc = {{0,7,8,3,6,11},{0,7,8,4,6,11},{0,7,8,6,11},{0,1},{0,2}};
	pSFC[2] = vdiosfc;

	vector<vector<int> > olgmsfc = {{0,1,3,4,5},{0,1,2,4,5},{0,1,4,5}};
	//vector<vector<int> > olgmsfc = {{0,6,5,6,11},{0,6,3,6,11},{0,6,6,11}};
	pSFC[3] = olgmsfc;

}


// DAG
void InitializeSFCMap(){

	//cout<<"Initialize the SFC maps"<<endl;
	vector<int> webs = {0,6,5,1,2,11};
	vector<int> voip = {0,6,9,3,1,10,6,11};
	vector<int> vdio = {0,1,2,7,8,3,4,6,11};
	vector<int> olgm = {0,6,3,5,6,11};
	sSFC = {webs,voip,vdio,olgm};

	vector<vector<int> > webpsfc = {{0,1,2,5},{0,1,5},{0,1,3,5},{0,1,4,5}};
	//vector<vector<int> > webpsfc = {{0,6,5,11},{0,6,11},{0,6,1},{0,6,2}};
	pSFC[0] = webpsfc;

	vector<vector<int> > voipsfc = {{0,1,2,5,6,7},{0,1,2,3,7},{0,1,2,4,7}};
	//vector<vector<int> > voipsfc = {{0,6,9,10,6,11},{0,6,9,3,11},{0,6,9,1}};
	pSFC[1] = voipsfc;

	vector<vector<int> > vdiosfc = {{0,3,4,5,7,8},{0,3,4,6,7,8},{0,3,4,7,8},{0,1,8},{0,2,8}};
	//vector<vector<int> > vdiosfc = {{0,7,8,3,6,11},{0,7,8,4,6,11},{0,7,8,6,11},{0,1},{0,2}};
	pSFC[2] = vdiosfc;

	vector<vector<int> > olgmsfc = {{0,1,3,4,5},{0,1,2,4,5},{0,1,4,5}};
	//vector<vector<int> > olgmsfc = {{0,6,5,6,11},{0,6,3,6,11},{0,6,6,11}};
	pSFC[3] = olgmsfc;



	map<pair<int, int>, int> webcpmg={ {{3,5},1},  {{2,5},1},  {{4,5}, 1} };
	mpmap[0] = webcpmg;


	map<pair<int, int>, int> voipcpmg={ {{4,7},1},  {{3,7},1},  {{6,7}, 1} };
	mpmap[1] = voipcpmg;

	map<pair<int, int>, int> vdiocpmg={ {{1,8},1},  {{2,8},1},  {{7,8}, 1}, {{4,7},1},  {{5,7},1},  {{6,7}, 1} };
	mpmap[2] = vdiocpmg;

	map<pair<int, int>, int> olgmcpmg={ {{1,4},1},  {{3,4},1},  {{2,4}, 1} };
	mpmap[3] = olgmcpmg;
}


void InitializeSeqSFCMap(){

	//cout<<"Initialize the SFC maps"<<endl;
	vector<int> webs = {0,6,5,1,2,11};
	vector<int> voip = {0,6,9,3,1,10,6,11};
	vector<int> vdio = {0,1,2,7,8,3,4,6,11};
	vector<int> olgm = {0,6,3,5,6,11};
	sSFC = {webs,voip,vdio,olgm};

	vector<vector<int> > webpsfc = {{0,1,2,3,4,5}};
	//vector<vector<int> > webpsfc = {{0,6,5,1,2,11}};
	pSFC[0] = webpsfc;

	vector<vector<int> > voipsfc = {{0,1,2,3,4,5,6,7}};
	//vector<vector<int> > voipsfc = {{0,6,9,3,1,10,6,11}};
	pSFC[1] = voipsfc;

	vector<vector<int> > vdiosfc = {{0,1,2,3,4,5,6,7,8}};
	//vector<vector<int> > vdiosfc = {{0,1,2,7,8,3,4,6,11}};
	pSFC[2] = vdiosfc;

	vector<vector<int> > olgmsfc = {{0,1,2,3,4,5}};
	//vector<vector<int> > olgmsfc = {{0,6,3,5,6,11}};
	pSFC[3] = olgmsfc;

}

void InitializeSFCMap_cplex(){

	//cout<<"Initialize the SFC maps"<<endl;
	vector<int> webs = {0,6,5,1,2,11};
	vector<int> voip = {0,6,9,3,1,10,6,11};
	vector<int> vdio = {0,1,2,7,8,3,4,6,11};
	vector<int> olgm = {0,6,3,5,6,11};
	sSFC = {webs,voip,vdio,olgm};

	vector<vector<int> > webpsfc = {{0,1,2,5},{0,1,5},{0,1,3,5},{0,1,4,5}};
	pSFC[0] = webpsfc;

	vector<vector<int> > voipsfc = {{0,1,2,5,6,7},{0,1,2,3,7},{0,1,2,4,7}};
	pSFC[1] = voipsfc;

	vector<vector<int> > vdiosfc = {{0,3,4,5,7,8},{0,3,4,6,7,8},{0,3,4,7,8},{0,1,8},{0,2,8}};
	pSFC[2] = vdiosfc;

	vector<vector<int> > olgmsfc = {{0,1,3,4,5},{0,1,2,4,5},{0,1,4,5}};
	pSFC[3] = olgmsfc;


}


void InitializeSeqSFCMap_cplex(){

	//cout<<"Initialize the SFC maps"<<endl;
	vector<int> webs = {0,6,5,1,2,11};
	vector<int> voip = {0,6,9,3,1,10,6,11};
	vector<int> vdio = {0,1,2,7,8,3,4,6,11};
	vector<int> olgm = {0,6,3,5,6,11};
	sSFC = {webs,voip,vdio,olgm};

	vector<vector<int> > webpsfc = {{0,1,2,3,4,5}};
	pSFC[0] = webpsfc;

	vector<vector<int> > voipsfc = {{0,1,2,3,4,5,6,7}};
	pSFC[1] = voipsfc;

	vector<vector<int> > vdiosfc = {{0,1,2,3,4,5,6,7,8}};
	pSFC[2] = vdiosfc;

	vector<vector<int> > olgmsfc = {{0,1,2,3,4,5}};
	pSFC[3] = olgmsfc;

}


#if 0

void pathlister(int curnode, vector <int> path){

	path.push_back(curnode);
	if((vnfnext.find(curnode)) == vnfnext.end()){
		pathmap[pathmap.size()] = path;
		//for(int i=0; i < path.size(); i++)
		//	cout<<path[i]<<" ";
		//cout<<endl;
		return;
	}
	set<int> nextnodes = vnfnext[curnode];
	for(auto nextnode: nextnodes){
		pathlister(nextnode, path);
	}
}



void pathfinder(vector<int> paravec, vector<int> mbox_sequence){


	int prenod = 0, curnod = 0;
	int seqlen = mbox_sequence.size();
	for(int tmpi = 0; tmpi < seqlen; tmpi++){

		//cout<<"***tmpi: "<<tmpi<<endl;//<<" "<<mbox_sequence[tmpi]<<endl;
		curnod++;


		if (paravec[tmpi] == 0) {

			//cout<<"first node: "<<prenod<<" "<<curnod<<endl;
			vnfnext[prenod].insert(curnod);
			vnfpre[curnod].insert(prenod);

		} else if (paravec[tmpi] == 1) {

			//cout<<"single node: "<<prenod<<" "<<curnod<<endl;
			vnfnext[prenod].insert(curnod);
			vnfpre[curnod].insert(prenod);


			if((vnfpre.find(prenod)) != vnfpre.end()){
				set<int> preprenodes = vnfpre[prenod];
				for(auto preprenode: preprenodes){
					if((vnfnext.find(preprenode)) != vnfnext.end()){
						set<int> prenextnodes = vnfnext[preprenode];
						for(auto prenextnode: prenextnodes){
							vnfnext[prenextnode].insert(curnod);
							vnfpre[curnod].insert(prenextnode);
						}
					}
				}
			}

		}
		else if (paravec[tmpi] == 2){
			//cout<<"para node: "<<prenod<<" "<<curnod<<endl;

			if((vnfpre.find(prenod)) != vnfpre.end()){
				set<int> preprenodes = vnfpre[prenod];
				for(auto preprenode: preprenodes){
					vnfnext[preprenode].insert(curnod);
					vnfpre[curnod].insert(preprenode);
				}
			}

		}
		prenod = curnod;
	}


	//cout<<"call pathlister"<<endl;

	pathmap.clear();
	curnod = 0;
	vector<int> path;
	pathlister(curnod, path);
	vnfnext.clear();
	vnfpre.clear();
}

#endif

void InitializeTrafficRequests_rob(const char *filename, unsigned long testbw) {

	string reqfile = filename;
	//cout << "filename: " << reqfile << endl;

	int totalrequirement = 0;
	int servicerequirement[4] = {22,28,38,19};


	string parafile = filename;
	ifstream infile(parafile.c_str());
	auto mat_vector = ReadCSVFile(parafile.c_str());
	int matpos = 0;
	traffic_requests.clear();
	int time_slot = -1, pre_time = -1;
	auto csv_vector = ReadCSVFile(reqfile.c_str());
	//cout<<csv_vector->size()<<" "<<max_time<<endl;

	int mcid = 0;
	for (mcid = 0; mcid < csv_vector->size() && time_slot < max_time; mcid++) {  //!!!!!!!!!!!!!!!!!!!!control time slot
	//for (mcid = 0; mcid < csv_vector->size() && mcid < max_time; mcid++) {	//!!!!!!!!!!!!!!!!!!!!!control mcid
		std::vector<std::string> &row = (*csv_vector)[mcid];

		//the first column is time slot
		int cur_time = atoi(row[0].c_str());

		//cout<<cur_time<<" "<<pre_time<<endl;
		if (cur_time != pre_time) {
			pre_time = cur_time;
			time_slot++;
			//cout<<"timeslot and mcid: "<<time_slot<<" "<<mcid<<endl;
			if (time_slot >= max_time) break;     //!!!!!!!!!!!!!!!!!!!!control time slot
		}

		vector<int> tmpvec;
		mcIDscIDmap[mcid] = tmpvec;

		int sid = atoi(row[6].c_str());
		std::vector<int> mbox_sequence = sSFC[sid];
		vector<vector<int> > paths = pSFC[sid];

		totalrequirement += servicerequirement[sid];


		//traffic_requests.emplace_back(row[0], row[1], row[2], row[3], row[4], row[5], mbox_sequence);
		traffic_requests.emplace_back(row[0], row[1], row[2], row[3], row[4], row[5], sid, mbox_sequence, paths);
		traffic_requests[mcid].mcID = mcid;
		traffic_requests[mcid].originalID = mcid;
		reqSFCLen.push_back(mbox_sequence.size());
		//cout<<mcid<<": "<<traffic_requests[mcid].min_bandwidth<<" "<<traffic_requests[mcid].sid<<endl;
		if(testbw != -1){
			traffic_requests[mcid].min_bandwidth = testbw;
		}

		//cout<<"to load sub chains"<<endl;

#if 1
		//if(mcid < 10)cout<<mcid<<": "<<traffic_requests[mcid].middlebox_sequence.size()<<" ";
		for(int scid = 0; scid < paths.size(); scid++)
		{
			vector<int> subseq = paths[scid];
			//if(mcid < 10)cout<<subseq.size()<<" ";
			traffic_requests_subchains[mcid].emplace_back(row[0], row[1], row[2], row[3], row[4], row[5], sid, subseq, paths);
			//cout<<traffic_requests_subchains[mcid][scid].middlebox_sequence.size()<<" ";
			if(testbw != -1){
				traffic_requests_subchains[mcid][scid].min_bandwidth = testbw;
			}
		}
		//if(mcid < 10)cout<<endl;
#endif



/*
		if(mcid == 11)
		{
			vector<int> seqq = traffic_requests[mcid].middlebox_sequence;
			for(int iii = 0; iii < seqq.size(); iii++)
			{
				cout<<seqq[iii]<<" ";
			}
			cout<<endl;

			vector<vector<int> > ps = traffic_requests[mcid].psfcpaths;
			for(int iii = 0; iii < ps.size(); iii++){
				vector<int> ppp = ps[iii];
				for(int jjj = 0; jjj< ppp.size(); jjj++){
					cout<<ppp[jjj]<<" ";
				}
				cout<<endl;
			}
		}
		*/


	}

	//cout<<"to calculate the duration for each request"<<endl;
	int last_time_stamp = max_time;
	int current_time = traffic_requests.back().arrival_time;
	for (int i = traffic_requests.size() - 1; i >= 0; i--) {

		while (current_time > max_time)  // logic correction by Rob
		{
			current_time = traffic_requests[i].arrival_time;
			i--;
		}

		if (current_time != traffic_requests[i].arrival_time) {
			last_time_stamp = current_time;
			current_time = traffic_requests[i].arrival_time;
		}
		traffic_requests[i].duration = (last_time_stamp - current_time) * 60;
		traffic_requests[i].accepted = true;
	}

	cout<<"Total resource requirement: "<<totalrequirement<<" ";
	cout<<"traffic_requests size: "<<traffic_requests.size()<<endl;

}




void InitializeTopology(const char *filename, double bwscale) {

	//printf("[Parsing %s]\n", filename);
	FILE *file_ptr = fopen(filename, "r");
	int node_count, edge_count;
	int ret_val = fscanf(file_ptr, "%d %d", &node_count, &edge_count);
	//printf(" node_count = %d, edge_count = %d\n", node_count, edge_count);
	graph.resize(node_count);
	nodes.resize(node_count);
	deployed_mboxes.resize(node_count);

	//cout<<"fetch nodes, including node ID and node capacity"<<endl;
	for (int i = 0; i < node_count; ++i) {
		int _nodeid, _nodecore;
		ret_val = fscanf(file_ptr, "%d %d", &_nodeid, &_nodecore);
		//cout<<"nodeid: "<<_nodeid<<" _nodecore: "<<_nodecore<<endl;
		//ret_val = fscanf(file_ptr, "%d %d", &nodes[i].node_id, &nodes[i].num_cores);
		nodes[i].node_id = _nodeid;
		nodes[i].num_cores = _nodecore*topoScale;//!!!!!!!!!!!!!!!!
		nodes[i].residual_cores = nodes[i].num_cores;
		graph[i].clear();
		shortest_path[i][i] = 0.0;
		shortest_edge_path[i][i] = 0;
		bw[i][i] = 0;
		sp_pre[i][i] = NIL;
		for (int j = i + 1; j < node_count; j++) {
			shortest_path[i][j] = shortest_path[j][i] = INF;
			shortest_edge_path[i][j] = shortest_edge_path[j][i] = INF;
			sp_pre[i][j] = sp_pre[j][i] = NIL;
			bw[i][j] = bw[j][i] = 0;
		}
	}
	//printf("nodes.size() = %u\n", nodes.size());
	for (auto &n : nodes) {
		DEBUG("%s\n", n.GetDebugString().c_str());
	}

	//cout<<"fetch edges, including bandwidth and delay of edges"<<endl;
	for (int j = 0; j < edge_count; ++j) {
		int source, destination, delay;
		unsigned long bandwidth;
		double latency;
		ret_val = fscanf(file_ptr, "%d %d %lu %lf", &source, &destination, &bandwidth, &latency);
		//source = source-1;
		//destination = destination -1;
		//bandwidth = bandwidth*bwscale;
		bandwidth = bwscale;
		//cout<<"bw and scale: "<<bandwidth<<" "<<bwscale<<endl;


		delay = latency;
		//cout<<source<<" "<<destination<<" "<<bandwidth<<" "<<latency<<" "<<delay<<endl;
		graph[source].emplace_back(&nodes[destination], bandwidth, delay);
		graph[destination].emplace_back(&nodes[source], bandwidth, delay);


		bw[source][destination] = bw[destination][source] = bandwidth;
		shortest_edge_path[source][destination] = 1;
		shortest_edge_path[destination][source] = 1;
		shortest_path[source][destination] = shortest_path[destination][source] = delay;
		sp_pre[source][destination] = source;
		sp_pre[destination][source] = destination;
	}

	// ??????, does it compute the shortest path as well as delay? If so, does it make sense?
	//cout<<"here should use the dij. alg. to calculate the result!!!!!!"<<endl;
	bool again = true;
	int maxiter = 10*node_count*node_count;
	int tmpiter = 0;
	while(again)
	{
		again = false;
		for (int k = 0; k < node_count; ++k) {
			for (int i = 0; i < node_count; ++i) {
				for (int j = 0; j < node_count; ++j) {
					if (i == j)
						continue;

					//finding shortest path measured by hop
//					int relaxed_cost = shortest_edge_path[i][k] + shortest_edge_path[k][j];
//					if (shortest_edge_path[i][j] > relaxed_cost) {
//						shortest_edge_path[i][j] = relaxed_cost;
//						shortest_path[i][j] = shortest_path[i][k] + shortest_path[k][j];
//						sp_pre[i][j] = sp_pre[k][j];
//						if(sp_pre[i][j] == NIL)
//							again = true;
//					}

					//finding shortest path measured by delay
					int relaxed_cost = shortest_path[i][k] + shortest_path[k][j];
					if (shortest_path[i][j] > relaxed_cost) {
						shortest_path[i][j] = relaxed_cost;
						shortest_edge_path[i][j] = shortest_edge_path[i][k] + shortest_edge_path[k][j];
						sp_pre[i][j] = sp_pre[k][j];
						if(sp_pre[i][j] == NIL || shortest_path[i][j] >= INF)
							again = true;
					}
				}
			}
		}

		//if(tmpiter++>maxiter)
		//	again = false;
	}

	//cout<<"compute the closeness"<<endl;
	closeness.resize(node_count);
	for (int i = 0; i < node_count; ++i) {
		double farness = 0.0;
		for (int j = 0; j < node_count; ++j) {
			if (i != j) {
				farness += shortest_edge_path[i][j];
			}
		}
		closeness[i] = 1.0 / farness;
	}
	fclose(file_ptr);


//	for(int i = 0; i < node_count; i++){
//		for(int j = 0; j < node_count; j++)
//			cout<<bw[i][j]<<" ";
//		cout<<endl;
//	}
//	cout<<"****************"<<endl;
//	for(int i = 0; i < node_count; i++){
//		for(int j = 0; j < node_count; j++)
//			cout<<shortest_path[i][j]<<" ";
//		cout<<endl;
//	}

}

#endif  // MIDDLEBOX_PLACEMENT_SRC_IO_H_
