#!/bin/bash


topolist="i2 india35 germany50"
alglist="parasfc"

g++ -g -std=c++0x -L/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio129/concert/lib/x86-64_linux/static_pic -I/opt/ibm/ILOG/CPLEX_Studio129/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio129/concert/include src/mapper_cplex.cc -lilocplex -lconcert -lcplex -lm -lpthread -ldl -DIL_STD -o mapper_cplex
	
g++ -g -std=c++0x src/mapper_parasfc.cc -o mapper_parasfc



###################### parallelized SFC #######################
if((0))
then
	for i in 100
	do
	    scal=$(echo "$i*0.015"|bc)
	    #echo $scal
		echo "=================== "$i" ==================="
		for alg in $alglist
		do
			topo="i2"
			#./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=true --max_time=$i --outPath=res/$alg.$topo.paraSFC --TopoScale=3 --bwexpan=0.3 #$scal
			
			topo="india35"
			./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=true --max_time=$i --outPath=res/$alg.$topo.paraSFC --TopoScale=12 --bwexpan=0.3  #$scal
			
			topo="germany50"
			#./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=true --max_time=$i --outPath=res/$alg.$topo.paraSFC --TopoScale=17  #$scal
		done
	done
fi



###################### sequential SFC #######################
if((1))
then
	for i in 100
	do
		echo "=================== "$i" ==================="
		for alg in $alglist
		do
			topo="i2"
			./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=false --max_time=$i --outPath=res/$alg.$topo.sequSFC --TopoScale=3 --bwexpan=0.5
			
			topo="india35"
			./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=false --max_time=$i --outPath=res/$alg.$topo.sequSFC --TopoScale=12 --bwexpan=0.5
			
			topo="germany50"
			./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=false --max_time=$i --outPath=res/$alg.$topo.sequSFC --TopoScale=17 --bwexpan=0.5
		done
	done
fi


###################### cplex #######################

if((0))
then
	g++ -g -std=c++0x -L/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio129/concert/lib/x86-64_linux/static_pic -I/opt/ibm/ILOG/CPLEX_Studio129/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio129/concert/include src/mapper_cplex.cc -lilocplex -lconcert -lcplex -lm -lpthread -ldl -DIL_STD -o mapper_cplex

	for i in 100
	do
	    scal=$(echo "$i*0.015"|bc)
	    #echo $scal
		echo "=================== "$i" ==================="
		for alg in $alglist
		do
			for topo in $topolist
			do
				./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=true --max_time=$i --outPath=res/$alg.$topo.paraSFC --TopoScale=$scal
				./mapper_$alg --topology_file=topo/$topo.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-$topo --parallel=false --max_time=$i --outPath=res/$alg.$topo.sequSFC --TopoScale=$scal
			done
		done
	done
fi

