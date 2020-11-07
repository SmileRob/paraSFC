# Instructions


##compiling the source codes:

General command:
```shell
g++ -g -std=c++0x src/<algorithm>.cc -o <algorithm>
```
<algorithm> is mapper_parasfc, mapper_greed or mapper_coordvnf
 
Example:
```shell
g++ -g -std=c++0x src/mapper_parasfc.cc -o mapper_parasfc
```

##deploy sequential/parallelized SFC
General command:
```shell
./<algorithm> --topology_file=topo/<topo_name>.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-<topo_name> --parallel=true --outPath=res/<algorithm>.<topo_name>.paraSFC 
```




