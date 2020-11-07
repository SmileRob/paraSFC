
## Requirements
-[Ubuntu 16.04.7 LTS]
-[gcc version >= 5.4.0]


## Installation of IBM ILOG CPLEX Optimization Studio on Linux 


The installer for IBM ILOG CPLEX Optimization Studio is distributed as a .bin file cplex_studioXXX.linux-x86.bin. This installer can be used on both 32-bit and 64-bit Linux platforms. .bin files are Linux self extracting files. The procedure is similar on a Mac OS platform. Once you download the installer, follow the steps below:

Make sure the .bin file is executable. If necessary, change its permission using the chmod command from the directory where the .bin is located:
chmod +x cplex_studioXXX.linux-x86.bin
Enter the following command to start the installation process:
./cplex_studioXXX.linux-x86.bin


Note: The XXX in the cplex_studioXXX.linux-x86.bin installer stands for the version number.


# Instructions


## compiling the source codes:

General command:
```shell
g++ -g -std=c++0x src/<algorithm>.cc -o <algorithm>
```
<algorithm> is mapper_parasfc, mapper_greed or mapper_coordvnf
 
Example:
```shell
g++ -g -std=c++0x src/mapper_parasfc.cc -o mapper_parasfc
```

## deploy sequential/parallelized SFC
General command:
```shell
./<algorithm> --topology_file=topo/<topo_name>.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-<topo_name> --parallel=true --outPath=res/<algorithm>.<topo_name>.paraSFC 
```




