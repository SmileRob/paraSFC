
## Prerequisites
- Ubuntu 16.04.7 or above.
- [CPLEX studio](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) version >= 12.09 based on 64-bit Linux platforms. Academics can obtain it via the [IBM Academic Initiative](https://developer.ibm.com/academic/).
- GCC 4.6 and above, that supports at least C++11.



## Installation of IBM ILOG CPLEX Optimization Studio on Linux 

Download IBM ILOG CPLEX Optimization Studio V12.9.0 or later versions from https://www.ibm.com/products/software .

The installer for IBM ILOG CPLEX Optimization Studio is distributed as a .bin file cplex_studioXXX.linux-x86.bin. This installer is used on 64-bit Linux platforms. .bin files are Linux self extracting files. Once you download the installer, follow the steps below:

Make sure the .bin file is executable. If necessary, change its permission using the chmod command from the directory where the .bin is located:
chmod +x cplex_studioXXX.linux-x86.bin

Enter the following command to start the installation process:
./cplex_studioXXX.linux-x86.bin


The XXX in the cplex_studioXXX.linux-x86.bin installer stands for the version number.

More information about the installation of cplex studio 12.9 can be found at [IBM Knowledge Center](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html)

Note: CPLEX Optimizer performed slightly better than [Google OR-tools](https://developers.google.com/optimization) on the classic benchmark, but was absolutely superior on the large-scale one. In fact, CP Optimizer was able to optimally solve 66% of the large-scale instances, against 29% of OR-Tools (quad core). By exploiting multi cores, OR-Tools was also able to find optimal solutions for the large-scale instances, showing that nowadays CP solvers in general could be successfully applied to real-world industrial problems.[1] 


# Instructions


## compiling the source codes:

General command:
```shell
g++ -g -std=c++0x src/<algorithm>.cc -o <algorithm>
```
<algorithm> can be one of the following options:
- mapper_parasfc
- mapper_greed
- mapper_coordvnf
 
Example:
```shell
g++ -g -std=c++0x src/mapper_parasfc.cc -o mapper_parasfc
```

## deploy sequential/parallelized SFC
General command:
```shell
./<algorithm> --topology_file=topo/<topo_name>.topo --middlebox_spec_file=middlebox-spec --traffic_request_file=dataset/request-traffic-<topo_name> --parallel=true --outPath=res/<algorithm>.<topo_name>.paraSFC 
```

[1]: Col, Giacomo Da and E. Teppan. “Google vs IBM: A Constraint Solving Challenge on the Job-Shop Scheduling Problem.” ArXiv abs/1909.08247 (2019): 259-265.

