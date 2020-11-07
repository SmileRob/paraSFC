
## Prerequisites
- Ubuntu 16.04.7 or above.
- [CPLEX studio](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) version >= 12.09 based on 64-bit Linux platforms. Academics can obtain it via the [IBM Academic Initiative](https://developer.ibm.com/academic/).
- GCC 4.6 and above, that supports at least C++11.



## Installation of IBM ILOG CPLEX Optimization Studio on Linux 

Download IBM ILOG CPLEX Optimization Studio V12.9.0 or later versions from https://www.ibm.com/products/software .

IBM ILOG CPLEX Optimization Studio provides the most efficient way of building models for mathematical programming, constraint programming and constraint-based scheduling, in order to tackle complex optimization problems such as planning and scheduling.

The installer for IBM ILOG CPLEX Optimization Studio is distributed as a .bin file cplex_studioXXX.linux-x86.bin. This installer is used on 64-bit Linux platforms. .bin files are Linux self extracting files. Once you download the installer, follow the steps below:

Make sure the .bin file is executable. If necessary, change its permission using the chmod command from the directory where the .bin is located:
chmod +x cplex_studioXXX.linux-x86.bin

Enter the following command to start the installation process:
./cplex_studioXXX.linux-x86.bin


The XXX in the cplex_studioXXX.linux-x86.bin installer stands for the version number.

More information about the installation of cplex studio 12.9 can be found at https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html


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




