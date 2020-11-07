# paraSFC
This repository contains codes of several algorithms for SFC deployment, including paraSFC, Greedy, CoordVNF and ILP (Integer Linear Programming)-based algorithm.

## Prerequisites
- Ubuntu 16.04.7 or above.
- GCC 4.6 and above, that supports at least C++11.
- [CPLEX studio](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) (version >= 12.09) based on 64-bit Linux platforms. Academics can obtain it via the [IBM Academic Initiative](https://developer.ibm.com/academic/).

### Installation of IBM ILOG CPLEX Optimization Studio on Linux 

Download IBM ILOG CPLEX Optimization Studio V12.9.0 or later versions from https://www.ibm.com/products/software.

The installer for IBM ILOG CPLEX Optimization Studio is distributed as a .bin file cplex_studioXXX.linux-x86.bin. This installer is used on 64-bit Linux platforms. .bin files are Linux self extracting files. Once you download the installer, follow the steps below:

Make sure the .bin file is executable. If necessary, change its permission using the chmod command from the directory where the .bin is located:
chmod +x cplex_studioXXX.linux-x86.bin

Enter the following command to start the installation process:
./cplex_studioXXX.linux-x86.bin


The XXX in the cplex_studioXXX.linux-x86.bin installer stands for the version number.

More information about the installation of cplex studio 12.9 can be found at [IBM Knowledge Center](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html)

Note: CPLEX Optimizer performed slightly better than [OR-tools](https://developers.google.com/optimization) on the classic benchmark, but was absolutely superior on the large-scale one. In fact, CPLEX Optimizer was able to optimally solve 66% of the large-scale instances, against 29% of OR-Tools (quad core). By exploiting multi cores, OR-Tools was also able to find optimal solutions for the large-scale instances, showing that nowadays CPLEX solvers in general could be successfully applied to real-world industrial problems.[1] 



## compile and run

General command:
```shell
cd parasfc-github
./run.sh
```

[1]: Col, Giacomo Da and E. Teppan. “Google vs IBM: A Constraint Solving Challenge on the Job-Shop Scheduling Problem.” ArXiv abs/1909.08247 (2019): 259-265.

