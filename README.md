# paraSFC
This repository contains the codes of the paper, "On the Effective Parallelization and Near-Optimal Deployment of Service Function Chains", which is currently under review for IEEE Transactions on Parallel and Distributed Systems. The codes include several algorithms for Service Function Chains (SFCs) deployment: paraSFC, Greedy, CoordVNF, and the ILP-based algorithm.

## 1. Prerequisites
- Ubuntu 16.04.7, or above.
- GCC 5.4, or above.
- [IBM ILOG CPLEX Optimization Studio 12.9](https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/), for 64-bit Linux.


## 2. Install IBM ILOG CPLEX

Step 1: Download IBM ILOG CPLEX. Academics can obtain it via the [IBM Academic Initiative](https://developer.ibm.com/academic/). The installer is distributed as a .bin file *cplex_studio129.linux-x86.bin*.

Step 2: Make sure the .bin file is executable. If necessary, change its permission using the chmod command:
``` shell
chmod +x cplex_studio129.linux-x86.bin
```

Step 3: Enter the following command to start the installation process:
``` shell
./cplex_studio129.linux-x86.bin
```

More information about the installation of CPLEX can be found at [IBM Knowledge Center](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_installing.html).

Comparison to Google [OR-tools](https://developers.google.com/optimization): CPLEX performs slightly better on the classic benchmark, but is absolutely superior on large-scale problems. In fact, CPLEX is able to optimally solve 66% of the large-scale instances, against 29% of OR-tools (quad core). CPLEX in general can be successfully applied to real-world industrial problems. [1]


## 3. Compile and Run

```shell
cd parasfc-github
./run.sh
```

Please find the evaluation results in our paper.

## Acknowledgments

We would like to thank Bari et al. for sharing their implementation of NFV orchistration simulation [2].
The data structures, IO functions and util functions in this  are developed based on Bari's implementation.

## References
[1] Col, Giacomo Da and E. Teppan. “Google vs IBM: A Constraint Solving Challenge on the Job-Shop Scheduling Problem.” ArXiv abs/1909.08247 (2019): 259-265.

[2] F. Bari, S. R. Chowdhury, R. Ahmed, R. Boutaba, and O. C. M. B. Duarte, “Orchestrating virtualized network functions,” IEEE Transactions on Network and Service Management, vol. 13, no. 4, pp. 725–739, 2016.
