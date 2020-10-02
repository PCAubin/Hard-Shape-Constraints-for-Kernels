# Hard-Shape-Constraint-for-Kernels
The experiments in the NeurIPS 2020 article "Hard Shape-Constrained Kernel Machines", Pierre-Cyril Aubin-Frankowski and Zoltan Szabo, https://arxiv.org/abs/2005.12636 can be reproduced with the Matlab code provided here.

The only dependency is to CVXGEN to compute the solution of the second-order-cone (SOC) constrained optimization problems (i.e. ridge regression and joint quantile regression). Please download the necessary files at http://cvxr.com/cvx/download/ The download will generate a 50Mo file, no other library manipulation should be required. We used Matlab 2018a. We recommend using the Mosek solver is available for an average two-fold speedup.

There is no Python version available. Nevertheless the code is straightforward to implement, especially when using a modeling language for convex optimization problems, such as https://www.cvxpy.org/. 

The .m files correspond to these figures in the article:

ToyQuadratic_MonotonousLeastSquares.m->Fig.2 and 3
Engel_MonotonousJointQuantReg.m->Fig.4a
Engel_ConcaveMonotonousJointQuantReg.m->Fig.4b
ENAC_MonotonousJointQuantReg.m->Fig.5

Each of these files was re-ran before submission and allows to reproduce faithfully the results of the article. However, for the tester' convenience, some datasets were cropped to keep short computation times (less than a minute on a laptop, an indicative value is provided in each .m file). If you wish to do the full experiment, uncomment the code lines as suggested in the respective .m files.

Pierre-Cyril Aubin 02/10/2020
