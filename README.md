# Hard-Shape-Constraints-for-Kernels
The experiments in the NeurIPS 2020 article "Hard Shape-Constrained Kernel Machines", Pierre-Cyril Aubin-Frankowski and Zoltan Szabo, https://arxiv.org/abs/2005.12636 can be reproduced with the Matlab code provided here.

*Shape constraints (such as non-negativity, monotonicity, convexity) play a central role in a large number of applications, as they usually improve performance for small sample size and help interpretability. However enforcing these shape requirements in a hard fashion is an extremely challenging problem. Classically, this task is tackled (i) in a soft way (without out-of-sample guarantees), (ii) by specialized transformation of the variables on a case-by-case basis, or (iii) by using highly restricted function classes, such as polynomials or polynomial splines. In this paper, we prove that hard affine shape constraints on function derivatives can be encoded in kernel machines which represent one of the most flexible and powerful tools in machine learning and statistics. Particularly, we present a tightened second-order cone constrained reformulation, that can be readily implemented in convex solvers. We prove performance guarantees on the solution, and demonstrate the efficiency of the approach in joint quantile regression with applications to economics and to the analysis of aircraft trajectories, among others.*

The only dependency of the code is to CVXGEN to compute the solution of the second-order-cone (SOC) constrained optimization problems (i.e. ridge regression and joint quantile regression). Please download the necessary files at http://cvxr.com/cvx/download/ The download will generate a 50Mo file, no other library manipulation should be required. We used Matlab 2018a. We recommend using the Mosek solver if available for an average two-fold speedup with respect to the default CVX solver (SeDuMi or SDPT3).

There is no Python version available. Nevertheless the code is straightforward to implement, especially when using a modeling language for convex optimization problems, such as CVXPY https://www.cvxpy.org/. 

The .m files correspond to the corresponding figures in the article:

ToyQuadratic_MonotonousLeastSquares.m->Fig.2 and 3

Engel_MonotonousJointQuantReg.m->Fig.4a

Engel_ConcaveMonotonousJointQuantReg.m->Fig.4b

ENAC_MonotonousJointQuantReg.m->Fig.5

Each of these files was re-ran before submission and allows to reproduce faithfully the results of the article. However, for the testers' convenience, some datasets were cropped to keep short computation times (less than a minute on a laptop, an indicative value is provided in each .m file). If you wish to do the full experiment, uncomment the code lines as suggested in the respective .m files.

Pierre-Cyril Aubin 02/10/2020
