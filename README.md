# Being Robust (in High Dimensions) Can Be Practical
A MATLAB implementation of [Being Robust (in High Dimensions) Can Be Practical](https://arxiv.org/abs/1703.00893) from ICML 2017.

Prerequisites 
===
This project requires installation of the following packages:
-[AgnosticMeanAndCovarianceCode](https://github.com/kal2000/AgnosticMeanAndCovarianceCode): An implementation of the algorithms from the paper [Agnostic Estimation of Mean and Covariance](https://arxiv.org/abs/1604.06968) from FOCS 2016, authored by [Kevin A. Lai](https://www.cc.gatech.edu/~klai9/), [Anup B. Rao](https://sites.google.com/site/anupraob/), and [Santosh Vempala](https://www.cc.gatech.edu/~vempala/). Used to compare the statistical accuracy of our algorithms against theirs.
-[CVX](http://cvxr.com/cvx/): A MATLAB library for convex optimization. Used for implementing methods from [Robust Principal Component Analysis?](https://dl.acm.org/citation.cfm?id=1970395) and [Robust PCA via Outlier Pursuit](https://arxiv.org/abs/1010.4237).
-[Fast SVD and PCA](https://www.mathworks.com/matlabcentral/fileexchange/47132-fast-svd-and-pca): Used for faster computation of SVD for some large matrices.
-[Novembre_etal_2008_misc](https://github.com/NovembreLab/Novembre_etal_2008_misc): A repository containing data from [Genes Mirror Geography in Europe](https://www.nature.com/articles/nature07331). Used for our semi-synthetic experiments. Code from `files/` should be extracted to the folder `robust-filter/Novembre_etal_2008_misc/`.

Explanation of Files
===
Code for estimating the mean of a Gaussian:
-`filterGaussianMean.m`: Our algorithm
-`ransacGaussianMean.m`: A RANSAC-based method
-`geoMedianGaussianMean.m`: Geometric median
-`pruneGaussianMean.m`: Coordinate-wise median, followed by naive pruning of distant points

Code for estimating the covariance of a Gaussian:
-`filterGaussianCov.m`: Our algorithm
    -`findMaxPoly.m`: Finds the structured degree-two polynomial which is maximized by the data
    -`flatten.m` and `sharpen.m`: Convert between matrix and vector representations
-`filterGaussianCovTuned.m`: A version of our algorithm which is tuned to select hyperparameters 
-`pruneGaussianCov.m`: Naive pruning of distant points
-`ransacMVE.m`: A RANSAC-based method
    -`MVE.m`: Approximating the MVE for a small dataset
-`ADPCP.m`: Principal Component Pursuit by Alternating Directions, from [Robust Principal Component Analysis?](https://dl.acm.org/citation.cfm?id=1970395)
    -`specThresh.m` and `shrinkage.m`: Singular value thresholding and shrinkage operators
    -`norm_nuc.m`: Compute the nuclear norm of a matrix
-`mahalanobis.m`: Compute Mahalanobis rescaling of a matrix

Test files, used for comparing methods for mean and covariance estimation
-`testMeanEstimators.m`: Compares mean estimation algorithms
-`testGeoMedian.m`: Demonstrates that the geometric median incurs an O(sqrt(d)) loss in accuracy
-`testCovEstimators.m`: Compares covariance estimation algorithms
-`processGenomicData.m`: Applies covariance estimation algorithms to semi-synthetic genome dataset

Reproducibility
===
Figures in the paper can be reproduced by running the following scripts:
-Figure 1: `testMeanEstimators.m`
-Figure 2: `testCovEstimators.m`
-Figures 3 and 4: `processGenomicData.m`

Reference
===
This repository is an implementation of our paper [Being Robust (in High Dimensions) Can Be Practical](https://arxiv.org/abs/1703.00893) from ICML 2017, authored by [Ilias Diakonikolas](http://www.iliasdiakonikolas.org/), [Gautam Kamath](http://www.gautamkamath.com/), [Daniel M. Kane](https://cseweb.ucsd.edu/~dakane/), [Jerry Li](http://www.mit.edu/~jerryzli/), [Ankur Moitra](http://people.csail.mit.edu/moitra/), and [Alistair Stewart](http://www.alistair-stewart.com/).

See also our original theory paper, [Robust Estimators in High Dimensions without the Computational Intractability](https://arxiv.org/abs/1604.06443), which appeared in FOCS 2016.
