# Being Robust (in High Dimensions) Can Be Practical
A MATLAB implementation of [Being Robust (in High Dimensions) Can Be Practical](https://arxiv.org/abs/1703.00893) from ICML 2017.

Prerequisites 
===
This project requires installation of the following packages:
* [Fast SVD and PCA](https://www.mathworks.com/matlabcentral/fileexchange/47132-fast-svd-and-pca): Used for faster computation of SVD for some large matrices.
* [Novembre_etal_2008_misc](https://github.com/NovembreLab/Novembre_etal_2008_misc): A repository containing data from [Genes Mirror Geography in Europe](https://www.nature.com/articles/nature07331). Used for our semi-synthetic experiments. Extract the following three files to `robust-filter/genomicData`: 
    * `POPRES_08_24_01.EuroThinFinal.LD_0.8.exLD.out0-PCA.eigs`
    * `POPRES_08_24_01.EuroThinFinal.LD_0.8.exLD.out0-PCA.eval`
    * `POPRESID_Color.txt`

Optional
---
The following packages are optional, and are only used for comparison of our algorithm with alternative estimators:
* [AgnosticMeanAndCovarianceCode](https://github.com/kal2000/AgnosticMeanAndCovarianceCode): An implementation of the algorithms from the paper [Agnostic Estimation of Mean and Covariance](https://arxiv.org/abs/1604.06968) from FOCS 2016, authored by [Kevin A. Lai](https://www.cc.gatech.edu/~klai9/), [Anup B. Rao](https://sites.google.com/site/anupraob/), and [Santosh Vempala](https://www.cc.gatech.edu/~vempala/). Used to compare the statistical accuracy of our algorithms against theirs.
* [CVX](http://cvxr.com/cvx/): A MATLAB library for convex optimization. Used for implementing methods from [Robust Principal Component Analysis?](https://dl.acm.org/citation.cfm?id=1970395) and [Robust PCA via Outlier Pursuit](https://arxiv.org/abs/1010.4237).

Explanation of Files
===
This repository contains several algorithms -- we identify files relevant to our algorithm, and those which are only used for comparison with alternatives.

Our Filter algorithms' files
---
The following are files which are used in our algorithms' implementations, and can all be found in the `filter_code` subdirectory.
* `filterGaussianMean.m`: Our algorithm for estimating the mean of a Gaussian
* `filterGaussianCov.m`: Our algorithm for estimating the covariance of a Gaussian
    * `findMaxPoly.m`: Finds the structured degree-two polynomial which is maximized by the data
    * `flatten.m` and `sharpen.m`: Convert between matrix and vector representations
* `filterGaussianCovTuned.m`: A version of our covariance estimation algorithm which is tuned to select hyperparameters 
* `mahalanobis.m`: Compute Mahalanobis rescaling of a matrix

Test files, for demonstrating performance of our estimators for mean and covariance. These can bef ound in the `test_code` subdirectory:
* `testGaussianMean.m`: Tests our mean estimation algorithm
* `testGaussianCov.m`: Tests our covariance estimation algorithm
* `testGenomicData.m`: Tests our covariance estimation algorithm on semi-synthetic genome dataset

Other algorithms' files
---
The following are files which are used for comparison with other competing algorithms, and can be found in the `comparison_code` subdirectory:

Algorithms for estimating the mean of a Gaussian:
* `ransacGaussianMean.m`: A RANSAC-based method
* `geoMedianGaussianMean.m`: Geometric median
* `pruneGaussianMean.m`: Coordinate-wise median, followed by naive pruning of distant points

Algorithms for estimating the covariance of a Gaussian:
* `pruneGaussianCov.m`: Naive pruning of distant points
* `ransacMVE.m`: A RANSAC-based method
    * `MVE.m`: Approximating the MVE for a small dataset
* `ADPCP.m`: Principal Component Pursuit by Alternating Directions, from [Robust Principal Component Analysis?](https://dl.acm.org/citation.cfm?id=1970395)
    * `specThresh.m` and `shrinkage.m`: Singular value thresholding and shrinkage operators
    * `norm_nuc.m`: Compute the nuclear norm of a matrix

Comparison files, for evaluating and comparing the performance of estimators for mean and covariance. These can be found in the `comparison_code` subdirectory:
* `testGeoMedian.m`: Demonstrates that the geometric median incurs an O(sqrt(d)) loss in accuracy
* `compareMeanEstimators.m`: Compares mean estimation algorithms
* `compareCovEstimators.m`: Compares covariance estimation algorithms
* `compareGenomicData.m`: Compares covariance estimation algorithms to semi-synthetic genome dataset

Reproducibility
===
Figures in the paper can be approximately reproduced by running the following scripts:
* Figure 1: `compareMeanEstimators.m`
* Figure 2: `compareCovEstimators.m`
* Figures 3 and 4: `compareGenomicData.m`

Reference
===
This repository is an implementation of our paper [Being Robust (in High Dimensions) Can Be Practical](https://arxiv.org/abs/1703.00893) from ICML 2017, authored by [Ilias Diakonikolas](http://www.iliasdiakonikolas.org/), [Gautam Kamath](http://www.gautamkamath.com/), [Daniel M. Kane](https://cseweb.ucsd.edu/~dakane/), [Jerry Li](http://www.mit.edu/~jerryzli/), [Ankur Moitra](http://people.csail.mit.edu/moitra/), and [Alistair Stewart](http://www.alistair-stewart.com/).

See also our original theory paper, [Robust Estimators in High Dimensions without the Computational Intractability](https://arxiv.org/abs/1604.06443), which appeared in FOCS 2016.
