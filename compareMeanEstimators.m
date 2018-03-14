%Generates Figure 1 from paper
clear
eps = 0.1;
tau = 0.1;
cher = 2.5;

filterErr = [];
medianErr = [];
ransacErr = [];
LRVErr = [];
sampErr = [];
noisySampErr = []; 
prunedErr = [];
ds = 100:50:400;

for d = ds
    fprintf('d = %d\n', d)
    N = 10*floor(d/eps^2);
    sumFilterErr = 0;
    sumMedErr = 0;
    sumRansacErr = 0;
    sumLRVErr = 0;
    sumSampErr = 0;
    sumNoisySampErr = 0;
    sumPrunedErr = 0;

    X =  mvnrnd(zeros(1,d), eye(d), round((1-eps)*N)) + ones(round((1-eps)*N), d);

    fprintf('Sampling Error w/o noise\n');
    sumSampErr = sumSampErr + norm(mean(X) - ones(1,d));

    Y1 = randi([0 1], round(0.5*eps*N), d); 
    Y2 = [12*ones(round(0.5*eps*N),1), -2 * ones(round(0.5*eps*N), 1), zeros(round(0.5 * eps * N), d-2)];
    X = [X; Y1; Y2];

    fprintf('Sampling Error w noise\n');
    sumNoisySampErr = sumNoisySampErr + norm(mean(X) - ones(1,d));

    fprintf('Pruning\n');
    [prunedMean, ~] = pruneGaussianMean(X, eps);
    sumPrunedErr = sumPrunedErr + norm(prunedMean - ones(1, d));

    fprintf('Median\n')
    gm = geoMedianGaussianMean(X);
    sumMedErr = sumMedErr + norm(gm - ones(1, d));

    fprintf('Ransac\n')
    sumRansacErr = sumRansacErr + norm(ransacGaussianMean(X, eps, tau) - ones(1, d));

    fprintf('LRV\n')
    sumLRVErr = sumLRVErr + norm(agnosticMeanGeneral(X, eps) - ones(1,d));

    fprintf('Filter\n')
    sumFilterErr = sumFilterErr + norm(filterGaussianMean(X, eps, tau, cher) - ones(1, d));

    medianErr = [medianErr sumMedErr];
    ransacErr = [ransacErr sumRansacErr];
    filterErr = [filterErr sumFilterErr];
    LRVErr = [LRVErr sumLRVErr];
    sampErr = [sampErr sumSampErr];
    noisySampErr = [noisySampErr sumNoisySampErr];
    prunedErr = [prunedErr sumPrunedErr];
end

noisySampErr = noisySampErr - sampErr;
prunedErr = prunedErr - sampErr;
medianErr = medianErr - sampErr;
ransacErr = ransacErr - sampErr;
LRVErr = LRVErr - sampErr;
filterErr = filterErr - sampErr;

plot(ds, noisySampErr, ds, prunedErr, ds, medianErr, '-ro', ds, ransacErr, ds, LRVErr, ds, filterErr, '-.b', 'LineWidth', 2)
xlabel('Dimension')
ylabel('Excess L2 error')
legend('Sampling Error (with noise)', 'Naive Pruning', 'Geometric Median', 'RANSAC', 'LRV', 'Filter')