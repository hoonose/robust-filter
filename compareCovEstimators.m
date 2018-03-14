clear;
eps = 0.05;
tau = 0.1;

sampErr = [];
noisyEmpErr = [];
ransacErr = [];
filterErr = [];
pruneErr = [];
LRVErr = [];
ds = 10:10:100;

%Set to 0 for isotropic covariance, generates left half of Figure 2
%Set to 1 for spiked covariance, generates right half of Figure 2
spikedCovariance = 1;

if spikedCovariance
    spike = 100;
else
    spike = 1;
end


for d = ds
    N =  0.5*d / eps^2;
    fprintf('Training with dimension = %d, number of samples = %d \n', d, round(N, 0))
    sumEmpErr = 0;
    sumNoisyEmpErr = 0;
    sumFilterErr = 0;
    sumMVEErr = 0;
    sumLRVErr = 0;
    sumPruneErr = 0;

    covar = eye(d);
    if spikedCovariance
        covar(1,1) = spike;
    end

    X =  mvnrnd(zeros(1,d), covar, round((1-eps)*N)); 
    if spikedCovariance
        U1 = orth(randn(d, d));
        Y = [ 0.5 * randi([-1 1], round(eps *N), d / 2) 0.8 * randi([-2 2], round(eps *N), d / 2 - 1) randi([-spike spike], round(eps *N), 1)] * U1; 
    else
        Y = zeros(round(eps * N), d);
    end
    Z = [X; Y];

    fprintf('Sampling Error w/o noise...');
    empCov = cov(X);
    sumEmpErr = sumEmpErr + norm(mahalanobis(empCov, covar) - eye(d), 'fro');
    fprintf('done\n')

    fprintf('Sampling Error with noise...');
    empCov = cov(Z);
    sumNoisyEmpErr = sumNoisyEmpErr + norm(mahalanobis(empCov, covar) - eye(d), 'fro');
    fprintf('done\n')

    fprintf('Filter...')
    [ourCov, filterPoints, ~] = filterGaussianCovTuned(Z, zeros(size(Z)),  eps, tau, false);
    sumFilterErr = sumFilterErr + norm(mahalanobis(ourCov, covar) - eye(d), 'fro');
    fprintf('done\n')

    fprintf('Prune...')
    [pruneCov, prunePoints, ~] = pruneGaussianCov(Z, zeros(size(Z)),  tau, false);
    sumPruneErr = sumPruneErr + norm(mahalanobis(pruneCov, covar) - eye(d), 'fro');
    fprintf('done\n')


    fprintf('RANSAC...')
    sumMVEErr = sumMVEErr + norm(mahalanobis(ransacMVE(Z, eps, ceil(d / eps), 1000),  covar)  -  eye(d), 'fro');
    fprintf('done\n')

    fprintf('LRV...')
    [~, lrvCov, ~] = agnosticCovarianceGeneral(Z, eps);
    sumLRVErr = sumLRVErr + norm(mahalanobis(lrvCov, covar) - eye(d), 'fro') ;
    fprintf('done\n')


    sampErr = [sampErr sumEmpErr];
    noisyEmpErr = [noisyEmpErr sumNoisyEmpErr];
    ransacErr = [ransacErr sumMVEErr];
    filterErr = [filterErr sumFilterErr];
    pruneErr = [pruneErr sumPruneErr]; 
    LRVErr = [LRVErr sumLRVErr];
end

noisyEmpErr = noisyEmpErr - sampErr;
pruneErr = pruneErr - sampErr;
ransacErr = ransacErr - sampErr;
LRVErr = LRVErr - sampErr;
filterErr = filterErr - sampErr;

figure(1);
plot(ds, noisyEmpErr, '-gx',  ds,  ransacErr, ds, LRVErr, ds, pruneErr, ds, filterErr, '-.b', 'LineWidth', 2)
xlabel('Dimension')
ylabel('Excess Frobenius error')
legend('Samping Error (with corruption)', 'RANSAC', 'LRV', 'Prune', 'Filter')

figure(2);
plot(ds, LRVErr, ds, filterErr, '-.b', 'LineWidth', 2)
xlabel('Dimension')
ylabel('Excess Frobenius error')
legend('LRV', 'Filter')