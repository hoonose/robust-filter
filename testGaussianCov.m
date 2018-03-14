clear;
eps = 0.05;
tau = 0.1;

sampErr = [];
noisyEmpErr = [];
filterErr = [];
ds = 10:10:50;

%Set to 0 for isotropic covariance
%Set to 1 for spiked covariance
spikedCovariance = 1;

if spikedCovariance
    spike = 100;
else
    spike = 1;
end


for d = ds
    fprintf('d = %d\n', d)
    N =  0.5*d / eps^2;
    sumEmpErr = 0;
    sumNoisyEmpErr = 0;
    sumFilterErr = 0;

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

    fprintf('Sampling error without noise\n')
    empCov = cov(X);
    sumEmpErr = sumEmpErr + norm(mahalanobis(empCov, covar) - eye(d), 'fro');

    fprintf('Sampling error with noise\n')
    empCov = cov(Z);
    sumNoisyEmpErr = sumNoisyEmpErr + norm(mahalanobis(empCov, covar) - eye(d), 'fro');

    fprintf('Filter\n')
    [ourCov, filterPoints, ~] = filterGaussianCovTuned(Z, zeros(size(Z)),  eps, tau, false);
    sumFilterErr = sumFilterErr + norm(mahalanobis(ourCov, covar) - eye(d), 'fro');

    sampErr = [sampErr sumEmpErr];
    noisyEmpErr = [noisyEmpErr sumNoisyEmpErr];
    filterErr = [filterErr sumFilterErr];
end

noisyEmpErr = noisyEmpErr - sampErr;
filterErr = filterErr - sampErr;

figure(1);
plot(ds, noisyEmpErr, '-gx',  ds, filterErr, '-.b', 'LineWidth', 2)
xlabel('Dimension')
ylabel('Excess Frobenius error')
legend('Samping Error (with corruption)', 'Filter')