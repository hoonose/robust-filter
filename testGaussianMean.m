clear;
eps = 0.1;
tau = 0.1;
cher = 2.5;

filterErr = [];
sampErr = [];
noisySampErr = []; 
ds = 100:50:400;

for d = ds
    N = 10*floor(d/eps^2);
    fprintf('Training with dimension = %d, number of samples = %d \n', d, round(N, 0))
    sumFilterErr = 0;
    sumSampErr = 0;
    sumNoisySampErr = 0;

    X =  mvnrnd(zeros(1,d), eye(d), round((1-eps)*N)) + ones(round((1-eps)*N), d);

    fprintf('Sampling Error w/o noise...');
    sumSampErr = sumSampErr + norm(mean(X) - ones(1,d));
    fprintf('...done\n')
    
    Y1 = randi([0 1], round(0.5*eps*N), d); 
    Y2 = [12*ones(round(0.5*eps*N),1), -2 * ones(round(0.5*eps*N), 1), zeros(round(0.5 * eps * N), d-2)];
    X = [X; Y1; Y2];

    fprintf('Sampling Error with noise');
    sumNoisySampErr = sumNoisySampErr + norm(mean(X) - ones(1,d));
    fprintf('...done\n')
    
    fprintf('Filter')
    sumFilterErr = sumFilterErr + norm(filterGaussianMean(X, eps, tau, cher) - ones(1, d));
    fprintf('...done\n')
    
    filterErr = [filterErr sumFilterErr];
    sampErr = [sampErr sumSampErr];
    noisySampErr = [noisySampErr sumNoisySampErr];
end

noisySampErr = noisySampErr - sampErr;
filterErr = filterErr - sampErr;

plot(ds, noisySampErr, ds, filterErr, '-.b', 'LineWidth', 2)
xlabel('Dimension')
ylabel('Excess L2 error')
legend('Sampling Error (with noise)', 'Filter')
