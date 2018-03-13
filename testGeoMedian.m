%Demonstrates that the geometric median incurs a sqrt(d) loss in accuracy
eps = 0.5;
dims = floor(linspace(1,700,10));
err = [];
for d = dims
    fprintf('In dimension %f\n', d);
    N = ceil(10*(d/eps^2));
    sumErr = 0;
    for i = 1:3
        X =  mvnrnd(zeros(1,d), eye(d), round((1-eps)*N)) + ones(round((1-eps)*N), d); 
        Y = zeros(round(eps * N), d);
        X = [X; Y];
        sumErr = sumErr + norm(geoMedianGaussianMean(X) - ones(1, d));
    end
    err = [err sumErr / 3];
end
plot(dims, err)