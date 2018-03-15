function [ estMean ] = filterGaussianMean(data, eps, tau, cher)
%filterGaussianMean Run the filter algorithm on a Gaussian
%Suggested value of cher = 2.5
    [N, d] = size(data);
    empiricalMean = mean(data);
    threshold = eps*log(1/eps);
    centeredData = bsxfun(@minus, data, empiricalMean)/sqrt(N);

    [U, S, ~] = svdsecon(centeredData', 1);

    lambda = S(1,1)^2;
    v = U(:,1);

    %If the largest eigenvalue is about right, just return
    if lambda < 1 + 3 * threshold
       estMean = empiricalMean;
    %Otherwise, project in direction of v and filter
    else
        delta = 2*eps;
        projectedData1 = data * v;
        med = median(projectedData1);

        projectedData = [abs(data*v - med) data];

        sortedProjectedData = sortrows(projectedData);
        for i = 1:N
            T = sortedProjectedData(i,1) - delta;
            if (N - i) > cher * N * (erfc(T / sqrt(2)) / 2 + eps/(d*log(d*eps/tau)))
                break
            end
        end
        if i == 1 || i == N
            estMean = empiricalMean;
        else 
            estMean = filterGaussianMean(sortedProjectedData(1:i, 2:end), eps, tau, cher);
        end
    end
end