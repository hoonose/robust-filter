function [ estMean ] = ransacGaussianMean( data, eps, tau)
%Use RANSAC to estimate the mean of Gaussian data
    [N, d] = size(data);
    empiricalMean = mean(data);
    ransacN = ceil(2*(d*log(4) + log(2/tau))/eps^2);
    if ransacN > N
        estMean = empiricalMean;
        return;
    end

    numIters = 100;
    thresh = d + 2*(sqrt(d * log(N/tau)) + log(N/tau)) + eps^2*(log(1/eps))^2;

    %As a baseline, see how many inliers the mean has
    bestInliers = 0;
    bestMean = empiricalMean;
    for j = 1:N
        if norm(empiricalMean - data(j))^2 <= thresh
            bestInliers = bestInliers + 1;
        end
    end

    for i = 1:numIters
        ransacData = data(randsample(1:N, ransacN),:);
        ransacMean = mean(ransacData);
        curInliers = 0;
        for j = 1:N
            if norm(ransacMean - data(j))^2 <= thresh
                curInliers = curInliers + 1;
            end
        end
        if curInliers > bestInliers
            bestMean = ransacMean;
            bestInliers = curInliers;
        end
    end

    estMean = bestMean;
end