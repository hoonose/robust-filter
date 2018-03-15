function [ estMean ] = geoMedianGaussianMean(data)
%Approximately compute the geometric median using Weiszfeld's algorithm
    numIters = 100;
    [N, ~] = size(data);
    curEstimate = mean(data);
    for i = 1:numIters
        num = 0;
        den = 0;
        for j = 1:N
            distToEstimate = norm(data(j,:) - curEstimate);
            num = num + data(j,:)/distToEstimate;
            den = den + 1/distToEstimate;
        end
        curEstimate = num/den;
    end
    estMean = curEstimate;
end