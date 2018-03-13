function [ estMean, filteredPoints] = pruneGaussianMean(data, eps)
%Use a naive method to prune points which are "obviously far" from mean
    [N, d] = size(data);
    coordMed = zeros(1,d);
    filteredPoints = zeros(size(data));
    numFilteredPoints = 0;
    
    for i = 1:d
        coordMed(i) = median(data(:,i));
    end
    
    for i = 1:N
        if norm(data(i,:) - coordMed) < 0.33*sqrt(d*log(N/eps))
            numFilteredPoints = numFilteredPoints + 1;
            filteredPoints(numFilteredPoints,:) = data(i,:);
        end
    end
    filteredPoints = filteredPoints(1:numFilteredPoints,:);
    estMean = mean(filteredPoints);
end