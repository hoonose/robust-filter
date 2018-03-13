function [ estCov, filteredPoints, filteredMetadata] = pruneGaussianCov( data, metadata, tau, debug )
% Remove any points which are sufficiently far away from the empirical  covariance
    [N, d] = size(data);
    C = 0.2;
    
    empCov = data' * data / N;
    restart = false;
    remove = [];
    for i = 1:N
        if mod(i, 10000) == 0 && debug 
            fprintf('Initial pruning iteration %d\n',i);
        end
        x = data(i, :);
        if x * empCov^(-1) * x' > C * d * log(N / tau)
            remove = [i remove];
            restart = true;
        end
    end
    if restart
        if debug
            fprintf('Pruning %d points\n', length(remove)); 
        end
        data(remove, :) = [];
        metadata(remove, :) = [];
        [estCov, filteredPoints, filteredMetadata] = pruneGaussianCov(data, metadata, tau, debug);
        return;
    end
    estCov = empCov;
    filteredPoints = data;
    filteredMetadata = metadata;
    if debug
        fprintf('Remaining points: %d\n', size(filteredPoints, 1));
    end
end