function [ estCov, filteredPoints, filteredMetadata] = filterGaussianCov( data, metadata, eps, tau, expConst, debug )
    MAX_COND = 10000;

    [N, d] = size(data);
    threshold = eps*(log(1/eps))^2;
    C1 = 0.4; 
	C2 = 0;
    
    empCov = data' * data / N;
    condition = cond(empCov);
    
    % Our empirical covariance is too ill-conditioned--we probably threw
    % away too many points, so let's just throw everything away and restart
    if condition > MAX_COND || N < d
        if debug
            fprintf('Ill conditioned %d %d %d\n', condition, N, d);
        end
        estCov = zeros(d);
        filteredPoints = [];
        filteredMetadata = [];
        return
    end
        
    
    empCovInv = empCov^(-1);
    % Remove any points which are very far away from the empirical
    % covariance
    restart = false;
    remove = [];
    for i = 1:N
        if mod(i, 10000) == 0 && debug 
            fprintf('Initial pruning iteration %d\n',i);
        end
        x = data(i, :);
        if x * empCovInv * x' > C1 * d * log(N / tau)
            remove = [i remove];
            restart = true;
        end
    end
    if restart
        data(remove, :) = [];
        metadata(remove, :) = [];
        [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(data, metadata, eps, tau, expConst, debug);
        return
    end
        
    if debug
        fprintf('After first for loop\n')
    end
    
    
    rotData = data * empCovInv^(1/2);
    if debug
        fprintf('before findMaxPoly\n');
    end
    [~, M, lambda] = findMaxPoly(rotData);
    if debug
        fprintf('after findMaxPoly\n');
    end
    
    if debug
        fprintf('lambda = %d, threshold = %d\n', lambda, 1 + C2 * threshold);
        if d < 12
            svd(M)
        end
    end
    
    %If top eigenvalue isn't too large, the dataset is OK
    %Otherwise, find points which violate the tail bound
    if lambda < 1 + C2 * threshold
        estCov = empCov;
        filteredPoints = data;
        filteredMetadata = metadata;
    else
        projectedData = zeros(N, 1);
        for i = 1:N
            projectedData(i) = rotData(i, :) * M * rotData(i, :)';
        end
                
        med = mean(projectedData);
		
        if debug
            fprintf('mean %d, trace %d, variance %d , eigenvalue %d \n', med, trace(M), var(projectedData), lambda);
            fprintf('N %d \n',N);
        end
		
        indices = abs(projectedData - med);
        sortedProjectedData = sortrows([indices data]);
        sortedMetadata = sortrows([indices metadata]);
        
        for i = 1:N
            T = sortedProjectedData(i,1);
            rhs = N * (12 * exp(-expConst*T) + 3 * eps/(d*log(N/tau))^2);
            if  debug && ((N - i) > rhs || mod(N - i, 1000) == 10)
                fprintf('N - i =  %d \n', N - i);
                fprintf('other thing =  %d \n', rhs);
                fprintf('T = %d\n', T);
                fprintf('first term RHS %d\n', N * 3 * exp(-T));
                fprintf('second term RHS %d\n', 3 * N * eps/(d*log(N/tau))^2);
            end
            if (N - i) > rhs
                break
            end
        end
        if i == 1 || i == N
            if debug
                fprintf('No data points removed, i = %d\n', i);
            end
            estCov = empCov;
            filteredPoints = data;
            filteredMetadata = metadata;
            return;
        else 
            [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(sortedProjectedData(1:i, 2:end), sortedMetadata(1:i,2:end), eps, tau, expConst, debug);
        end
    end
end