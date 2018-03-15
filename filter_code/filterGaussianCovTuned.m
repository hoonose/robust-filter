function [ estCov, filteredPoints, filteredMetadata] = filterGaussianCovTuned( data, metadata, eps, tau, debug )
    expConstUB = inf;
    expConstLB = 0;
    expConstCand = 2;
    N = size(data, 1);    
    targetNoise = eps*N;
    
    if debug
        fprintf('Number of points: %d\n', N);
        fprintf('Target number of noise points is %d\n', targetNoise)
        fprintf('Trying to filter with constant %d\n\n\n', expConstCand);
    end
    
    
    [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(data, metadata, eps, tau, expConstCand, debug);
    outN = size(filteredPoints, 1);
    
    if debug
        fprintf('Number of points output by filter is %d\n', outN);
    end
    
    %Find upper and lower bounds for removing the right number of points
    if outN < (N - 1.5*targetNoise)
        while outN < (N - 1.5*targetNoise)
            expConstUB = expConstCand;
            expConstCand = expConstCand / 2;
            
            if debug
                fprintf('Trying to filter with constant %d\n\n\n', expConstCand);
            end
            [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(data, metadata, eps, tau, expConstCand, debug);
            outN = size(filteredPoints, 1);
            
            if debug
                fprintf('Number of points output by filter is %d\n', outN);
            end
        end
        expConstLB = expConstCand;
    elseif outN > (N - 0.5*targetNoise)
        while outN > (N - 0.5*targetNoise)
            expConstLB = expConstCand;
            expConstCand = expConstCand * 2;
            if debug
                fprintf('Trying to filter with constant %d\n\n\n', expConstCand);
            end
            [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(data, metadata, eps, tau, expConstCand, debug);
            outN = size(filteredPoints, 1);
            if debug
                fprintf('Number of points output by filter is %d\n', outN);
            end
        end
        expConstUB = expConstCand;
    end
    if outN > (N - 1.5*targetNoise) && outN < (N - 0.5*targetNoise)
        if debug
            fprintf('Succeeded with constant %d\n', expConstCand);
        end
        return;
    end
    
    %Binary search within the interval to remove the right number of points
    numTrials = 0;
    while 1
        numTrials = numTrials + 1;
        if numTrials >= 10
            if debug
                fprintf('Quitting due to too many trials\n');
            end
            [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(data, metadata, eps, tau, expConstUB, debug);
            return;
        end
        
        expConstCand = (expConstUB + expConstLB)/2;
        if debug
            fprintf('Trying to filter with constant %d\n\n\n', expConstCand);
        end
        [estCov, filteredPoints, filteredMetadata] = filterGaussianCov(data, metadata, eps, tau, expConstCand, debug);
        outN = size(filteredPoints, 1);
        
        if debug
            fprintf('Number of points output by filter is %d\n', outN);
        end
        
        if outN < (N - 1.5*targetNoise)
           expConstUB = expConstCand; 
        elseif outN > (N - 0.5*targetNoise)
            expConstLB = expConstCand;
        else
            if debug
                fprintf('Succeeded with constant %d\n', expConstCand);
            end
            return;
        end
    end
end