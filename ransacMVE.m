function [ estCov ] = ransacMVE( data, eps, sampSize, iter )
% Use RANSAC to estimate the Minimum Volume Ellipsoid (MVE) as per 
% http://webmining.spd.louisville.edu/wp-content/uploads/2014/05/A-Brief-Overview-of-Robust-Statistics.pdf

    [~, d] = size(data);
    
    estCov = zeros(d, d);
    bestVol = Inf;
    
    for i = 1:iter
        if mod(i, 10000) == 0
            fprintf('ransac iteration i = %d\n', i)
        end
        samp = datasample(data, sampSize);
        [currentCov, currentVol] = MVE(samp, eps);
        if currentVol < bestVol
            estCov = currentCov;
        end
    end

    estCov = estCov / chi2inv(1 - eps, d);
    
end