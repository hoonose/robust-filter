function [ estCov, vol ] = MVE( samp, eps )
    [N, d] = size(samp);

    empCov = samp' * samp / N;
    empCovInv = empCov^(-1);
    mDists = zeros(N, 1);
    for j = 1:N
        mDists(j) = samp(j, :) * empCovInv * samp(j, :)';
    end
    mDistsSorted = sort(mDists);
    lambda = mDistsSorted(N - floor(eps * N));
    estCov = lambda * empCov;
    vol = det(empCov) * lambda^d;
end