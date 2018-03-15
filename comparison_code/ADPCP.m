function [ L, S ] = ADPCP( M, lambda )
%Principal Component Pursuit by Alternating Directions
%Based on Algorithm 1 of "Robust Principal Component Analysis?" by
%Candes, Li, Ma, and Wright
    [m, n] = size(M);
    
    mu = m * n / (4 * norm_nuc(M));
    
    L = zeros(m, n);
    S = zeros(m, n);
    Y = zeros(m, n);
    
    
    delta = 1e-7;

    while norm(M - L - S, 'fro') > delta * norm(M, 'fro')
        L = specThresh(M - S + Y / mu, 1 / mu);
        S = shrinkage(M - L + Y / mu, lambda / mu);
        Y = Y + mu * (M - L - S);
    end
end