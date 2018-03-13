function [ c, M, lambda ] = findMaxPoly(data)
% Given a dataset and an empirical covariance matrix, outputs the
% polynomial 
%           p*(x) = x^T M x - c
% with unit variance

    [N, d] = size(data);
    v = flatten(eye(d));
    
    dataKron = zeros(N, d * d);
    for i = 1:N
         dataKron(i, :) = kron(data(i, :), data(i, :));
    end
    
    empFourth = dataKron' * dataKron / N - v' * v;
    [U, S, ~] = svdsecon(empFourth, 1);
    
    lambda = S(1, 1) / 2;
    Mflat = U(:, 1);

    M = sharpen(Mflat, d) / sqrt(2);
    c = trace(M) / sqrt(2);
	
end
