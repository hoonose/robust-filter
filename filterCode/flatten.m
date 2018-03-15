function [ v ] = flatten( mat )
% Given a d by d matrix M, convert to a d^2 dimensional vector v
    [d, ~] = size(mat);
    v = zeros(1, d * d);
    for i = 1:d
        for j = 1:d
            v((i - 1) * d + j) = mat(i, j);
        end
    end
end