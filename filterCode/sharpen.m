function M = sharpen(v,d)
% Given a d^2 dimensional vector v, convert to a d by d matrix M
    M = zeros(d, d);
    for i = 1:d
        for j = 1:d
            M(i, j) = v((i-1)*d + j);
        end
    end
end