function [ norm ] = norm_nuc( A )
%norm_nuc Compute the nuclear norm (or Schatten-1 norm) of a matrix A
    [~, S, ~] = svdecon(A);
    norm = trace(abs(S));
end