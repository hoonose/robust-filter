function [ Mthres ] = specThresh( M, tau )
%Singular Value Thresholding operator, as defined in CandesLMW'11
    [U, S, V] = svd(M);
    Mthres = U * shrinkage(S, tau) * V';
end