function [Sscaled] = mahalanobis(Shat, S)
%Computes Mahalanobis rescaling of Shat with respect to S
    Sscaled = inv(sqrtm(S))*Shat*inv(sqrtm(S));
end