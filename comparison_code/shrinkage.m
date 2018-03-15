function [ M ] = shrinkage( M, tau )
%Shrinkage operator, as defined in CandesLMW'11
    [ilist,jlist,vlist] = find(M);
    for i = 1:length(ilist)
        v = vlist(i);
        M(ilist(i), jlist(i)) = sign(v) * max(abs(v) - tau, 0);
    end  
end