

function [cond,beta] = condition(h)
    beta=h./180;
    L=2*pi ./w_to_k(2*pi/8,h);
    cond= (1+6.4 .*beta) .* h./L;
    return cond,beta
end