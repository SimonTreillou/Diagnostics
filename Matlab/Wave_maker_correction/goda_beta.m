function [res] = goda_beta(gamma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
res = 0.06238/(0.23 + 0.0336*gamma - 0.185*(1.9+gamma)^-1);
res = res * (1.094 - 0.01915*log(gamma));
end