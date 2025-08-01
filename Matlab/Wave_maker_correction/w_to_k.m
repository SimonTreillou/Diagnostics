function [k] = w_to_k(w,h)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    g=9.81;
    K1=0.666666666;K2=0.3555555555;K3=0.1608465608;
    K4=0.063209876;K5=0.0217540484;K6=0.0065407983;
    
    khd=h.*(w).^2/g;
    kh=sqrt(khd.*khd+khd/(1.+khd.*(K1+khd.*(K2+khd.*(K3+khd.*(K4+khd.*(K5+K6.*khd)))))));
    k=kh./h;
end