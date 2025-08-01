function [w] = k_to_w(k,h)  
    g=9.81;
    w=sqrt(g.*k.*tanh(k.*h));
end