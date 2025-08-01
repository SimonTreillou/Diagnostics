function [S,f] = mycspd(x,y,nfft,fs)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    % Length of analyzed signal
    len = length(x);

    % Window function definition
    ham=hamming(fix(nfft/2)*2);
    nw = length(ham);

    % Segmentation
    nseg=fix(len/nw);
    Stmp = zeros(nfft,1);
    for iseg=0:nseg-1
        ind=nw*iseg+[1:nw];
        xw = ham.*x(ind);
        yw = ham.*y(ind);
        Px = fft(xw,nfft);
        Py = fft(yw,nfft);
        Pxy = Py.*conj(Px);
        Stmp = Stmp + Pxy;
    end
    nfac=(fs*nseg*norm(ham)^2);
    S=[Stmp(1); 2*Stmp(2:nfft/2); Stmp(fix(nfft/2)+1)]/nfac;
    f=(fs/nfft)*[0:nfft/2]';
end