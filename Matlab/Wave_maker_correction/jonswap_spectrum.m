function [w,S]=jonswap_spectrum(alpha,Nf,T,gamma,useRad)
    %alpha=parameter
    %Nf=number of frequencies
    %T=peak period
    %gamma=enhancement factor
    g=9.81;
    fm=1/T;
    wm=(2*pi)*fm;
    %useRad = 1; % 0=Hz or 1=rad/s
    if useRad
        wmin=wm*0.2; %CROCO parameters
        wmax=wm*2;   %CROCO parameters
        w=linspace(wmin,wmax,Nf);
        dw=(wmax-wmin)/Nf;
        sigma=0.07*(w<wm) + 0.09*(w>wm);
        b=exp(-1./(2.*sigma.^2) .* (w./wm - 1).^2);
        S= 0.3119 .* (wm)^4 .* (1./w.^5) .* exp((-5/4) .* (w./wm).^-4) .* gamma.^b;
        S=16*S*dw;
        %S= (alpha*g^2) .* (1./w.^5) .* exp((-5/4) .* (w./wm).^-4) .* gamma.^b;
    else
        fmin=fm*0.2;
        fmax=fm*2;
        f=linspace(fmin,fmax,Nf);
        dw=(fmax-fmin)/Nf;
        sigma=0.07*(f<fm) + 0.09*(f>fm);
        b=exp(-1./(2.*sigma.^2) .* (f./fm - 1).^2);
        S= 0.3119.* (fm)^4 .* (1./f.^5) .* exp((-5/4) .* (f./fm).^-4) .* gamma.^b;
        %S= (alpha*g^2).* (1./f.^5) .* exp((-5/4) .* (f./fm).^-4) .* gamma.^b;
        %S= (alpha*g^2)/(16*pi^4) .* (1./f.^5) .* exp((-5/4) .* (f./fm).^-4) .* gamma.^b;
        w=f;
        %S=16*S*dw;
    end
    %dw=w(2)-w(1);
    %S=S*16*(w(2)-w(1));
    %S=sqrt(S/sum(S));
    %S=S/sum(S*(dw));
end
