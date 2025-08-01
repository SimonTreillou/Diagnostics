function [theta,AG] = dirspread_salatin(thetam,sigma,Ndir,rad)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~rad
    thetam=thetam*pi/180;
    sigma=sigma*pi/180;
end

Ns=20/sigma;

theta=linspace(-pi/2,pi/2,Ndir)+thetam;
AG=zeros(Ndir,1);
for kf=1:Ndir
    AG(kf)=1/(2*pi);
    for kn=1:Ns
        AG(kf)=AG(kf)+1/pi*exp(-0.5*(kn*sigma)^2)*cos(kn*(theta(kf)-thetam*pi/180));
    end
end
end