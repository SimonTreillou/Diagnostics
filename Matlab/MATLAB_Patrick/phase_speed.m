% program spectrum_latband_quick
%
clear all
close all
%=====================================
%
indir='../';
fname='rip_his_3D_9h.nc';
fname='rip_his.nc';
%fname='rip_his_mu006.nc';
%
imin0=110; imax0=imin0+20;

model3D=0;
%
%=====================================
fname=[indir,fname];
nc=netcdf(fname);

h=nc{'h'}(:);
xl=nc{'xl'}(:);
x=mean(nc{'x_rho'}(:,imin0:imax0),2)-xl+93; 
y=mean(nc{'y_rho'}(:,imin0:imax0),2);
N=length(nc('s_rho'));

tstr=30;
tend=length(nc{'scrum_time'}(:));

if model3D,
  v=squeeze(mean(squeeze(nc{'v'}(tstr:tend,N,:,imin0:imax0)),3));
else
  v=squeeze(mean(nc{'vbar'}(tstr:tend,:,imin0:imax0),3));
end
close(nc)
[M,L]=size(x);

v2=v.^2; 
s=zeros(size(v2));
s(v2<0.3*max(max(v2)))=1;
figure; imagesc(s); shading flat; colorbar;

dtheta=0.2; % 1
theta = 1:dtheta:179;
[R,xp] = radon(s,theta);

figure; pcolor(theta,xp,R); shading flat; colorbar

Rt = sum(R.^2, 1);
figure; plot(theta, Rt); axis([90 100 -Inf Inf]); grid on;

dt=60;
dy=3;
mint=80/dtheta; maxt=110/dtheta;
iang=find(Rt==max(Rt(mint:maxt)));
ang=(iang-1)*dtheta+1;
cp=tan(ang*pi/180)*dy/dt;
disp(['Phase speed = ',num2str(cp),' ---  Angle = ',num2str(ang)])
