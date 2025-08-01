%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023// ';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/';   % croco history file name
fname      = strcat(dirpath,'rip_avgtot.nc');

makemovie = 0;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
model3D   = 1;
%
%======================================================================

yindex = 101;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record

time=nc{'scrum_time'}(:);
tend=length(time);

N=10;
xr=squeeze(nc{'x_rho'}(:));
yr=squeeze(nc{'y_rho'}(:));
[~,ix0] = min(abs(xr(1,:)-300+50));

if model3D,
 u=squeeze(nc{'u'}(:,N,:,1:ix0));
 v=squeeze(nc{'v'}(:,N,:,1:ix0+1));
else
 u=nc{'ubar'}(:,:,1:ix0);
 v=nc{'vbar'}(:,:,1:ix0);
end

for it=1:tend
  ur(it,:,:)=u2rho_2d(squeeze(u(it,:,:)));
  vr(it,:,:)=v2rho_2d(squeeze(v(it,:,:)));
end
mu=mean(ur);
mv=mean(vr);

for it=1:tend
 eke(it,:,:)=0.5*((ur(it,:,:)).^2+(vr(it,:,:)).^2);
end

hf = figure('position',[1000 500 400 800]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

time=nc{'scrum_time'}(:);
tend=length(time);
pm=squeeze(nc{'pm'}(:));
pn=squeeze(nc{'pn'}(:));

KEdy=trapz(yr(:,1),eke,2);
KEdxdy=trapz(xr(1,1:ix0+1),KEdy,3);

plot(time,KEdxdy);
title('Mean KE');
xlabel('t (s)');
ylabel('KE ($m^5s^{-2}$)','Interpreter','latex');