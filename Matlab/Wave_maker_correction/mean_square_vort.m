%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
fname      = strcat(dirpath,'rip_his_out.nc');

makemovie = 0;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
%
%======================================================================

yindex = 101;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record

if makemovie,
 movObj = QTWriter('swash.mov');
 tstr=1;
 tend=tindex;
else,
 tstr=tindex;
 tend=tstr;
end

hf = figure('position',[1000 500 400 800]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

time=nc{'scrum_time'}(:);
tend=length(time);
pm=squeeze(nc{'pm'}(:));
pn=squeeze(nc{'pn'}(:));
xr=squeeze(nc{'x_rho'}(:));
yr=squeeze(nc{'y_rho'}(:));

for it=1:tend
    ubar=squeeze(nc{'ubar'}(it,:,:));
    vbar=squeeze(nc{'vbar'}(it,:,:));
    w(it,:,:)=vorticity(ubar,vbar,pm,pn);
end

[P,L]=size(yr);
w(:,:,2:L-1)=0.5*(w(:,:,1:L-2)+w(:,:,2:L-1));
w(:,:,1)=w(:,:,2);w(:,:,L)=w(:,:,L-1);
w(:,2:P-1,:)=0.5*(w(:,1:P-2,:)+w(:,2:P-1,:));
w(:,1,:)=w(:,2,:);w(:,P,:)=w(:,P-1,:);
Z=w.^2;
Zdy=trapz(yr(:,1),Z,2);
Zdxdy=trapz(xr(1,:),Zdy,3);

plot(time,Zdxdy);
title('Mean-square vorticity');
xlabel('t (s)');
ylabel('Z ($m^2s^{-2}$)','Interpreter','latex');