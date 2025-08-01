%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SHOREFACE test case
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
%       Application of the CROCO embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%  Patrick Marchesiello, IRD 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname3= '/Users/simon/Code/IB09/IB09_fini/rip_avg.nc';    % croco file name
fname1='/Users/simon/Code/IB09/IB09PATPSdx1bathynewang7strat/rip_avg.nc';
fname1='/Users/simon/Code/testAVDV/testVADVintweno/rip_avg.nc';

fname2='/Users/simon/Code/CONFIGS/CALMIP/IB09_dx1_2D/IB09_dx1_tr1000_w011_sansinit_DB10.2_suite2/rip_avg.nc';
fname1='/Users/simon/Code/IB09/IB09_PSbienposang17deb3sansspread2sansang/rip_avg.nc';
fname1     = '/Users/simon/Code/CONFIGS/IB09_randomphase/rip_histot.nc';   % croco history file name

fname1 = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname1 = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

makemovie = 0;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
%
%======================================================================

yindex = 160;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------
fname=fname1;
nc=netcdf(fname);
tstr=5;
tend=length(nc{'scrum_time'}(:)); % reads last record
%tindex=5; 

hf = figure('position',[1000 500 800 400]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

hr=squeeze(nc{'h'}(yindex,:));
L=length(hr);
xr=squeeze(nc{'x_rho'}(yindex,:));
yr=squeeze(nc{'y_rho'}(yindex,:));

dx=xr(2)-xr(1);
xmin=min(xr);
xmax=max(xr);
Dcrit=nc{'Dcrit'}(:);

N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(end,yindex,:));
zr=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'r',2));
zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
dzr=zr(2:end,:)-zr(1:end-1,:);         % ---> zw(2:N,:)
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
dzru=zru(2:end,:)-zru(1:end-1,:);      % ---> zwu(2:N,:)
dzw=zw(2:end,:)-zw(1:end-1,:);         % ---> zr
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dzwu=zwu(2:end,:)-zwu(1:end-1,:);      % ---> zru
%
xr2d=repmat(xr,[N 1]);
D=hr+zeta;
D2d=repmat(D,[N 1]);
%
%xmin= min(xr); xmax=90;
zmin=-max(hr); zmax=0.2*max(hr);


u=squeeze(nc{'u'}(tstr:tend,:,:,:));
v=squeeze(nc{'v'}(tstr:tend,:,:,:));

for it=1:(tend-tstr+1)
  ur(it,:,:,:)=u2rho_3d(squeeze(u(it,:,:,:)));
  vr(it,:,:,:)=v2rho_3d(squeeze(v(it,:,:,:)));
end

mu=mean(ur);
mv=mean(vr);

for it=1:(tend-tstr+1)
    eke(it,:,:,:)=0.5*((ur(it,:,:,:)-mu).^2+(vr(it,:,:,:)-mv).^2);
end

meke=squeeze(mean(mean(eke,3)));

cmin=min(min(meke)); cmax=max(max(meke)); nbcol=20;
cint=(cmax-cmin)/nbcol;

Dcrit=Dcrit+1.e-6;
var(D2d<Dcrit)=NaN;

ztop=squeeze(zw(end,:,:));
ztop(D<Dcrit+0.01)=NaN;

map=colormap(parula(nbcol));
%map(nbcol/2  ,:)=[1 1 1];
%map(nbcol/2+1,:)=[1 1 1];
%map=[1 0 0
%0 0 1];
%map=colormap(redwhiteblue(-1.5,3,100));
cmax=0.026
time=nc{'scrum_time'}(end);
contourf(xr2d-xr(1,end)+50,zr,meke,[cmin:cint:cmax],'LineStyle','none'); hold on
colorbar;
plot(xr-xr(1,end)+50,-hr,'color','k','LineWidth',3);
hn=plot(xr-xr(1,end)+50,ztop,'color','r','LineWidth',2);
grid on
axis([-100 15 -3 0.5])
caxis([cmin cmax])
tmin=floor(time/60);
title(['IB09: time and longshore-mean EKE at ',num2str(time/3600),' h'])
hold off
set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');

if makepdf
 export_fig -transparent kumar.pdf
end

return