%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot of (u,w) horizontal vorticity.
%  To clear and adapt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = '/Users/simon/Code/IB09/IB09_PSbienposT2/rip_avg.nc';   % croco history file name
fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
%fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';
%
%======================================================================

yindex = 50;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
time=nc{'scrum_time'}(:);
tindex=length(nc{'scrum_time'}(:)); % reads last record
yindex=400;
%tindex=5; 


hf = figure('position',[1000 500 800 400]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');


hr=squeeze(nc{'h'}(yindex,:));
L=length(hr);
pm=nc{'pm'}(:);
xr=squeeze(nc{'x_rho'}(yindex,:));
yr=squeeze(nc{'y_rho'}(yindex,:));
dx=xr(2)-xr(1);
xmin=min(xr);
xmax=max(xr);
Dcrit=nc{'Dcrit'}(:);
% vertical grid
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
zr=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'r',2));
zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
dzr=zr(2:end,:)-zr(1:end-1,:);         % ---> zw(2:N,:)
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
dzru=zru(2:end,:)-zru(1:end-1,:);      % ---> zwu(2:N,:)
dzw=zw(2:end,:)-zw(1:end-1,:);         % ---> zr
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dzwu=zwu(2:end,:)-zwu(1:end-1,:);      % ---> zru
xr2d=repmat(xr,[N 1]);
D=hr+zeta;
D2d=repmat(D,[N 1]);
zmin=-max(hr); zmax=0.2*max(hr);


u=squeeze(nc{'u'}(tindex,:,yindex,:));
dudz=u./dzwu;


wtmp=squeeze(nc{'w'}(tindex,:,yindex,:));
wtmp2=(wtmp(2:end,:)+wtmp(1:end-1,:))/2;
w=(wtmp2(:,2:end)+wtmp2(:,1:end-1))/2;

dx=1/pm(1,1);
dwdx=w/dx;
[P,Q]=size(dwdx);
%dwdx=dwdx(2:P,2:Q)-dwdx(1:P-1,1:Q-1);

vort=dudz-dwdx;

l=max(min(min(vort)),max(max(vort)));
cmin=-l;
cmax=l;
nbcol=20;
cint=(cmax-cmin)/nbcol;

Dcrit=Dcrit+1.e-6;
vort(D2d(:,2:end)<Dcrit)=NaN;

ztop=squeeze(zw(end,:,:));
ztop(D<Dcrit+0.01)=NaN;

map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
map=[1 0 0
0 0 1];

contourf(xr2d(:,2:end)-xr(1,end)+50,zr(:,2:end),vort,[cmin:cint:cmax],'LineStyle','none'); hold on
colorbar;
plot(xr-xr(1,end)+50,-hr,'color','k','LineWidth',3);
hn=plot(xr-xr(1,end)+50,ztop,'color','r','LineWidth',2);
grid on
axis([xmin-xr(1,end)+50 xmax-xr(1,end)+50 zmin zmax])
caxis([cmin cmax])
tmin=floor(time/60);
title(['IB09: (u,w) vorticity at ',num2str(time(tindex)/3600),' h']);
xlabel('x (m)','Interpreter','latex');
ylabel('z (m)','Interpreter','latex');
hold off
set(gca,'fontsize',15);
set(gcf,'PaperPositionMode','auto');

return