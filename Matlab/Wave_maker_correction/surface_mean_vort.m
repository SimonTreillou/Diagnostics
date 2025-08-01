%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = '/Users/simon/Code/IB09/IB09PATPSdx1bathynewang11stratradiatif/rip_avg.nc';   % croco history file name
%fname = '/Users/simon/Code/IB09/testVADVC2periodok/rip_avg.nc';
%fname     = '/Users/simon/Code/CONFIGS/IB09_randomphase/rip_avgtot.nc';   % croco history file name
%fname = '/Users/simon/Code/CONFIGS/IB09_21AVR_spinup/rip_his.nc';
fname = '/Users/simon/Code/CONFIGS/testLD_corr_testforper/rip_his.nc';

fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
%fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

%fname     = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx1_tr1000_w011_sansinit_suite3/rip_avg.nc';   % croco history file name
%fname='/Users/simon/Code/testAVDV/testVADVsans/rip_avg.nc';

makemovie = 1;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
%
%======================================================================

yindex = 101;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tend=length(nc{'scrum_time'}(:)); % reads last record
tstr=5;

hf = figure('position',[1000 500 400 800]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');

rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   255
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;


ubar=squeeze(nc{'ubar'}(tstr,:,:,10));
vbar=squeeze(nc{'vbar'}(tstr,:,:,10));
pm=squeeze(nc{'pm'}(:));
pn=squeeze(nc{'pn'}(:));

w=vorticity(ubar,vbar,pm,pn);
xr=squeeze(nc{'x_rho'}(:))-300;
yr=squeeze(nc{'y_rho'}(:));
[P,L]=size(yr);
w(:,2:L-1)=0.5*(w(:,1:L-2)+w(:,2:L-1));
w(:,1)=w(:,2);w(:,L)=w(:,L-1);
w(2:P-1,:)=0.5*(w(1:P-2,:)+w(2:P-1,:));
w(1,:)=w(2,:);w(P,:)=w(P-1,:);
wTot=w;


for tindex=tstr+1:tend % ---------------------------------------------
    ubar=squeeze(nc{'ubar'}(tindex,:,:,10));
    vbar=squeeze(nc{'vbar'}(tindex,:,:,10));
    pm=squeeze(nc{'pm'}(:));
    pn=squeeze(nc{'pn'}(:));

    w=vorticity(ubar,vbar,pm,pn);
    xr=squeeze(nc{'x_rho'}(:))-300;
    yr=squeeze(nc{'y_rho'}(:));
    [P,L]=size(yr);
    w(:,2:L-1)=0.5*(w(:,1:L-2)+w(:,2:L-1));
    w(:,1)=w(:,2);w(:,L)=w(:,L-1);
    w(2:P-1,:)=0.5*(w(1:P-2,:)+w(2:P-1,:));
    w(1,:)=w(2,:);w(P,:)=w(P-1,:);

    wTot=w+wTot;
end

mw = wTot/ (tend-tstr);

cmin=-0.04;
%cmin=-max(abs(mw(:)));
cmax=-cmin;
cint=0.005;
h=pcolor(xr(:,:),yr(:,:),mw(:,:));
set(h,'EdgeColor','none');
caxis([cmin cmax]);
colormap(rgb);
colorbar();
ylim([min(yr(:,1)),max(yr(:,1))]);
xlim([min(xr(1,:)),max(xr(1,:))]);

%%
plot(mean(mw(:,:),1))
