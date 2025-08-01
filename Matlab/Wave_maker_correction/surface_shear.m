%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = '/Users/simon/Code/IB09/IB09PATPSdx1bathynewang11stratradiatif/rip_avg.nc';   % croco history file name
fname     = '/Users/simon/Code/IB09/IB09PATPSdx1bathynewang7strat/rip_avg.nc';   % croco history file name
fname='/Users/simon/Code/testAVDV/testVADVintweno/rip_avg.nc';
fname='/Users/simon/Code/IB09/IB09_PSbienposang17deb3sansspread2sansang/rip_avg.nc';

fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';


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


vbar=squeeze(nc{'vbar'}(tindex,:,:));
ubar=squeeze(nc{'ubar'}(tindex,:,:));
pm=squeeze(nc{'pm'}(:));
pn=squeeze(nc{'pn'}(:));
w=vorticity(ubar,vbar,pm,pn);
xr=squeeze(nc{'x_rho'}(:));
yr=squeeze(nc{'y_rho'}(:));
[P,L]=size(yr);
w(:,2:L-1)=0.5*(w(:,1:L-2)+w(:,2:L-1));
w(:,1)=w(:,2);w(:,L)=w(:,L-1);

xi = shear(ubar,vbar,pm,pn);
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
colormap(rgb);

w(2:P-1,:)=0.5*(w(1:P-2,:)+w(2:P-1,:));
w(1,:)=w(2,:);w(P,:)=w(P-1,:);

xi(:,2:L-1)=0.5*(xi(:,1:L-2)+xi(:,2:L-1));
xi(:,1)=xi(:,2);xi(:,L)=xi(:,L-1);
xi(2:P-1,:)=0.5*(xi(1:P-2,:)+xi(2:P-1,:));
xi(1,:)=xi(2,:);xi(P,:)=xi(P-1,:);

cmin=-0.01;
cmax=0.01;
cint=0.001;
h=pcolor(xr(:,:)-300,yr(:,:),w(:,:));
set(h,'EdgeColor','none');
caxis([cmin cmax]);
colorbar();
