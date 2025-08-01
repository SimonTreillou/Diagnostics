%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023// ';
fname      = strcat(dirpath,'rip_avg_out.nc');

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

zz=squeeze(nc{'zeta'}(:,:,1:ix0).^2);

PEdy=trapz(yr(:,1),0.5*g*zz,2);
PEdxdy=trapz(xr(1,1:ix0),PEdy,3);


hf = figure('position',[1000 500 400 800]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');
plot(time,PEdxdy);
title('Mean PE');
xlabel('t (s)');
ylabel('KE ($m^5s^{-2}$)','Interpreter','latex');