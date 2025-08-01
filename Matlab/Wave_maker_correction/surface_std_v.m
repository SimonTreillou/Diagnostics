clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%
%fname     = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx1_3D/rip_tot_avg.nc';    % croco file name
fname     = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/rip_avg_out.nc'
%fname     = '/Users/simon/Code/CONFIGS/IB09_randomphase/rip_avgtot.nc';   % croco history file name
%fname     = '/Users/simon/Code/CONFIGS/IB09_21AVR_spinup/rip_avg.nc';   % croco history file name
%fname     = '../rip_avg_2D_SC.nc';
%fname     = '../rip_avg_2D_LC.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

makepdf   = 0;                       % make pdf file
%
%===============================================================
model3D=1;

nc=netcdf(fname);
tstr=5 ;
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1;  end;
%tend=88;

h=nc{'h'}(:);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:)-xl;
y=nc{'y_rho'}(:);
N=length(nc('s_rho'));



v=nc{'vbar'}(tstr:tend,:,:);
for it=1:tend-tstr+1
  vr(it,:,:)=v2rho_2d(squeeze(v(it,:,:)));
end

mv=squeeze(std(vr));

figure
contourf(x,y,mv,10)
%caxis([-0.01 0.01])
% New colors for the COLORMAP example:
colormap(brewermap([],"RdPu"))
%axis([-200 0 100 400])
shading flat;
caxis
colorbar
title('Longshore velocity standard deviation')
xlabel('x (m)','Interpreter','latex')
ylabel('y (m)','Interpreter','latex')

std(std(mv))



