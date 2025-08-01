%clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%
fname='/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang7_spread10_init016_trc500/rip_avg.nc'
fname     = '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/rip_avgtot.nc';   % croco history file name
fname     = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_avgtot.nc';   % croco history file name
%fname = '/Users/simon/Code/CONFIGS/IB09_21AVR_spinup/rip_avg.nc';
%fname     = '../rip_avg_2D_SC.nc';
%fname     = '../rip_avg_2D_LC.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

makepdf   = 0;                       % make pdf file
%
%===============================================================
model3D=1;

nc=netcdf(fname);
tstr=5;
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1;  end;
%tend=88;

h=nc{'h'}(:);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:)-xl;
y=nc{'y_rho'}(:);
N=length(nc('s_rho'));

if model3D,
 u=squeeze(nc{'u'}(tstr:tend,N,:,:));
 v=squeeze(nc{'v'}(tstr:tend,N,:,:));
else
 u=nc{'ubar'}(tstr:tend,:,:);
 v=nc{'vbar'}(tstr:tend,:,:);
end

for it=1:tend-tstr+1
  ur(it,:,:)=u2rho_2d(squeeze(u(it,:,:)));
  vr(it,:,:)=v2rho_2d(squeeze(v(it,:,:)));
end
mu=mean(ur);
mv=mean(vr);

for it=1:tend-tstr+1
 eke(it,:,:)=0.5*((ur(it,:,:)-mu).^2+(vr(it,:,:)-mv).^2);
end
meke=squeeze(mean(eke));

figure
pcolor(x,y,meke)
%caxis([0 0.08])
%axis([-200 0 100 400])
shading flat;
colorbar
title('Mean longshore EKE')
xlabel('x (m)','Interpreter','latex')
ylabel('y (m)','Interpreter','latex')

mean(mean(meke))

%%
Jstr=1;
plot(mean(meke(Jstr:end,:)))
hold on
plot(mean(meke2(Jstr:end,:)))
ylim([0 0.2]);


