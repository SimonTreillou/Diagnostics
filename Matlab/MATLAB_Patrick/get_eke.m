clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%
fname     = '/Users/simon/Code/IB09/IB09PATPSdx1bathynewang11stratradiatif/rip_avg.nc';    % croco file name
%fname     = '../rip_avg_2D_SC.nc';
%fname     = '../rip_avg_2D_LC.nc';

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
%caxis([0 0.02])
%axis([-200 0 100 400])
shading flat;
colorbar
title('Mean longshore EKE')
xlabel('x (m)','Interpreter','latex')
ylabel('y (m)','Interpreter','latex')

mean(mean(meke))



