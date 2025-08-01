%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Differences between spread at dx=2m
% Simon Treillou, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

dirpath    = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang7_spread10/'
fname      = strcat(dirpath,'rip_avg.nc');
nc=netcdf(fname,'r');
xr=squeeze(nc{'x_rho'}(1,:));
x=xr-(450-50);
y=squeeze(nc{'y_rho'}(:,1));
h=squeeze(nc{'h'}(1,:));
[dyr,ixyr] = min(abs(y(:,1)-248));
[dyr,ix0] = min(abs(x(1,:)));
[dyr,ixb] = min(abs(x(1,:)+181));
[dyr,ix150] = min(abs(x(1,:)+150));
load HR16-Figures/HR16_fig3_Vmod.mat
load HR16-Figures/HR16_fig3_Vobs.mat
V=squeeze(mean(mean(nc{'vbar'}(10:end,:,:),1),2));
Vstd=squeeze(std(mean(nc{'vbar'}(10:end,:,:),2),1));
mean(Vstd)
figure(1);
%plot(x(ix150:ix0)',V(ix150:ix0));
f = fill([x(ix150:ix0)';flipud(x(ix150:ix0)')],[V(ix150:ix0)-Vstd(ix150:ix0); ...
    flipud(V(ix150:ix0)+Vstd(ix150:ix0))],[.5 .9 .9],'linestyle','none');
hold on
alpha(f,0.5)




dirpath    = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang7_spread30/'
fname      = strcat(dirpath,'rip_avg.nc');
nc=netcdf(fname,'r');
xr=squeeze(nc{'x_rho'}(1,:));
x=xr-(450-50);
y=squeeze(nc{'y_rho'}(:,1));
h=squeeze(nc{'h'}(1,:));
[dyr,ixyr] = min(abs(y(:,1)-248));
[dyr,ix0] = min(abs(x(1,:)));
[dyr,ixb] = min(abs(x(1,:)+181));
[dyr,ix150] = min(abs(x(1,:)+150));
load HR16-Figures/HR16_fig3_Vmod.mat
load HR16-Figures/HR16_fig3_Vobs.mat
V=squeeze(mean(mean(nc{'vbar'}(10:end,:,:),1),2));
Vstd=squeeze(std(mean(nc{'vbar'}(10:end,:,:),2),1));
mean(Vstd)
%plot(x(ix150:ix0)',V(ix150:ix0));
f = fill([x(ix150:ix0)';flipud(x(ix150:ix0)')],[V(ix150:ix0)-Vstd(ix150:ix0); ...
    flipud(V(ix150:ix0)+Vstd(ix150:ix0))],[.9 .5 .9],'linestyle','none');
alpha(f,0.5);


ylim([-0.1,0.5]);
scatter(HR16_fig3_Vobs(:,1),HR16_fig3_Vobs(:,2),'Marker','*');
ylabel('V ($m.s^{-1}$)','Interpreter','latex');
xlabel('$x$ (m)','Interpreter','latex');
legend('DS=10','DS=30','observations', ...
    'Interpreter','latex','Location','northwest');

