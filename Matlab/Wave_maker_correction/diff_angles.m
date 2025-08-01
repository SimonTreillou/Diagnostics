%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Differences between wave incidence angles, dx=2m
% Simon Treillou, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

dirpath    = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang5_spread30/'
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
maxims(1)=max(V);
Vstd=squeeze(std(mean(nc{'vbar'}(10:end,:,:),2),1));

figure(1);
plot(x(ix150:ix0)',V(ix150:ix0));
hold on
ylim([-0.1,0.5]);
scatter(HR16_fig3_Vobs(:,1),HR16_fig3_Vobs(:,2),'Marker','*');
ylabel('V ($m.s^{-1}$)','Interpreter','latex');
xlabel('$x$ (m)','Interpreter','latex');


dirpath    = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang7_spread30/'
fname      = strcat(dirpath,'rip_avg.nc');
nc=netcdf(fname,'r');
V=squeeze(mean(mean(nc{'vbar'}(10:end,:,:),1),2));
maxims(2)=max(V);
plot(x(ix150:ix0)',V(ix150:ix0));

dirpath    = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang10_spread30/'
fname      = strcat(dirpath,'rip_avg.nc');
nc=netcdf(fname,'r');
V=squeeze(mean(mean(nc{'vbar'}(10:end,:,:),1),2));
maxims(3)=max(V);
plot(x(ix150:ix0)',V(ix150:ix0));
legend('observations','Ang=5','Ang=7','Ang=10', ...
    'Interpreter','latex','Location','northwest');
ylabel('V ($m.s^{-1}$)','Interpreter','latex');
xlabel('$x$ (m)','Interpreter','latex');

figure(2);
plot([5,7,10],maxims,'-o');
ylabel('$V_{max}$ ($m.s^{-1}$)','Interpreter','latex');
xlabel('Incidence angle (°)','Interpreter','latex');
