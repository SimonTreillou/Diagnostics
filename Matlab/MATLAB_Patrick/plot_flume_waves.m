close all
clear all
%
%======================================================================
%  COMPUTE LINEAR SUMMATION OF WAVE SPECTRUM
%======================================================================
%
% Vyzikas et al., 2018: The evolution of free and bound waves during 
% dispersive focusing in a numerical and physical flume. Coastal 
% Engineering, 132, 95–109.
%
scale=0.154/0.05;  % scaling for flume experiment
h=1;               % Tank depth
g=9.8;

fname='flume_his_A00.nc';
%======================================================================
%
%
% Read linear waves     -------------
%
load datwaves.txt
amp=datwaves(:,1);
frq=datwaves(:,2);
pha=datwaves(:,3);
kw =datwaves(:,4);
%
t0=64;        % focal time
x0=14.1;      % and position
tl=[0:0.02:128];
x=x0;         % waves at focal point
for i=1:length(tl);
  zl(i)=scale*sum(amp.*cos(kw*(x-x0) -frq*(tl(i)-t0) -pha));
end
maxzl=max(zl);
[itl]=find(zl==maxzl);
maxtl=tl(itl);
tl=tl-maxtl;

%
% Read CROCO waves     ------------------
%
nc=netcdf(fname);
z=squeeze(nc{'zeta'}(:,2,:));
x=squeeze(nc{'x_rho'}(2,:));
t=nc{'scrum_time'}(:);
close(nc)
[LT LX]=size(z);
maxz=max(max(z));
[it ix]=find(z==maxz);
maxt=t(it);
maxx=x(ix);
maxz=z(it,ix);
t=t-maxt;
disp(['max level ',num2str(maxz),' at time ',num2str(maxt), ...
                           ' and position ',num2str(maxx)]);
%
%  Flume experiment
%
load Data_flume
tf=Data_flume(:,1);
zf=Data_flume(:,2);
zfi=interp1(tf,zf,t,'cubic');
zfi=smooth(zfi,7);
%
%  MAKE PLOT
%
figure
hold on
plot(t,zfi,'color',[0.7 0.7 0.7],'linewidth',3);
plot(t,z(:,ix),'color','k','linewidth',2);
plot(tl,zl,'color','k','linewidth',1,'linestyle','--');
legend('Experiment','CROCO','Linear theory')
grid on
axis([-2 +2 -0.15 0.25])
xlabel('Time [sec]')
ylabel('Water level [m]')
title('Rogue wave')
set(gca,'fontsize',15)
hold off

export_fig -transparent rogue_CROCO.pdf

return

for k=1:LT
 zo(k)=interp1(x,z(k,:),14.1);
end
figure
plot(t-64,zo);
grid on
axis([-2 +2 -0.15 0.25])
xlabel('Time [sec]')
ylabel('Water level [m]')
title('water level at x=14.1')


