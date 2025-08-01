clear all
close all
% ------ DATA -------
% xZ/Z: topo
% xB/B: breaking
% xV/V: longshore current
load('DATA/GPP2014_topo_currents_breaking.mat')
x0=75;
xZ=-xZ+x0;
Z=Z-0.8;

fname     = 'rip_his_FUNC02.nc';
nc=netcdf(fname);
xl=nc{'xl'}(:);
xr=squeeze(nc{'x_rho'}(1,:));
yr=squeeze(nc{'y_rho'}(:,1));
close(nc);

Length_XI=xl;
xs=85;        % inner surf zone 
db=50;        % distance from xs to sand bar
xx=Length_XI-xr;

alpha=0.025;
h1=-4.5  ...
      -1.7*exp(-6*(((xx-xs-db)./db).^2)) ...
      +3.1*(1+tanh(0.025*(xx-xs))) ...
      +0.014*(xx+log(cosh(alpha*(xx-xs))./cosh(alpha*xs))./alpha);

alpha=0.025;
h2=-4.5  ...
      -1.4*exp(-6*(((xx-xs-db)./db).^2)) ...
      +2.7*(1+tanh(0.028*(xx-xs))) ...
      +0.0185*(xx+log(cosh(alpha*(xx-xs))./cosh(alpha*xs))./alpha);

h3=11.3-0.88e-5*(abs(xx-700)).^2.17;

h4=min(h2,h3);

xr=xr-xl+93;

figure('position',[1500 500 500 500])
%plot(xr,-h,'r'); hold on
%plot(xr,-h2,'b'); hold on
plot(xr,-h4,'k'); hold on
plot(xZ,Z,'r--','linewidth',2);
axis([-600 +100 -14 3])
grid on

return
yy=yr;
yper=0;
eps=0.01;
for iper=1:3;
  yper = yper + eps*cos(2*pi*iper*yy./(xs+db) ...
                       +2*pi*tanh(iper/10));
end
figure
plot(yr,yper)

