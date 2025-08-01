%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
%fname     = '../rip_his.nc';  % CROCO file name
%fname     = ['/Users/simon/Code/IB09/IB09_fini/rip_his.nc'];
%dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
dirpath    = '/Users/simon/Code/IB09/IB09_PSbienposT2/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_randomphase/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/';   % croco history file name
fname      = strcat(dirpath,'rip_histot.nc');
model3D   = 1;
itrac=290;
jtrac=60;
makepdf   = 0;             % make pdf file
Dcrit=0.2;
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);

dx=3;

tstr=30;
tend=length(nc{'scrum_time'}(:));
time=squeeze(nc{'scrum_time'}(:));
if tend<tstr; tstr=1; end;

h=nc{'h'}(:,1:252);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:,:)-xl+50;
y=nc{'y_rho'}(:,:);
pm=nc{'pm'}(:,:);
pn=nc{'pn'}(:,:);
N=length(nc('s_rho'));
x0=252;

if model3D,
  ui=squeeze(nc{'u'}(1,N,:,:));
  vi=squeeze(nc{'v'}(1,N,:,:));
  u=squeeze(nc{'u'}(tstr:tend,N,:,:));
  v=squeeze(nc{'v'}(tstr:tend,N,:,:));
  ub=squeeze(nc{'u'}(tstr:tend,1,:,:));
  vb=squeeze(nc{'v'}(tstr:tend,1,:,:));
  ubar=squeeze(nc{'ubar'}(tstr:tend,:,:));
  vbar=squeeze(nc{'vbar'}(tstr:tend,:,:));
  t=squeeze(nc{'tpas01'}(:,N,:,:));
  zeta=squeeze(nc{'zeta'}(:,:,:));
  time=squeeze(nc{'scrum_time'}(:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end

%%

[dyr,y82] = min(abs(y(:,1)-82-y(jtrac,1)));
[dyr,y248] = min(abs(y(:,1)-248-y(jtrac,1)));
[dyr,y546] = min(abs(y(:,1)-546-y(jtrac,1)));
[dyr,ix0] = min(abs(x(1,:)+40));

%%

subplot(3,1,1);
plot(time./1e4-0.2, t(:,y82,ix0));
ylabel('$D$ (ppb)','Interpreter','latex');
%ylim([0,300])
%xticklabels({'0','0.2','0.4','0.6','0.8','1','1.2','1.4','1.6','1.8 \times 10^4'});
title("(a) $y =$"+string(y(y82,1))+" m",'Interpreter','latex');

subplot(3,1,2);
plot(time./1e4-0.2,t(:,y248,ix0));
ylabel('$D$ (ppb)','Interpreter','latex');
%ylim([0,200])
title("(b) $y =$"+string(y(y248,1))+" m",'Interpreter','latex');


subplot(3,1,3);
plot(time./1e4-0.2,t(:,y546,ix0));
%ylim([0,60])
xlabel('$t$ (s)','Interpreter','latex');
title("(c) $y =$"+string(y(y546,1))+" m",'Interpreter','latex');
ylabel('$D$ (ppb)','Interpreter','latex');