%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
%fname     = '../rip_his.nc';  % CROCO file name
dirpath    = '/Users/simon/Code/IB09/IB09_PSbienposT2/';
%dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
%dirpath    = '/Users/simon/Code/Matlab/';
%dirpath    = '/Users/simon/Code/CONFIGS/IB09/';
%dirpath    = '/Users/simon/Code/CONFIGS/Compile_Patrick/'
fname      = strcat(dirpath,'rip_his.nc');
%fname     = '../../../IB09_exp_Qnew_wind/rip_his.nc';
fname     = '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/rip_histot.nc';   % croco history file name
fname     = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_histot.nc';   % croco history file name

model3D   = 1;

makepdf   = 0;             % make pdf file
Dcrit=0.2;
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);

dx=3;

tstr=10;
time=nc{'scrum_time'}(:);
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1; end;
net_v=2.13e6;

h=nc{'h'}(1,:);
beach_l = 50;%abs(min(h)/0.02); % 0.02 is beach slope (cf HR16)
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:,:)-xl+beach_l;
y=nc{'y_rho'}(:,:);
pm=nc{'pm'}(:,:);
pn=nc{'pn'}(:,:);
N=length(nc('s_rho'));

if model3D,
  ui=squeeze(nc{'u'}(1,N,:,:));
  vi=squeeze(nc{'v'}(1,N,:,:));
  u=squeeze(nc{'u'}(tstr:tend,N,:,:));
  v=squeeze(nc{'v'}(tstr:tend,N,:,:));
  ub=squeeze(nc{'u'}(tstr:tend,1,:,:));
  vb=squeeze(nc{'v'}(tstr:tend,1,:,:));
  ubar=squeeze(nc{'ubar'}(tstr:tend,:,:));
  vbar=squeeze(nc{'vbar'}(tstr:tend,:,:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end


%% 

[~,xn] = min(abs(y(:,1)-350));
xinit = y(xn,1);
t=squeeze(nc{'tpas01'}(:,10,:,:));
[dyr,ix0] = min(abs(x(1,:)+8));
[dyr,ixmax] = min(abs(x(1,:)+356));
time=time-time(1);
decal=900;

[dyr,x82] = min(abs(y(:,1)-82-xinit));
tstart=4300;
tstop=19000;
[~,tinit] = min(abs(time-tstart));
[~,tend1] = min(abs(time-tstart-1*decal));
[~,tend2] = min(abs(time-tstart-2*decal));
[~,tend3] = min(abs(time-tstart-3*decal));
[~,tend4] = min(abs(time-tstart-4*decal));
[~,tend] = min(abs(time-tstop));
t821 = squeeze(mean(t(tinit:tend1,x82,ixmax:ix0)));
t822 = squeeze(mean(t(tinit:tend2,x82,ixmax:ix0)));
t823 = squeeze(mean(t(tinit:tend3,x82,ixmax:ix0)));
t824 = squeeze(mean(t(tinit:tend4,x82,ixmax:ix0)));
t82 = squeeze(mean(mean(t(tinit:tend,x82-20:x82+20,ixmax:ix0),2)));

[dyr,x248] = min(abs(y(:,1)-248-xinit));
tstart=1500;
tstop=23000;
tstart=6000;
tstop=24500;
[~,tinit] = min(abs(time-tstart));
[~,tend1] = min(abs(time-tstart-1*decal));
[~,tend2] = min(abs(time-tstart-2*decal));
[~,tend3] = min(abs(time-tstart-3*decal));
[~,tend4] = min(abs(time-tstart-4*decal));
[~,tend] = min(abs(time-tstop));
t2481 = squeeze(mean(t(tinit:tend1,x248,ixmax:ix0)));
t2482 = squeeze(mean(t(tinit:tend2,x248,ixmax:ix0)));
t2483 = squeeze(mean(t(tinit:tend3,x248,ixmax:ix0)));
t2484 = squeeze(mean(t(tinit:tend4,x248,ixmax:ix0)));
t248 = squeeze(mean(mean(t(tinit:tend,x248-20:x248+20,ixmax:ix0),2)));

[dyr,x1069] = min(abs(y(:,1)-1069-xinit));
tstart=6000;
tstop=24500;
[~,tinit] = min(abs(time-tstart));
[~,tend1] = min(abs(time-tstart-1*decal));
[~,tend2] = min(abs(time-tstart-2*decal));
[~,tend3] = min(abs(time-tstart-3*decal));
[~,tend4] = min(abs(time-tstart-4*decal));
[~,tend] = min(abs(time-tstop));
t10691= squeeze(mean(t(tinit:tend1,x1069,ixmax:ix0)));
t10692= squeeze(mean(t(tinit:tend2,x1069,ixmax:ix0)));
t10693= squeeze(mean(t(tinit:tend3,x1069,ixmax:ix0)));
t10694= squeeze(mean(t(tinit:tend4,x1069,ixmax:ix0)));
t1069 = squeeze(mean(mean(t(tinit:tend,x1069-30:x1069+30,ixmax:ix0),2)));

%% Loading data from HR15
load Data/HR16_fig10_D82mod.mat
load Data/HR16_fig10_D82obs.mat
load Data/HR16_fig10_D248mod.mat
load Data/HR16_fig10_D248obs.mat
load Data/HR16_fig10_D1069mod.mat
load Data/HR16_fig10_D1069obs.mat

%% Plotting
figure;

fig=subplot(3,1,1);
plot(x(1,ixmax:ix0),t82);
hold on
scatter(H16_fig10_D82obs(:,1),H16_fig10_D82obs(:,2));
plot(x(1,ixmax:ix0),t821);
plot(x(1,ixmax:ix0),t822);
plot(x(1,ixmax:ix0),t823);
plot(x(1,ixmax:ix0),t824);
xline(-81,lineStyle="--",lineWidth=2);
legend('Full time','obs','1h','2h','3h','4h', ...
    'Location','northwest','Interpreter','latex');
ylabel('$\langle D \rangle$  (ppb)','Interpreter','latex');
title(strcat("(a) y = ",string(82),' m'),'Interpreter','latex');

subplot(3,1,2);
plot(x(1,ixmax:ix0),t248);
hold on
%plot(HR16_fig10_D248mod(:,1),HR16_fig10_D248mod(:,2));
scatter(HR16_fig10_D248obs(:,1),HR16_fig10_D248obs(:,2));
plot(x(1,ixmax:ix0),t2481);
plot(x(1,ixmax:ix0),t2482);
plot(x(1,ixmax:ix0),t2483);
plot(x(1,ixmax:ix0),t2484);
xline(-81,lineStyle="--",lineWidth=2);
legend('Full time','obs','1h','2h','3h','4h', ...
    'Location','northwest','Interpreter','latex');
ylabel('$\langle D \rangle$ (ppb)','Interpreter','latex');
title(strcat("(b) y = ",string(248),' m'),'Interpreter','latex');

subplot(3,1,3);
[d,ix] = min(abs(y(:,1)-1069));
plot(x(1,ixmax:ix0),t1069);
hold on
xlim([-350,0]);
%plot(HR16_fig10_D1069mod(:,1),HR16_fig10_D1069mod(:,2));
scatter(HR16_fig10_D1069obs(:,1),HR16_fig10_D1069obs(:,2));
plot(x(1,ixmax:ix0),t10691);
plot(x(1,ixmax:ix0),t10692);
plot(x(1,ixmax:ix0),t10693);
plot(x(1,ixmax:ix0),t10694);
xline(-81,lineStyle="--",lineWidth=2);
legend('Full time','obs','1h','2h','3h','4h', ...
    'Location','northwest','Interpreter','latex');
ylab = ylabel(' D  (ppb)','Interpreter','latex');
title(strcat("(c) y = ",string(1069),' m'),'Interpreter','latex');

