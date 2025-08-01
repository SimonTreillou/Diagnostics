%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure 3 from Hally-Rosendahl & Feddersen, 2016
% Significant wave height Hs and alongshore average velocity V verification
% Simon Treillou, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_randomphase/';
dirpath = '/Users/simon/Code/CONFIGS/IB09_21AVR_spinup/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_S102/';   % croco history file name
fname      = strcat(dirpath,'rip_avg.nc');
hname      = strcat(dirpath,'rip_his.nc');
dianame    = strcat(dirpath,'rip_diags_eddy_avg.nc');
stname  = strcat(dirpath,'stations.nc');
ncs=netcdf(stname,'r');
model3D = 1;
%
%======================================================================
%% --------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname,'r');
ncd=netcdf(dianame,'r');
isdiag = true;
if isempty(ncd); isdiag=false; end
nch=netcdf(hname,'r');

tstr=1;
time=squeeze(nc{'scrum_time'}(:));
tend=length(time);
if tend<tstr; tstr=1; end;

% horizontal grid
yr=squeeze(nc{'y_rho'}(:,:));
[~,yindex] = min(abs(yr(:,1)-248));
hr=squeeze(nc{'h'}(yindex,:));
xr=squeeze(nc{'x_rho'}(yindex,:));
xu=0.5*(xr(1:end-1)+xr(2:end));
L=length(hr);
xl=nc{'xl'}(:);
y=nc{'y_rho'}(:,:);
pm=nc{'pm'}(:,:);
pn=nc{'pn'}(:,:);
dx=unique(pm);
dy=unique(pn);

% vertical grid
h=hr;
beach_l=50;
x=xr-xl+beach_l;
xu=xu-xl+beach_l;
h0=hr;
N=length(nc('s_rho'));
theta_s=nc.theta_s(:); 
theta_b=nc.theta_b(:); 
hc=nc.hc(:); 
zeta=squeeze(nc{'zeta'}(tend,yindex,:));
Dcrit=nc{'Dcrit'}(:)*1.0;
zeta(h<Dcrit)=zeta(h<Dcrit)-h(h<Dcrit); % add land topo
zr=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2));
zw=squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,'w',2));
zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
dz1=zr(1,:)-zw(1,:);
dzu1=zru(1,:)-zwu(1,:);
dzu3=zru(3,:)-zwu(1,:);
%
xr2d=repmat(xr,[N 1]);
xu2d=repmat(xu,[N 1]);
xw2d=repmat(xr,[N+1 1]);
D   =zw(N+1,:)-zw(1,:);
D2d =repmat(D,[N 1]);
Du  =zwu(N+1,:)-zwu(1,:);
Du2d=repmat(Du,[N 1]);
%
[~,ix150] = min(abs(x+150));
[~,ix0] = min(abs(x));
[dyr,ixb] = min(abs(x+81));

%%
% "Loading" positions
HR16_fig3_Vobs = [-135.4128    0.1629;
                  -100.4587    0.1389;
                   -80.9174    0.1217;
                   -59.1743    0.1937;
                   -34.9541    0.2349];
HR16_fig3_Hsobs = [-136.2132    0.7438;
                   -100.9191    0.8277;
                    -81.0662    0.8213;
                    -59.8346    0.7636];

% ---------------------------------------------------------------------
% --- read/compute 3D model fields (tindex) ---
% --------------------------------------------------------------------
% ... zonal velocity ...                         ---> xu,zu
u=squeeze(nc{'u'}(tend,:,yindex,:));

% ... vertical velocity ...                      ---> xr,zw
w=squeeze(nc{'w'}(tend,:,yindex,:));

% ... total viscosity
Akv=squeeze(nc{'AKv'}(tend,:,yindex,:));


if isdiag
    % ---------------------------------------------------------------------
    % --- read/compute significant wave height ---
    % --------------------------------------------------------------------
    % ... wave setup ...  
    % time-averaged and alongshore <zeta>^2 
    t0 = 10;
    sup=squeeze(mean(mean(nc{'zeta'}(t0:end,:,:),1),2))'; %  time-averaged
    z0=sup;
    z0(h0<Dcrit)=z0(h0<Dcrit)-h0(h0<Dcrit); % add slope
    

    % ... Hrms ...
    % time-averaged and alongshore <zeta^2>
    zz=squeeze(mean(ncd{'zz'}(t0:end,:,:),1)); % time-averaged
    zz=squeeze(mean(zz,1)); % average alongshore
    % see Svendsen (2005) p.126
    hrms=4.0083*sqrt(zz-z0.^2)/sqrt(2.0011); %max(0,...) ou abs
    hs=squeeze(sqrt(2.0011)*hrms);
    hs=squeeze(hs);
end

%%
if model3D,
  ui=squeeze(nch{'u'}(1,N,:,:));
  vi=squeeze(nch{'v'}(1,N,:,:));
  u=squeeze(nch{'u'}(tstr:tend,N,:,:));
  v=squeeze(nch{'v'}(tstr:tend,N,:,:));
  ub=squeeze(nch{'u'}(tstr:tend,1,:,:));
  vb=squeeze(nch{'v'}(tstr:tend,1,:,:));
  ubar=squeeze(nch{'ubar'}(tstr:tend,:,:));
  vbar=squeeze(nch{'vbar'}(tstr:tend,:,:));
  t=squeeze(nch{'tpas01'}(:,N,:,:));
  zeta=squeeze(nch{'zeta'}(:,:,:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end

%% Calculating time-averaged alongshore current V at yr=248m
V=squeeze(mean(mean(nc{'vbar'}(5:end,:,:),1),2));
Vstd=squeeze(std(mean(nc{'vbar'}(5:end,:,:),2),1));

%% RMSE: Calculating RMSE between HR16 data points and CROCO simulation
% For Hs
if isdiag
    for i=1:length(HR16_fig3_Hsobs(:,1))
        [dyr,indexData] = min(abs(x-HR16_fig3_Hsobs(i,1)));
        estimate(i) = abs(hs(indexData));
    end
    rmseHs = rmse(HR16_fig3_Hsobs(:,2),estimate');
end
rmseHsHR16=0.03;
% For V
estimate = [];
for i=1:length(HR16_fig3_Vobs(:,1))
    [dyr,indexData] = min(abs(x-HR16_fig3_Vobs(i,1)));
    estimate(i) = abs(V(indexData));
end
rmseV = rmse(HR16_fig3_Vobs(:,2),estimate');
rmseHR15 = 0.02;

%% Plotting
% Plotting Hs
fig=subplot(2,1,1);
ylim([0,1.0]);
hold on;
scatter(HR16_fig3_Hsobs(:,1),HR16_fig3_Hsobs(:,2),'Marker','*');
if isdiag
    plot(x,abs(hs));
    str={'RMSE = ',rmseHs};
end
plot(ncs{'Xgrid'}(7:43)-300,squeeze(4*std(ncs{'zeta'}(100:end,7:43))),'o');
legend('observations', ...
    'CROCO model RMSE='+string(round(rmseHs,3)),'Interpreter', ...
    'latex','Location','southwest');
ylabel('$H_s$ (m)','Interpreter','latex');
xlim([-150, 0])


% Plotting alongshore drift
subplot(2,1,2);
fill([x(ix150:ix0)';flipud(x(ix150:ix0)')],[V(ix150:ix0)-Vstd(ix150:ix0); ...
    flipud(V(ix150:ix0)+Vstd(ix150:ix0))],[.9 .9 .9],'linestyle','none');
hold on;
line(x',V);
ylim([-0.0,0.3]);
scatter(HR16_fig3_Vobs(:,1),HR16_fig3_Vobs(:,2),'Marker','*');
plot(ncs{'Xgrid'}(7:43)-300,squeeze(mean(ncs{'v'}(100:end,7:43,1))),'o');
ylabel('V ($m.s^{-1}$)','Interpreter','latex');
str='RMSE='+string(round(rmseV,3));
legend('CROCO '+str,'observations', ...
    'Interpreter','latex','Location','best');
ylabel('V ($m.s^{-1}$)','Interpreter','latex');
xlim([-150, 0])