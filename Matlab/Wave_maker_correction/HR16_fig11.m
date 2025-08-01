%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
%fname     = '/Users/simon/rip_avg.nc';  % CROCO file name
%fname     = '/Users/simon/Code/CONFIGS/IB09_exp_Qnew/rip_his.nc';
dirpath    = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/';   % croco history file name
dirpath    = '/Users/simon/Code/CONFIGS/IB09_S102/';   % croco history file name
%dirpath    = '/Users/simon/Code/CONFIGS/IB09/';
%dirpath    = '/Users/simon/Code/CONFIGS/Compile_Patrick/'
fname      = strcat(dirpath,'rip_his.nc');
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

tstr=1;
tend=length(nc{'scrum_time'}(:));
time=squeeze(nc{'scrum_time'}(:));
time=time-time(1);
%time=time-time(1);
if tend<tstr; tstr=1; end;
net_v = 8.9e5; 

xl=nc{'xl'}(:);
h=nc{'h'}(1,:);
beach_l = abs(min(h)/0.02);
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
  t=squeeze(nc{'tpas01'}(:,N,:,:));
  zeta=squeeze(nc{'zeta'}(:,:,:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end

%% Modeled SZ+IS cumulative (time-integrated) alongshore dye transport at yr=248
[dyr,ixyr] = min(abs(y(:,1)-248-y(jtrac,1)));
[dyr,x0] = min(abs(x(1,:)-0));
[dyr,x140] = min(abs(x(1,:)+140));

h=nc{'h'}(:,x140:x0);
t=squeeze(sum(nc{'tpas01'}(:,:,:,:),2));

level=1;
V=squeeze(nc{'v'}(:,level,ixyr,x140:x0)); % Alongshore current V(t,x) at yr
D=squeeze(nc{'tpas01'}(:,level,ixyr,x140:x0)); % Dye concentration D(t,x) at yr
d=zeros(tend,length(y(:,1)),length(x(1,x140:x0)));
for i=1:tend
    d(i,:,:)=h+squeeze(zeta(i,:,x140:x0)); % total water depth
end
d=squeeze(mean(d(:,:,:),2));

integ=d.*D.*V;
tau=squeeze(trapz(x(1,x140:x0),integ,2));


cumul_tau=zeros(1,tend);

for i=2:tend
    cumul_tau(i)=trapz(time(1:i),tau(1:i));
end

%% Modeled SZ cumulative (time-integrated) alongshore dye transport at yr=248
[dyr,xb] = min(abs(x(1,:)+82));

VSZ=squeeze(nc{'v'}(:,level,ixyr,xb:x0)); % Alongshore current V(t,x) at yr
DSZ=squeeze(nc{'tpas01'}(:,level,ixyr,xb:x0)); % Dye concentration D(t,x) at yr
dSZ=zeros(tend,length(y(:,1)),length(xb:x0));
h=nc{'h'}(:,xb:x0);
for i=1:tend
    dSZ(i,:,:)=h+squeeze(zeta(i,:,xb:x0)); % total water depth
end
dSZ=squeeze(mean(dSZ(:,:,:),2));

integSZ=dSZ.*DSZ.*VSZ;
tauSZ=squeeze(trapz(x(1,xb:x0),integSZ,2));
cumul_tauSZ=zeros(1,tend);

for i=2:tend
    cumul_tauSZ(i)=trapz(time(1:i),tauSZ(1:i));
end

%% Modeled IS cumulative (time-integrated) alongshore dye transport at yr=248
VIS=squeeze(nc{'v'}(:,level,ixyr,x140:xb)); % Alongshore current V(t,x) at yr
DIS=squeeze(nc{'tpas01'}(:,level,ixyr,x140:xb)); % Dye concentration D(t,x) at yr
dIS=zeros(tend,length(y(:,1)),length(x140:xb));
h=nc{'h'}(:,x140:xb);
for i=1:tend
    dIS(i,:,:)=h+squeeze(zeta(i,:,x140:xb)); % total water depth
end
dIS=squeeze(mean(dIS(:,:,:),2));

integIS=dIS.*DIS.*VIS;
tauIS=squeeze(trapz(x(1,x140:xb),integIS,2));
cumul_tauIS=zeros(1,tend);

for i=2:tend
    cumul_tauIS(i)=trapz(time(1:i),tauIS(1:i));
end


%% Plotting
% Load collected data from HR16
load Data/HR16_fig11_modSZIS.mat;
Q=0.058
HR16_fig11_modSZIS(:,1)=sort(HR16_fig11_modSZIS(:,1));
HR16_fig11_modSZIS(:,2)=sort(HR16_fig11_modSZIS(:,2));
fig=plot(HR16_fig11_modSZIS(:,1)*1e4,HR16_fig11_modSZIS(:,2)*1e6);
%fig=plot(HR16_fig11_modSZIS(:,1),HR16_fig11_modSZIS(:,2)*1e6 /(512*20000));
hold on;
%plot(time./1e4,time*512);
%plot(time./1e4,cumul_tau /(max(time)*512));
plot(time,cumul_tau);
plot(time,cumul_tauSZ);
plot(time,cumul_tauIS);
legend('funwaveC modeled IS+SZ','CROCO modeled IS+SZ','CROCO modeled SZ', ...
    'CROCO modeled IS','Interpreter','latex','Location','northwest');
xlabel('$t$ (s)','Interpreter','latex');
ylabel('$(\int_0^t T^y d\tau)/QT_f$ (ppb m3)','Interpreter','latex');
xlim([0,2e4]);
ylim([0,10e6]);
grid("minor");
%saveas(fig,'HR16_fig11.pdf');

%%
[dyr,a] = min(abs(x(1,:)+182));
% %%
% dirpath    = '/Users/simon/';
% fname = strcat(dirpath,'rip_his_822206.nc');
% nc=netcdf(fname);
% t=nc{'tpas01'}(:);
% time=nc{'scrum_time'}(:);%-50000;
% zeta=nc{'zeta'}(:);
% h=nc{'h'}(:);
% vbar=nc{'vbar'}(:);
% v=nc{'v'}(:);
% yr=nc{'y_rho'}(:);
% xr=nc{'x_rho'}(:);
% y=nc{'y_rho'}(:,1);
% x=nc{'x_rho'}(1,:);
% xl=nc{'xl'}(:);
% pm=nc{'pm'}(:);
% pn=nc{'pn'}(:);
% net_v = 8.9e5; 
% Q=512; %ppb m3 s-1
% t0=0;
% dt=0.05;
% dx=1/unique(pm);
% dy=1/unique(pn);
% xT = -10;
% yT = 67.5;
% beach_l = abs(min(h(1,:))/0.02); % 0.02 is beach slope (cf HR16)
% x = nc{'x_rho'}(1,:)-xl+beach_l;
% [~,indexy] = min(abs(y-yT));
% [~,indexx] = min(abs(x-xT));
% hT = h(indexy,indexx);
% 
% [~,ix]=min(abs(time-t0));
% Hz = (squeeze(zeta(:,indexy,indexx)) + hT);
% cumul = (time-t0)*Q;
% cumul(1:ix) = 0.;
% cumul(ix+1:end) =  (time(ix+1:end)-t0).*Q;% .* net_v;% ./ (dx*dy.*Hz(ix+1:end));
% 
% 
% cumul_CROCO = time.*0;
% lost=0.;
% for i=1:length(time)
%      %lost=sum(sum(t(i-1,:,1,:)))+lost*0.98;
%     cumul_CROCO(i)=sum(sum(sum(t(i,:,:,:))))/20.;% *net_v * 1e-9;%+lost;
% end
% 
% cumul_outbounds = time.*0;
% perdu=0.;
% cumul_perdu=time.*0;
% for i=2:(length(time)-1)
%     perdu=0.;
%     v=vbar(i,2,:);
%     for j=1:length(v)
%         perdu=perdu-min(0,sign(v(j)))*sum(t(i,:,2,j));
%         perdu=perdu-min(0,sign(vbar(i,1,j)))*sum(t(i,:,1,j));
%         perdu=perdu+max(0,sign(vbar(i,end,j)))*sum(t(i,:,end,j));
%         perdu=perdu+max(0,sign(vbar(i,end-1,j)))*sum(t(i,:,end-1,j));
%     end
%     cumul_perdu(i)=perdu+cumul_perdu(i-1);
%     cumul_outbounds(i)=sum(sum((t(i,:,1,:))))+cumul_outbounds(i-1);
% end
% 
