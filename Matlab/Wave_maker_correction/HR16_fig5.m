%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
%fname     = '../rip_his.nc';  % CROCO file name
%dirpath    = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx2_ang8_spread30_init016_a023_w85_suite/';
dirpath    = '/Users/simon/Code/IB09/IB09_PSbienposT2/';
dirpath    = '/Users/simon/Code/CONFIGS/IB09_randomphase/';   % croco history file name
dirpath    = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/';   % croco history file name
dirpath    = '/Users/simon/Code/CONFIGS/IB09_S102/';   % croco history file name


fname      = strcat(dirpath,['rip_avg' ...
    '.nc']);
%fname     = '../../../IB09_exp_Qnew_wind_periodic/rip_his.nc';
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
%time = time - time(1);
if tend<tstr; tstr=1; end;

h=nc{'h'}(1,:);
beach_l = abs(min(h)/0.02); % 0.02 is beach slope (cf HR16)
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:,:)-xl+beach_l;
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
  zeta=squeeze(nc{'zeta'}(:,:,:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end

%% Calculating times
%time=squeeze(nc{'scrum_time'}(:));
time=squeeze(nc{'scrum_time'}(:));
t=squeeze(nc{'tpas01'}(:,:,:,:));

% t=0:42h
[~,ta] = min(abs(time(:)-(42*60)-time(1)));

% t=1:29h
[~,tb] = min(abs(time(:)-((60+29)*60)-time(1)));

% t=2:03h
[~,tc] = min(abs(time(:)-((2*60+3)*60)-time(1)));

% t=2:22h
[~,td] = min(abs(time(:)-((2*60+22)*60)-time(1)));

% t=4:29h
[~,te] = min(abs(time(:)-((4*60+29)*60)-time(1)));

% t=4:53h
[~,tf] = min(abs(time(:)-((4*60+53)*60)-time(1)));

[~,ixyr] = min(abs(y(:,1)-248));
[~,ix0] = min(abs(x(1,:)));
[~,ixb] = min(abs(x(1,:)+81));
[~,ix310] = min(abs(x(1,:)+310));

%% Calculating downstream maximum position

[~,xb] = min(abs(x(1,:)+81));
[~,xb40] = min(abs(x(1,:)+(81+40)));

ymax_positions=[];
int_times = [ta,tb,tc,td,te,tf];
for j=1:length(int_times)
    tj = squeeze(sum(t(int_times(j),:,:,xb40:xb)));
    imaxtmp = 1;
    for i=1:size(t,3)
        if sum(tj(i,:)>3)>0
            imaxtmp=i;
        end
    end
    ymax_positions(j)=y(imaxtmp,1)-y(jtrac,1);
end
%%
% f = figure;
lim = 15;

subplot(2,3,1);
pA = pcolor(x(1,ix310:ix0),y(:,1)-y(jtrac,1),squeeze(t(ta,end,:,ix310:ix0)));
set(pA, 'EdgeColor', 'none');
ylabel('$y$ (m)','Interpreter','latex');
s = seconds(time(ta)-time(1));
s.Format = 'hh:mm';
title(strcat("(a) t = ",string(s),' h'));
colorbar();
caxis([0 lim]);
hold on
plot(ones(length(y(:,1)),1).*x(1,ixb),y(:,1),'Color','red', ...
    'LineStyle','--');
plot(x(1,ix310:ix0),x(:,ix310:ix0)./x(:,ix310:ix0) .*y(ixyr,1), ...
    'Color','red','LineStyle','--');
colormap("parula");
ylim([0,1900]);
plot(-90,ymax_positions(1),'MarkerSize',11,'MarkerFaceColor','green','Marker','v');
[dyr,x82] = min(abs(y(:,1)-82-y(jtrac,1)));
plot(-20,y(x82)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x248] = min(abs(y(:,1)-248-y(jtrac,1)));
plot(-20,y(x248)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x546] = min(abs(y(:,1)-546-y(jtrac,1)));
plot(-20,y(x546)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1069] = min(abs(y(:,1)-1069-y(jtrac,1)));
plot(-20,y(x1069)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1662] = min(abs(y(:,1)-1662-y(jtrac,1)));
plot(-20,y(x1662)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
plot(x(1,itrac),0,'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','green');

subplot(2,3,2);
pB = pcolor(x(1,ix310:ix0),y(:,1)-y(jtrac,1),squeeze(t(tb,end,:,ix310:ix0)));
set(pB, 'EdgeColor', 'none');
colorbar();
s = seconds(time(tb)-time(1));
s.Format = 'hh:mm';
title(strcat("(b) t = ",string(s),' h'));
caxis([0 lim]);
hold on
plot(ones(length(y(:,1)),1).*x(1,ixb),y(:,1),'Color','red', ...
    'LineStyle','--');
plot(x(1,ix310:ix0),x(:,ix310:ix0)./x(:,ix310:ix0) .*y(ixyr,1), ...
    'Color','red','LineStyle','--');
ylim([0,1900]);
plot(-90,ymax_positions(2),'MarkerSize',11,'MarkerFaceColor','green','Marker','v');
[dyr,x82] = min(abs(y(:,1)-82-y(jtrac,1)));
plot(-20,y(x82)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x248] = min(abs(y(:,1)-248-y(jtrac,1)));
plot(-20,y(x248)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x546] = min(abs(y(:,1)-546-y(jtrac,1)));
plot(-20,y(x546)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1069] = min(abs(y(:,1)-1069-y(jtrac,1)));
plot(-20,y(x1069)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1662] = min(abs(y(:,1)-1662-y(jtrac,1)));
plot(-20,y(x1662)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
plot(x(1,itrac),0,'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','green');

subplot(2,3,3);
pC = pcolor(x(1,ix310:ix0),y(:,1)-y(jtrac,1),squeeze(t(tc,end,:,ix310:ix0)));
set(pC, 'EdgeColor', 'none');
colorbar();
s = seconds(time(tc)-time(1));
s.Format = 'hh:mm';
title(strcat("(c) t = ",string(s),' h'));
caxis([0 lim]);
hold on
plot(ones(length(y(:,1)),1).*x(1,ixb),y(:,1),'Color','red', ...
    'LineStyle','--');
plot(x(1,ix310:ix0),x(:,ix310:ix0)./x(:,ix310:ix0) .*y(ixyr,1), ...
    'Color','red','LineStyle','--');
ylim([0,1900]);
plot(-90,ymax_positions(3),'MarkerSize',11,'MarkerFaceColor','green','Marker','v');
[dyr,x82] = min(abs(y(:,1)-82-y(jtrac,1)));
plot(-20,y(x82)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x248] = min(abs(y(:,1)-248-y(jtrac,1)));
plot(-20,y(x248)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x546] = min(abs(y(:,1)-546-y(jtrac,1)));
plot(-20,y(x546)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1069] = min(abs(y(:,1)-1069-y(jtrac,1)));
plot(-20,y(x1069)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1662] = min(abs(y(:,1)-1662-y(jtrac,1)));
plot(-20,y(x1662)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
plot(x(1,itrac),0,'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','green');


subplot(2,3,4);
pD = pcolor(x(1,ix310:ix0),y(:,1)-y(jtrac,1),squeeze(t(td,end,:,ix310:ix0)));
set(pD, 'EdgeColor', 'none');
xlabel('$x$ (m)','Interpreter','latex');
ylabel('$y$ (m)','Interpreter','latex');
s = seconds(time(td)-time(1));
s.Format = 'hh:mm';
title(strcat("(d) t = ",string(s),' h'));
caxis([0 lim])
colorbar();
hold on
plot(ones(length(y(:,1)),1).*x(1,ixb),y(:,1),'Color','red', ...
    'LineStyle','--');
plot(x(1,ix310:ix0),x(:,ix310:ix0)./x(:,ix310:ix0) .*y(ixyr,1), ...
    'Color','red','LineStyle','--');
ylim([0,1900]);
plot(-90,ymax_positions(4),'MarkerSize',11,'MarkerFaceColor','green','Marker','v');
[dyr,x82] = min(abs(y(:,1)-82-y(jtrac,1)));
plot(-20,y(x82)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x248] = min(abs(y(:,1)-248-y(jtrac,1)));
plot(-20,y(x248)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x546] = min(abs(y(:,1)-546-y(jtrac,1)));
plot(-20,y(x546)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1069] = min(abs(y(:,1)-1069-y(jtrac,1)));
plot(-20,y(x1069)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1662] = min(abs(y(:,1)-1662-y(jtrac,1)));
plot(-20,y(x1662)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
plot(x(1,itrac),0,'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','green');


subplot(2,3,5);
pE = pcolor(x(1,ix310:ix0),y(:,1)-y(jtrac,1),squeeze(t(te,end,:,ix310:ix0)));
set(pE, 'EdgeColor', 'none');
xlabel('$x$ (m)','Interpreter','latex');
s = seconds(time(te)-time(1));
s.Format = 'hh:mm';
title(strcat("(e) t = ",string(s),' h'));
caxis([0 lim])
colorbar();
hold on
plot(ones(length(y(:,1)),1).*x(1,ixb),y(:,1),'Color','red', ...
    'LineStyle','--');
plot(x(1,ix310:ix0),x(:,ix310:ix0)./x(:,ix310:ix0) .*y(ixyr,1), ...
    'Color','red','LineStyle','--');
ylim([0,1900]);
plot(-90,ymax_positions(5),'MarkerSize',11,'MarkerFaceColor','green','Marker','v');
[dyr,x82] = min(abs(y(:,1)-82-y(jtrac,1)));
plot(-20,y(x82)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x248] = min(abs(y(:,1)-248-y(jtrac,1)));
plot(-20,y(x248)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x546] = min(abs(y(:,1)-546-y(jtrac,1)));
plot(-20,y(x546)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1069] = min(abs(y(:,1)-1069-y(jtrac,1)));
plot(-20,y(x1069)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1662] = min(abs(y(:,1)-1662-y(jtrac,1)));
plot(-20,y(x1662)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
plot(x(1,itrac),0,'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','green');


subplot(2,3,6);
pF = pcolor(x(1,ix310:ix0),y(:,1)-y(jtrac,1),squeeze(t(tf,end,:,ix310:ix0)));
ylim([0,1900]);
set(pF, 'EdgeColor', 'none');
colorbar();
caxis([0 lim])
s = seconds(time(tf)-time(1));
s.Format = 'hh:mm';
title(strcat("(f) t = ",string(s),' h'));
xlabel('$x$ (m)','Interpreter','latex');
hold on
plot(-90,ymax_positions(6),'MarkerSize',11,'MarkerFaceColor','green','Marker','v');
plot(ones(length(y(:,1)),1).*x(1,ixb),y(:,1),'Color','red', ...
    'LineStyle','--');
plot(x(1,ix310:ix0),x(:,ix310:ix0)./x(:,ix310:ix0) .*y(ixyr,1), ...
    'Color','red','LineStyle','--');
[dyr,x82] = min(abs(y(:,1)-82-y(jtrac,1)));
plot(-20,y(x82)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x248] = min(abs(y(:,1)-248-y(jtrac,1)));
plot(-20,y(x248)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x546] = min(abs(y(:,1)-546-y(jtrac,1)));
plot(-20,y(x546)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1069] = min(abs(y(:,1)-1069-y(jtrac,1)));
plot(-20,y(x1069)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
[dyr,x1662] = min(abs(y(:,1)-1662-y(jtrac,1)));
plot(-20,y(x1662)-y(jtrac,1),'Marker','o','MarkerSize',11,'MarkerFaceColor','yellow')
plot(x(1,itrac),0,'Marker','pentagram','MarkerSize',15,'MarkerFaceColor','green');
%text(-70,800,'SZ','FontSize',18,'Color','r');
%text(-250,800,'IS','FontSize',18,'Color','r');


map=colormap(redwhiteblue(-0.5,3,50));
map(7,:)=[0.9    0.9  1.0000];
map(6,:)=[0.8    0.8  1.0000];
map(5,:)=[0.7    0.7  1.0000];
map(4,:)=[0.5    0.5  1.0000];
map(3,:)=[0.3    0.3  1.0000];
map(2,:)=[0.2    0.2  1.0000];
map(1,:)=[0.     0.   1.0000];
colormap(map);