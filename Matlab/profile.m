%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Code by Simon Treillou
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
addpath(genpath("Tools"));
set(groot,'defaultAxesTickLabelInterpreter','latex');  
%================== User defined parameters ===========================
%
% --- model params ---
%
dirpath = '/scratch/users/treillou/';
fnames{1}='IB09_2024_2';
fnames{2}='IB09_2024_3';
fnames{3}='IB09_2024_4';
fnames{4}='IB09_2024_5';
fnames{5}='IB09_2024_5_repeat';
fnames{6}='IB09_2024_6';
fnames{7}='IB09_2024_7';
nbfiles=length(fnames);
makepdf   = 0;             % make pdf file
Dcrit=0.2;
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

for i=1:nbfiles
    fn=[dirpath,fnames{i},'/rip_avg.nc'];
    if i==1
        time=ncread(fn,'scrum_time');
        t=ncread(fn,'tpas01');
        t=squeeze(t(:,:,9,:));
    else
        time=cat(1,time,ncread(fn,'scrum_time'));
        tmp=ncread(fn,'tpas01');
        t=cat(3,t,squeeze(tmp(:,:,9,:)));
    end
end

% Grid
h=ncread(fn,'h');h=h(:,1);
xl=ncread(fn,'xl');
x=findx_IB09(fn);
y=ncread(fn,'y_rho');
pm=ncread(fn,'pm');
pn=ncread(fn,'pn');
N=length(ncread(fn,'s_rho'));


%% 

[~,xn] = min(abs(y(1,:)-60));
xinit = y(1,xn);
%t=squeeze(nc{'tpas01'}(:,9,:,:));
[dyr,ix0] = min(abs(x+8));
[dyr,ixmax] = min(abs(x+356));
time=time-time(1); %spin-up

[dyr,x82] = min(abs(y(1,:)-82-xinit));
[~,tinit] = min(abs(time-2300));
[~,tend] = min(abs(time-19000));
usedtime82=(time(tend)-time(tinit))/(19000-4300);
t82 = squeeze(mean(t(ixmax:ix0,x82,tinit:tend)));
std82 = squeeze(std(t(ixmax:ix0,x82,tinit:tend)));
t82h = t82 + std82./length(tinit:tend);
t82b = t82 - std82./length(tinit:tend);

[dyr,x248] = min(abs(y(1,:)-208-xinit));
[~,tinit] = min(abs(time-1000));
[~,tend] = min(abs(time-23000));
usedtime248=(time(tend)-time(tinit))/(23000-1500);
t248 = squeeze(mean(t(ixmax:ix0,x248,tinit:tend)));
std248 = squeeze(std(t(ixmax:ix0,x248,tinit:tend)));
t248h = t248 + std248/length(tinit:tend);
t248b = t248 - std248/length(tinit:tend);

[dyr,x546] = min(abs(y(:,1)-546-xinit));
[~,tinit] = min(abs(time-5000));
[~,tend] = min(abs(time-24500));
usedtime546=(time(tend)-time(tinit))/(24500-5000);
t546 = squeeze(mean(t(ixmax:ix0,x546,tinit:tend)));
std546 = squeeze(std(t(ixmax:ix0,x546,tinit:tend)));
t546h = t546 + std546/length(tinit:tend);
t546b = t546 - std546/length(tinit:tend);

[dyr,x1069] = min(abs(y(1,:)-1069-xinit));
[~,tinit] = min(abs(time-6800));
[~,tend] = min(abs(time-24000));
usedtime1069=(time(tend)-time(tinit))/(24500-6000);
%tinit=1;
tend=length(time);

t2=t(ixmax:ix0,x1069,tinit:tend);
t1069 = squeeze(mean(t2,3));
std1069 = squeeze(std(squeeze(t2)'))';
t1069h = t1069 + std1069/sqrt(length(tinit:tend));
t1069b = t1069 - std1069/sqrt(length(tinit:tend));


%% Loading data from HR15
load("HR16_fig10_D1069mod_2") % with plotdigitalizer
load("HR16_fig10_D1069obs_2")

%% Plotting
figure('Position',[100 100 1000 400],'PaperOrientation','landscape');
ccroco='blue';cObs='black';cFW='red';
alpha=1.0;
c1=[0,114,189,alpha*255]/255;
c2=[145,0,85,alpha*255]/255;
ccroco=c1;
cFW=c2;
wcroco=3;wFW=3;wObs=wFW;

x2 = [(x(ixmax:ix0)'), fliplr(x(ixmax:ix0)')];
inBetween = [(t1069h.*min(h(ixmax:ix0),2.7)./h(ixmax:ix0))', fliplr((t1069b.*min(h(ixmax:ix0),2.7)./h(ixmax:ix0))')];
fill(x2, inBetween,ccroco(1:3),'linestyle','none','FaceAlpha',0.3);
hold on;
plot(x(ixmax:ix0),t1069.*min(h(ixmax:ix0),2.7)./h(ixmax:ix0),'LineWidth',wcroco,'Color',ccroco);
hold on
plot(HR16_fig10_D1069mod_2(:,1),HR16_fig10_D1069mod_2(:,2),'LineWidth',wFW,'Color',cFW);
plot(HR16_fig10_D1069obs_2(:,1),HR16_fig10_D1069obs_2(:,2),'LineWidth',wObs,'Color',cObs,'LineStyle',':');
xline(-81,lineStyle="--",lineWidth=2);
xlim([-300,0]);
legend('','3D CROCO','2D funwaveC (HR16)', 'Data','Location','best','Interpreter','latex','FontSize',17);
ylabel('$C$  (ppb)','Interpreter','latex','FontSize',15);
grid("minor");
set(gca,'linewidth',1.5);
xlabel('$x$  ($m$)','Interpreter','latex','FontSize',15);

if makepdf
    print(gcf, '-dpdf', '-r600', ['./IB09_tracerprofil.pdf']);
end


rmsCROCO=RMSE_differentarrays(HR16_fig10_D1069obs_2(:,1),HR16_fig10_D1069obs_2(:,2),x(1,ixmax:ix0),t1069.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0));
rmsFUNWAVE=RMSE_differentarrays(HR16_fig10_D1069obs_2(:,1),HR16_fig10_D1069obs_2(:,2),HR16_fig10_D1069mod_2(:,1),HR16_fig10_D1069mod_2(:,2));
disp("------------------------------------------------------")
disp("  - RMSE between obs. and CROCO    : "+string(rmsCROCO)+" ppb")
disp("  - RMSE between obs. and funwaveC : "+string(rmsFUNWAVE)+" ppb")


%%
orange=[255, 193, 19]/255;
bleucyan=[0, 213, 235]/255;
bleufonce=[21, 129, 209]/255;
figure('Position',[100 100 1000 400],'PaperOrientation','landscape');
ccroco='blue';cObs='black';cFW='red';
alpha=1.0;
c1=[0,114,189,alpha*255]/255;
c2=[145,0,85,alpha*255]/255;
ccroco=bleucyan;
cFW='r';
cObs='w';
wcroco=5;wFW=4;wObs=wFW;

x2 = [x(1,ixmax:ix0), fliplr(x(1,ixmax:ix0))];
inBetween = [t1069h.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0), fliplr(t1069b.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0))];
fill(x2, inBetween,ccroco(1:3),'linestyle','none','FaceAlpha',0.3);
hold on;
plot(x(1,ixmax:ix0),t1069.*min(h(1,ixmax:ix0),2.7)./h(1,ixmax:ix0),'LineWidth',wcroco,'Color',ccroco);
plot(HR16_fig10_D1069mod_2(:,1),HR16_fig10_D1069mod_2(:,2),'LineWidth',wFW,'Color',cFW);
plot(HR16_fig10_D1069obs_2(:,1),HR16_fig10_D1069obs_2(:,2),'LineWidth',wObs,'Color',cObs,'LineStyle',':');
xline(-81,'LineStyle',"--",'LineWidth',2,'Color',orange);
xlim([-300,0]);
%legend('','3D CROCO','2D funwaveC (HR16)', 'Data','Location','best','Interpreter','latex','FontSize',17);
%ylabel('$C$  (ppb)','Interpreter','latex','FontSize',15);
grid();
set(gca,'linewidth',1.5);
%xlabel('$x$  ($m$)','Interpreter','latex','FontSize',15);
ax = gca;
ax.FontSize = 18;
ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';
box off
set(gcf,'Color',[5/255, 20/255, 23/255])
set(gca,'Color',[5/255, 20/255, 23/255])

%%
% Generate a curve
x = linspace(-10, 10, 100);  % X values
y = sin(x);  % Curve (e.g., sine function)

% Create upper and lower bounds (+1 and -1 around the curve)
y_upper = y + 1;
y_lower = y - 1;

% Plot the original curve
figure;
plot(x, y, 'k', 'LineWidth', 2); % Black curve

hold on;
% Fill the region between y_upper and y_lower
fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot the upper and lower bounds
plot(x, y_upper, 'r--', 'LineWidth', 1.5); % Upper bound (red dashed)
plot(x, y_lower, 'r--', 'LineWidth', 1.5); % Lower bound (red dashed)

hold off;
xlabel('X');
ylabel('Y');
title('Curve with ±1 Region');
grid on;


