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
fnames3D{1}='IB09_2024_2';
fnames3D{2}='IB09_2024_3';
fnames3D{3}='IB09_2024_4';
nbfiles=length(fnames3D);
makepdf   = 0;             % make pdf file
Dcrit=0.2;

%% funwaveC
xHR16 = [0,922.6594301221168,2985.0746268656717,4993.215739484396,5970.149253731343,7951.153324287653,9036.635006784261,9986.431478968792,10990.502035278156,11967.435549525102];
yHR16 = [0,174,774,1194,1146,1380,1548,1638,1860,1956];

%% COMPUTE

for i=1:nbfiles
    fn=[dirpath,fnames3D{i},'/rip_avg.nc'];
    if i==1
        time=ncread(fn,'scrum_time');
        t=ncread(fn,'tpas01');
        t=squeeze(sum(t(:,:,:,:),3));
    else
        time=cat(1,time,ncread(fn,'scrum_time'));
        tmp=ncread(fn,'tpas01');
        t=cat(3,t,squeeze(sum(tmp(:,:,:,:),3)));
    end
end
time=time-time(1);

% Grid
h=ncread(fn,'h');h=h(:,1);
xl=ncread(fn,'xl');
x=findx_IB09(fn);
y=ncread(fn,'y_rho');
pm=ncread(fn,'pm');
pn=ncread(fn,'pn');
N=length(ncread(fn,'s_rho'));

tracinit=60;
y=y(:,tracinit:end)-y(1,tracinit);
indx=1:int16(length(time)/100):length(time);
[~,x0] = min(abs(x(:,1)+81-40));
[~,xb] = min(abs(x(:,1)+81+40));
[~,xb40] = min(abs(x(:,1)+(81+40)));
%t=squeeze(nc{'tpas01'}(indx,:,tracinit:end,xb:x0));
ymax = max(y(1,:));
tinit = max(max(t(:,:,1)));
load("HR16_fig6.mat");

%
t2 = t(xb:x0,tracinit:end,indx);
ymax_positions=[];
for j=1:size(t2,3)
    tj = squeeze(t2(:,1:end-10,j));
    imaxtmp = 1;
    for i=1:size(tj,2)
        if sum(tj(:,i)>3)>0
            imaxtmp=i;
        end
    end
    ymax_positions(j)=y(1,imaxtmp);
end

time3D=time(indx);
ymax_positions3D=ymax_positions;

%% Plotting
figure('Position',[100 100 700 300],'PaperOrientation','landscape');
ccroco='blue';cObs='black';cFW='red';
alpha=1.0;
c1=[0,114,189,alpha*255]/255;
c2=[145,0,85,alpha*255]/255;
ccroco=c1;
cFW=c2;
wcroco=3;wFW=3;wObs=wFW;

%plot(time2D,ymax_positions2D,'LineWidth',3,'Color',cFW);
plot(xHR16,yHR16,'LineWidth',3,'Color',cFW);
hold on
plot(time3D,ymax_positions3D,'LineWidth',3,'Color',ccroco);
plot(HR16_fig6(:,1),HR16_fig6(:,2),'ok','MarkerSize',7);
plot(0:2e4,ones(length(0:2e4),1).*ymax,'LineStyle','--','LineWidth',1.5,'Color','k');
ylabel('$y_p$ (m)','Interpreter','latex','FontSize',15);
legend('2D funwaveC','3D CROCO','Obs.','Location','southeast','FontSize',12,'Interpreter','latex');
%plot(time,0.1936*time,'LineStyle','--','Color','red')
xlabel('$t$ (s)','Interpreter','latex','FontSize',15);
txt = '$y_{max}$';
ylim([0 1800])
xlim([0 10000]);
text(1000,ymax-50,txt,'Interpreter','latex','FontSize',15);
grid("minor");
set(gca,'linewidth',1);

if makepdf
    print(gcf, '-dpdf', '-r600', ['./Figures/IB09_tracerplumeedge.pdf']);
end