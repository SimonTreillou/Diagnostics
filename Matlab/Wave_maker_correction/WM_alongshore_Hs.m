%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Alongshore significant wave height (for offshore, inner-shelf and 
% surfzone sensors)
% Simon Treillou, 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
stname1 = '/Users/simon/Code/CONFIGS/BAKER_G1D/stations.nc';
stname2 = '/Users/simon/Code/CONFIGS/BAKER_G2B_default/stations.nc';

for i=1:2
    if i==1
        stname=stname1;
    elseif i==2
        stname=stname2;
    end
    nc=netcdf(stname,'r');
    xpos=nc{'Xgrid'}(:);
    ypos=nc{'Ygrid'}(:);
    %sta=find(((ypos==150)));
    sta=find(((xpos==40)));
    % x-positions
    xpos=nc{'Xgrid'}(:);
    xpos=xpos(sta);
    % y-positons
    ypos=nc{'Ygrid'}(:);
    ypos=ypos(sta);
    tstr=100;
    tend=length(nc{'scrum_time'}(:));
    zeta1=nc{'zeta'}(tstr:tend,:);
    zeta1=zeta1(:,sta);
    u=squeeze(nc{'u'}(tstr:tend,:,10));
    u=u(:,sta);
    v=squeeze(nc{'v'}(tstr:tend,:,10));
    v=v(:,sta);
    time=nc{'scrum_time'}(tstr:tend);
    hs=4*std(zeta1);

    if i==1
        hs1=hs;
    elseif i==2;
        hs2=hs;
    end
end

%% DATA

trial="G1D-IS";
%trial="G2B-SZ";
load('/Users/simon/Code/BAKER/'+trial+'/press1.txt')
load('/Users/simon/Code/BAKER/'+trial+'/press2.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press3.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press4.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press5.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press6.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press7.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press8.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press9.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press10.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press11.txt');
load('/Users/simon/Code/BAKER/'+trial+'/press12.txt');

Pobs=[press4 press5 press6 press7 press8 press9 press10 press11 press12];
%Pobs=[press11];
%Pobs=[press1 press2 press3 press4 press5 press6 press7 press8];
%Pobs=[press9 press10 press11 press12];

%Pobs=[press6];

Hsobs=[];
Hz = 100;
maxfac = 1.2;
offset = 0.05;
Fs=100; 

for i=1:size(Pobs,2)
    eta=BAKER_press2eta(squeeze(Pobs(:,i)),Fs,offset,maxfac);
    Hsobs(i)=4*std(eta);
end

%% PLOTTING

figure('Position',[500 500 1400 500]);
Hsobs(end)=nan;
plot(ypos/10-15,Hsobs,'o','Color','k','LineWidth',2);
hold on
xlabel('$y (m)$','Interpreter','latex','FontSize',15);
ylabel('$H_s/H_i$ $(m)$','Interpreter','latex','FontSize',15);
%title('Significant wave height','Interpreter','latex','FontSize',15)
grid("minor");
plot(ypos/10-15,hs2,'LineWidth',3,'Color',[111,160,135,0.7*255]/255);
plot(ypos/10-15,hs1,'LineWidth',3,'Color',[192,24,135,0.7*255]/255);
%ylim([0.5 1.5]);
%xlim([70 190])
legend('Data','Default','Corrected', ...
    'Interpreter','latex','FontSize',15,'Location','northwest');
print(gcf, '-dpdf','-bestfit', './Figures/Baker-hs.pdf');


%%

trial="G1D-IS";
%trial="G2B-SZ";
vr='wg';
ttt = BAKER_read_data('/Users/simon/Code/BAKER/'+trial+'/',vr);
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'1.txt')
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'2.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'3.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'4.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'5.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'6.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'7.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'8.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'9.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'10.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'11.txt');
load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'12.txt');

 path_directory='/Users/simon/Code/BAKER/'+trial+'/'; 
 original_files=dir([path_directory+'wg*']); 
 for k=1:length(original_files)
    filename=[path_directory+original_files(k).name];
    disp(filename)
end


WGobs=[wg1 wg2 wg3 wg4 wg5 wg6 wg7 wg8];
Hz=100;
WL=2^5;
OL=2^4;

for i=1:size(WGobs,2)
    %wg4=wg4(50000:end);
    [See,f,Seec] = pwelch(WGobs(:,i),WL*Hz,OL*Hz,[],Hz);
    [~,ixTp]=max(abs(See));
    Tpobs_off(i)=1/f(ixTp);
    Hsobs_off(i)=4*sqrt(trapz(f,See));
end

%% PLOTTING OFFSHORE ARRAY (X~19m)

figure('Position',[500 500 1400 500]);
%Hsobs_off(end)=nan;
plot(ypos(2:end)/10-15,Hsobs_off,'o','Color','k','LineWidth',2);
hold on
xlabel('$y (m)$','Interpreter','latex','FontSize',15);
ylabel('$H_s$ $(m)$','Interpreter','latex','FontSize',15);
%title('Significant wave height','Interpreter','latex','FontSize',15)
grid("minor");
plot(ypos(2:end)/10-15,hs2(2:end),'LineWidth',3,'Color',[111,160,135,0.7*255]/255);
plot(ypos(2:end)/10-15,hs1(2:end),'LineWidth',3,'Color',[192,24,135,0.7*255]/255);
%ylim([0.5 1.5]);
%xlim([70 190])
legend('Data','Default','Corrected', ...
    'Interpreter','latex','FontSize',15,'Location','northwest');
title("Offshore",'Interpreter','latex','FontSize',18);
print(gcf, '-dpdf','-bestfit', './Figures/Baker-hs_offshore.pdf');