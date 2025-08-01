%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Alongshore significant wave height (for inner-shelf sensors)
% Simon Treillou, 2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
stname1 = '/Users/simon/Code/CONFIGS/BAKER_G1D/stations.nc';
stname2 = '/Users/simon/Code/CONFIGS/BAKER_G2B_default/stations.nc';

tstr = 15*60;  % 15 minutes spin up


fname ='/Users/simon/Code/CONFIGS/BAKER_G1D/rip_his.nc';
nc=netcdf(fname);
h=nc{'h'}(1,:);
x=nc{'x_rho'}(1,:);

%% DATA 

repo="/Users/simon/Code/BAKER/G1D-IS/";
vr='press';
PRESS = BAKER_read_data(repo,vr);

ini=5;
clear PRESSobs Xobs Yobs Zobs
for k=ini:11
    Xobs(k-ini+1) = PRESS.("press"+string(k)).x;
    Yobs(k-ini+1) = PRESS.("press"+string(k)).y;
    Zobs(k-ini+1) = PRESS.("press"+string(k)).z;
    PRESSobs(1:size(PRESS.press1.data),k-ini+1) = PRESS.("press"+string(k)).data;
end
Xobs_IS2 = mean(Xobs);
Zobs_IS2 = mean(Zobs);

Hz = 100;
maxfac = 1.2;
[~,ix]=min(abs(x+15-Xobs_IS2));
offset= Zobs_IS2 - (1.07-h(ix));
Fs=100; 
WL=32;
OL=16;

clear Hsobs
for i=1:size(PRESSobs,2)
    eta=BAKER_press2eta(squeeze(PRESSobs(tstr*Hz:end,i)),Fs,offset,maxfac);
    %Hsobs(i)=4*std(eta);
    [See,f,Seec] = pwelch(eta,WL*Hz,OL*Hz,[],Hz);
    %[~,ixTp]=max(abs(See));
    %Tpobs_off(i)=1/f(ixTp);
    [~,f2]=min(abs(f-1.2));
    [~,f1]=min(abs(f-0.25));
    Hsobs(i)=4*sqrt(trapz(f(f1:f2),See(f1:f2)));
end

%% MODEL 
for i=1:2
    if i==1
        stname=stname1;
    elseif i==2
        stname=stname2;
    end
    nc=netcdf(stname,'r');
    xpos=nc{'Xgrid'}(:);
    ypos=nc{'Ygrid'}(:);
    sta=find(((xpos==120)));
    % x-positions
    xpos=nc{'Xgrid'}(:);
    xpos=xpos(sta);
    % y-positons
    ypos=nc{'Ygrid'}(:);
    ypos=ypos(sta);
    tend=length(nc{'scrum_time'}(:));
    zeta1=nc{'zeta'}(tstr:tend,:);
    zeta1=zeta1(:,sta);
    time=nc{'scrum_time'}(tstr:tend);
    clear hs
    for k=1:length(sta)
        [S,f] = mycspd(zeta1(:,k),zeta1(:,k),64,10);
        [~,f2]=min(abs(f-1.2));
        [~,f1]=min(abs(f-0.25));
        hs(k)=4*sqrt(trapz(f(f1:f2),S(f1:f2)));
    end
    %hs=4*std(zeta1);
    if i==1
        hs1=hs;
    elseif i==2;
        hs2=hs;
    end
end

%% PLOTTING OFFSHORE ARRAY (X~27m)

figure('Position',[500 500 1400 500]);
plot(Yobs,Hsobs,'o','Color','k','LineWidth',2);
hold on
xlabel('$y (m)$','Interpreter','latex','FontSize',15);
ylabel('$H_s$ $(m)$','Interpreter','latex','FontSize',15);
grid("minor");
ylim([0. 0.3])
plot(ypos/10-15,hs2,'*-','LineWidth',3,'Color',[111,160,135,0.7*255]/255);
plot(ypos/10-15,hs1,'*-','LineWidth',3,'Color',[192,24,135,0.7*255]/255);
legend('Data','Default','Corrected', ...
    'Interpreter','latex','FontSize',15,'Location','northwest');
title("Offshore",'Interpreter','latex','FontSize',18);
print(gcf, '-dpdf','-bestfit', './Figures/Baker-hs_innershelf.pdf');