%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure 3 from Hally-Rosendahl & Feddersen, 2016
% Significant wave height Hs and alongshore average velocity V
% Simon Treillou, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
stname1 = '/Users/simon/Code/CONFIGS/BAKER_G2B/stations.nc';
stname1 = '/Users/simon/Code/CONFIGS/BAKER_G1D/stations.nc';
stname2 = '/Users/simon/Code/CONFIGS/BAKER_G2B_default/stations.nc';

for i=1:1
    if i==1
        stname=stname1;
    elseif i==2
        stname=stname2;
    end
    nc=netcdf(stname,'r');
    xpos=nc{'Xgrid'}(:);
    ypos=nc{'Ygrid'}(:);
    sta_crossshore=find(((ypos==150)));
    %sta=find(((xpos==120)));
    time=nc{'scrum_time'}(:);
    dt=time(2)-time(1);
    [~,tstr]=min(abs(time<15*60));
    [~,tend]=min(abs(time<25*60));

    for j=1:length(sta_crossshore)
        xpos=nc{'Xgrid'}(:);
        sta=find((xpos==xpos(sta_crossshore(j))));

        
        % x-positions
        xpos=nc{'Xgrid'}(:);
        xpos=xpos(sta);
        % y-positons
        ypos=nc{'Ygrid'}(:);
        ypos=ypos(sta);

        tend=length(nc{'scrum_time'}(:));
        zeta1=nc{'zeta'}(tstr:tend,:);
        zeta1=zeta1(:,sta);
        clear hs
        for ii=1:length(sta)
            zz=squeeze(zeta1(:,ii));
            [S,f]=mycspd(zz,zz,2^6,1/dt);
            [~,f1]=min(abs(f-0.3));
            [~,f2]=min(abs(f-2));
            hs(ii)=4*sqrt(trapz(f(f1:f2),S(f1:f2)));
        end
        mhs=mean(hs);
        shs=std(hs);
    
        if i==1
            hs1(j)=mhs;
            hs1_std(j)=shs;
        elseif i==2;
            hs2(j)=mhs;
            hs2_std(j)=shs;
        end
    end
    xpos=nc{'Xgrid'}(:);
    ypos=nc{'Ygrid'}(:);
    sta=find(ypos==150);
    xpos=xpos(sta);

end


fname ='/Users/simon/Code/CONFIGS/BAKER_G1D/rip_his.nc';
nc=netcdf(fname);
h=nc{'h'}(1,:);
x=nc{'x_rho'}(1,:);
dx=1./nc{'pm'}(1,1);
dy=1./nc{'pn'}(1,1);
% zeta=nc{'zeta'}(end-4,200,:);
% load /Users/simon/Code/Matlab/Wave_maker_correction/BAKER_bathy
% plot(x+15,-h+1.07,'LineWidth',2,'Color','k');
% hold on
% plot(x+15,zeta+1.07,'LineWidth',1.5,'Color','red','LineStyle','--');
% plot(x+15,ones(size(x))*1.07,'LineStyle','--','LineWidth',1.5,'Color','k')
% ylim([-0.05 1.5]);
% xlim([19 35]);
% legend('CROCO bathy','Real bathy','Interpreter','latex','FontSize',15)
% xlabel('$x (m)$','Interpreter','latex','FontSize',15);
% ylabel('$z$ $(m)$','Interpreter','latex','FontSize',15);


%% DATA
% TR="G1D";
% timeobs=0:1:(45*60)*100;
% [~,tinit]=min(abs(timeobs-15*60*100));
% [~,tend]=min(abs(timeobs-25*60*100));
% 
% 
% trial=TR+"-IS";
% vr='wg';
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'1.txt')
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'2.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'3.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'4.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'5.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'6.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'7.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'8.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'9.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'10.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'11.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'12.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'13.txt')
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'14.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/'+vr+'15.txt');
% 
% WGobs1=[wg1 wg2 wg3 wg4 wg5 wg6 wg7];
% WGobs2=[wg9 wg10 wg11 wg12 wg13 wg14 wg15];
% Hz=100;
% WL=2^5;
% OL=2^4;
% 
% for i=1:size(WGobs1,2)
%     [See,f,Seec] = pwelch(WGobs1(tinit:tend,i),WL*Hz,OL*Hz,[],Hz);
%     [~,ixTp]=max(abs(See));
%     %Tpobs_off(i)=1/f(ixTp);
%     [~,f1]=min(abs(f-0.3));
%     [~,f2]=min(abs(f-2));
%     Hsobs_off1(i)=4*sqrt(trapz(f(f1:f2),See(f1:f2)));
% end
% Hsobs_off1_std=std(Hsobs_off1);
% Hsobs_off1=sum(Hsobs_off1)/length(Hsobs_off1);
% 
% for i=1:size(WGobs2,2)
%     [See,f,Seec] = pwelch(WGobs2(tinit:tend,i),WL*Hz,OL*Hz,[],Hz);
%     [~,ixTp]=max(abs(See));
%     %Tpobs_off(i)=1/f(ixTp);
%     [~,f1]=min(abs(f-0.3));
%     [~,f2]=min(abs(f-2));
%     Hsobs_off2(i)=4*sqrt(trapz(f(f1:f2),See(f1:f2)));
% end
% Hsobs_off2_std=std(Hsobs_off2);
% Hsobs_off2=sum(Hsobs_off2)/length(Hsobs_off2);
% 
% 
% trial=TR+"-IS";
% load('/Users/simon/Code/BAKER/'+trial+'/press1.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press2.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press3.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press4.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press5.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press6.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press7.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press8.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press9.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press10.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press11.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press12.txt');
% 
% PobsIS1=[press3];
% PobsIS2=[press4 press5 press6 press7 press8 press9 press10 press11 press12];
% 
% Hz = 100;
% maxfac = 1.2;
% offset = 0.05;
% Fs=100; 
% 
% for i=1:size(PobsIS1,2)
%     eta=BAKER_press2eta(squeeze(PobsIS1(tinit:tend,i)),Fs,offset,maxfac);
%     [S,f]=mycspd(eta,eta,2^8,Hz);
%     [~,f1]=min(abs(f-0.3));
%     [~,f2]=min(abs(f-2));
%     Hsobs_IS1(i)=4*sqrt(trapz(f(f1:f2),S(f1:f2)));
%    % Hsobs_IS1(i)=4*std(eta);
% end
% Hsobs_IS1_std=std(Hsobs_IS1);
% Hsobs_IS1=sum(Hsobs_IS1)/length(Hsobs_IS1);
% 
% for i=1:size(PobsIS2,2)
%     eta=BAKER_press2eta(squeeze(PobsIS2(tinit:tend,i)),Fs,offset,maxfac);
%     [S,f]=mycspd(eta,eta,2^8,Hz);
%     [~,f1]=min(abs(f-0.3));
%     [~,f2]=min(abs(f-2));
%     Hsobs_IS2(i)=4*sqrt(trapz(f(f1:f2),S(f1:f2)));
%     %Hsobs_IS2(i)=4*std(eta);
% end
% Hsobs_IS2_std=std(Hsobs_IS2);
% Hsobs_IS2=sum(Hsobs_IS2)/length(Hsobs_IS2);
% 
% trial=TR+"-SZ";
% load('/Users/simon/Code/BAKER/'+trial+'/press1.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press2.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press3.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press4.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press5.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press6.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press7.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press8.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press9.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press10.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press11.txt');
% load('/Users/simon/Code/BAKER/'+trial+'/press12.txt');
% 
% PobsSZ1=[press1 press2 press3 press4 press5 press6 press7 press8];
% PobsSZ2=[press9 press10 press11 press12];
% 
% for i=1:size(PobsSZ1,2)
%     eta=BAKER_press2eta(squeeze(PobsSZ1(tinit:tend,i)),Fs,offset,maxfac);
%     [S,f]=mycspd(eta,eta,2^8,Hz);
%     [~,f1]=min(abs(f-0.3));
%     [~,f2]=min(abs(f-2));
%     Hsobs_SZ1(i)=4*sqrt(trapz(f(f1:f2),S(f1:f2)));
%     %Hsobs_SZ1(i)=4*std(eta);
% end
% Hsobs_SZ1_std=std(Hsobs_SZ1);
% Hsobs_SZ1=sum(Hsobs_SZ1)/length(Hsobs_SZ1);
% 
% for i=1:size(PobsSZ2,2)
%     eta=BAKER_press2eta(squeeze(PobsSZ2(tinit:tend,i)),Fs,offset,maxfac);
%     [S,f]=mycspd(eta,eta,2^8,Hz);
%     [~,f1]=min(abs(f-0.3));
%     [~,f2]=min(abs(f-2));
%     Hsobs_SZ2(i)=4*sqrt(trapz(f(f1:f2),S(f1:f2)));
%     %Hsobs_SZ2(i)=4*std(eta);
% end
% Hsobs_SZ2_std=std(Hsobs_SZ2);
% Hsobs_SZ2=sum(Hsobs_SZ2)/length(Hsobs_SZ2);
% 
% Hsobs=[Hsobs_off1 Hsobs_off2 Hsobs_IS1 Hsobs_IS2 Hsobs_SZ1 Hsobs_SZ2];
% Hsobs_std=[Hsobs_off1_std Hsobs_off2_std Hsobs_IS1_std Hsobs_IS2_std Hsobs_SZ1_std Hsobs_SZ2_std];
%
%
%% PLOTTING

load ./BAKER_G1D_Hs_Lidar
load ./BAKER_G1D_Hs_InSitu

figure('Position',[500 500 1400 500]);
subplot(2,1,1);
%errorbar(xpos/10+15,hs2,hs2_std,'LineWidth',1.5,'Color',[111,160,135,0.7*255]/255);
hold on
%plot(xpos/10+15,hs2,'-*','LineWidth',3,'Color',[111,160,135,0.7*255]/255);
plot(xpos/10+15,hs1,'-*','LineWidth',3,'Color',[192,24,135,0.7*255]/255);
plot(BAKER_G1D_Hs_Lidar(:,1),BAKER_G1D_Hs_Lidar(:,2),'Color','r', ...
    'LineStyle','-.','LineWidth',2);
scatter(BAKER_G1D_Hs_InSitu(:,1),BAKER_G1D_Hs_InSitu(:,2),400,'.','Color','b');
%errorbar(xpos/10+15,hs1,hs1_std,'LineWidth',1.5,'Color',[192,24,135,0.7*255]/255);
%errorbar(xpos/10+15,Hsobs,Hsobs_std,'o','Color','k','MarkerFaceColor','k','LineWidth',1.5)
legend('Default','Corrected','Lidar','In situ','Interpreter','latex','FontSize',15);
%legend('Default','Corrected','In situ','Interpreter','latex','FontSize',15);
grid("minor");
ylim([0 0.35]);
xlim([19 35]);
%xlabel('$x (m)$','Interpreter','latex','FontSize',15);
ylabel('$H_s$ $(m)$','Interpreter','latex','FontSize',15);

subplot(2,1,2);
load /Users/simon/Code/Matlab/Wave_maker_correction/BAKER_bathy

plot(x+15,-h+1.07,'LineWidth',2,'Color','k');
hold on
plot(BAKER_bathy(:,1),BAKER_bathy(:,2),'LineWidth',1.5,'Color','red','LineStyle','--');
plot(x+15,ones(size(x))*1.07,'LineStyle','--','LineWidth',1.5,'Color','k')
ylim([-0.05 1.5]);
xlim([19 35]);
legend('CROCO bathy','Real bathy','Interpreter','latex','FontSize',15,'Location','best')
xlabel('$x (m)$','Interpreter','latex','FontSize',15);
ylabel('$z$ $(m)$','Interpreter','latex','FontSize',15);

print(gcf, '-dpdf','-bestfit', './Figures/Baker-hs_crossshore.pdf');





%%
fname = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_avgtot.nc';
vlevel = 1;
imin=100;
imax=150;
npts=1024;
[k,E]=get_e_spectra(fname,1,vlevel,imin,npts);

[k2,E2]=get_e_spectra2(fname,1,vlevel,imin,imax,npts);


loglog(k,smooth(E,10));
hold on
loglog(k,k.^(-5/3));
loglog(k2,E2);


