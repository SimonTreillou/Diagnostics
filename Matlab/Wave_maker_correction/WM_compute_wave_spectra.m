%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute wave spectra
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
% Wave spectrum input
wd=15;      % Incidence angle
wds=10;     % Directional spread
gamma=20.0; % Peak-enhancement parameter
wa=0.23;    % Wave amplitude
wp=13.0;    % Peak period

% Files
hname1="WM-Corrected-S10";
hname2="WM-Default-S10";

for it=1:2
    if it==1
        hname=hname1;
        S1=0;
    elseif it==2
        hname=hname2;
        S2=0;
    end

    stname  = '/Users/simon/Code/CONFIGS/'+hname+'/stations.nc';
    nc=netcdf(stname,'r');
    Tinit=100;
    
    %
    % --- extracting data ---
    %

    xpos=nc{'Xgrid'}(:);
    ypos=nc{'Ygrid'}(:);
    sta=find(xpos==20);
    tend=length(nc{'scrum_time'});
    %xpos=nc{'Xgrid'}(sta:sta+ext);
    for j=1:1
        xpos=nc{'Xgrid'}(sta(j));
        ypos=nc{'Ygrid'}(sta(j));
        
        %zeta1=nc{'zeta'}(100:tend,sta:sta+ext);
        zeta1=nc{'zeta'}(Tinit:tend,sta(j));
        %u1=squeeze(nc{'u'}(100:tend,sta:sta+ext,10));
        u1=squeeze(nc{'u'}(Tinit:tend,sta(j),10));
        %v1=squeeze(nc{'v'}(100:tend,sta:sta+ext,10));
        v1=squeeze(nc{'v'}(Tinit:tend,sta(j),10));
        time=nc{'scrum_time'}(Tinit:tend);
        depth=6.95;

        
        %
        % --- converting data ---
        %
        % DATA STRUCTURE
        ID.data=[zeta1 u1 v1];
        xtmp=zeros(size(ID.data,2),1);
        ytmp=repmat((ypos-ypos(1)),size(ID.data,2)/size(ypos,1),1);
        ztmp=zeros(size(ID.data,2),1);
        ID.layout=[xtmp'; ytmp'; ztmp'];
        ztmp={'elev'};utmp={'velx'};vtmp={'vely'};
        for i=1:size(zeta1,2)-1
            ztmp{end+1}='elev';utmp{end+1}='velx';vtmp{end+1}='vely';
        end
        ID.datatypes=[ztmp utmp vtmp];
        ID.depth=depth;
        ID.fs=1/(time(2)-time(1));
        
        check_data(ID,1)
        
        % SPECTRAL MATRIX STRUCTURE
        nf=1000;
        nd=1000;
        SM.freqs=linspace(0.01,0.25,nf);
        SM.dirs=linspace(-180,180,nd);
        SM.S=zeros(nf,nd);
        SM.axisdir=0.0;
        SM.funit='Hz';
        SM.dunit='cart';
        
        check_data(SM,2)
        
        % ESTIMATION PARAMETER STRUCTURE
        EP.method='EMEP';
        EP.dres=180;
        EP.iter=100;
        EP.smooth='ON';
        EP.nfft=256;
        
        check_data(EP,3)
        
        % COMPUTING SPECTRUM
        if it==1
            [SMout1_tmp,EPout1_tmp]=dirspec(ID,SM,EP,{'MESSAGE',0,'PLOTTYPE',0});
            freqs1=SMout1_tmp.freqs;
            dirs1=SMout1_tmp.dirs;
            S1=S1+SMout1_tmp.S;
        elseif it==2
            [SMout2_tmp,EPout2_tmp]=dirspec(ID,SM,EP,{'MESSAGE',0,'PLOTTYPE',0});
            freqs2=SMout2_tmp.freqs;
            dirs2=SMout2_tmp.dirs;
            S2=S2+SMout2_tmp.S;
        end
    end
        
end

%% PLOTTING
SMout=SMout1_tmp;
EPout=EPout1_tmp;
SP = true; % choose if one or 4 figures
if SP; figure("Name",hname,'Position', [10 10 900 600]); end;

% 3D frequency-direction spectrum
if SP; subplot(221); else; figure(1); end;
plotspec(SMout, 1);
ylim([-90,90]);
xlim([0.05,0.25]);
%zlim([0.0 0.3*pi/180]);
colorbar();
%caxis([0,0.3*pi/180]);
title('3D wave spectrum at X='+string(xpos)+' and Y='+string(ypos));

% Polar frequency-direction spectrum
if SP; subplot(222); else; figure(2); end;
plotspec(SMout, 2);
title('Polar wave spectrum at X='+string(xpos)+' and Y='+string(ypos));


% Directional spectrum
if SP; subplot(223); else; figure(3); end;
dmin=(wd-80)*(pi/180);
dmax=(wd+80)*(pi/180);
Ndir=400;
wd_bry=linspace(dmin,dmax,Ndir);
wd=15;
wd=wd*pi/180;
wds=10.0;
wds=wds*pi/180;
wa_bry_d=exp(-((wd_bry-wd)/max(1.5*wds,1.e-12)).^2);
%wa_bry_d=wa_bry_d/sum(wa_bry_d*(wd_bry(2)-wd_bry(1)));
wa_bry_d=wa_bry_d/max(wa_bry_d);
DSpecObs=trapz(SMout.freqs,SMout.S);
plot(wd_bry*180/pi,wa_bry_d,'LineWidth',2)
hold on
plot(SMout.dirs,DSpecObs/max(DSpecObs),'LineWidth',2);
title('Directional spectrum at X='+string(xpos)+' and Y='+string(ypos))
legend('Theoretical','Observed');
ylabel('Normalized power spectra','Interpreter','latex');
xlabel('$\theta$ ()','Interpreter','latex');

% Frequency spectrum
if SP; subplot(224); else; figure(4); end;
alpha=0.005;
Nf=50;
[w,S]=jonswap_spectrum(alpha,Nf,wp,gamma,0);
dw=(w(2)-w(1));
S=(wa*sqrt(8))^2/16*S/sum(S*dw);
disp(sum(wa*S)*dw);
FSpecObs=trapz(SMout.dirs,SMout.S');
%FSpecObs=sum(SMout.S'*(SMout.dirs(2)-SMout.dirs(1)));
disp(sum(FSpecObs)*(SMout.freqs(2)-SMout.freqs(1)));
plot(SMout.freqs,FSpecObs,'LineWidth',2);
hold on
plot(w,S,'LineWidth',2); %see Suanda 2016;
title('Frequency spectrum at X='+string(xpos)+' and Y='+string(ypos))
legend('Observed','Theoretical');
ylabel('Normalized power spectra','Interpreter','latex');
xlabel('$f$ ($Hz$)','Interpreter','latex');

% 
% fname  = '/Users/simon/Code/CONFIGS/'+hname+'/rip_his.nc';
% nc=netcdf(fname,'r');
% v=nc{'v'}(end,10,:,:);
% x=squeeze(nc{'x_rho'}(1,:));
% y=squeeze(nc{'y_rho'}(:,1));
% y=(y(2:end)+y(1:end-1))/2;
% nc.close();
% nc=netcdf(stname,'r');
% 
% figure('Position',[10 10 900 1200]);
% contourf(x-300,y,squeeze(v));
% hold on
% ylabel('Y-grid points');
% xlabel('X-grid points');
% colorbar();
% plot(nc{'Xgrid'}(:),nc{'Ygrid'}(:),'o','MarkerSize',8,'MarkerFaceColor','red','MarkerEdgeColor','red');
% text(nc{'Xgrid'}(:)-4,nc{'Ygrid'}(:)+40,string(1:length(nc{'Xgrid'}(:))),'BackgroundColor','white');
%     plot(nc{'Xgrid'}(1:5),nc{'Ygrid'}(1:5),'o','MarkerSize',8,'MarkerFaceColor','red','MarkerEdgeColor','red');
%     text(nc{'Xgrid'}(1)+7,nc{'Ygrid'}(1),"Station (150,50)",'BackgroundColor','white');
%     plot(nc{'Xgrid'}(6:10),nc{'Ygrid'}(6:10),'o','MarkerSize',8,'MarkerFaceColor','red','MarkerEdgeColor','red');
%     text(nc{'Xgrid'}(6)+7,nc{'Ygrid'}(6),"Station (150,200)",'BackgroundColor','white');
%     plot(nc{'Xgrid'}(11:15),nc{'Ygrid'}(11:15),'o','MarkerSize',8,'MarkerFaceColor','red','MarkerEdgeColor','red');
%     text(nc{'Xgrid'}(11)+7,nc{'Ygrid'}(11),"Station (150,350)",'BackgroundColor','white');
% title('Surface $v$','Interpreter','latex');

%% Directional spread comparison
figure('PaperPositionMode', 'auto');clear alpha
dmin=(wd-80)*(pi/180);
dmax=(wd+80)*(pi/180);
Ndir=400;
wd_bry=linspace(dmin,dmax,Ndir);
wd=15;
wd=wd*pi/180;
wds=10.0;
wds=wds*pi/180;
wa_bry_d=exp(-((wd_bry-wd)/max(1.5*wds,1.e-12)).^2);
%wa_bry_d=wa_bry_d/sum(wa_bry_d*(wd_bry(2)-wd_bry(1)));
wa_bry_d=wa_bry_d/max(wa_bry_d);
DSpecObs1=trapz(SMout1.freqs,SMout1.S);
DSpecObs2=trapz(SMout2.freqs,SMout2.S);
scatter(wd_bry*180/pi,wa_bry_d,'o','k')
hold on
plot(SMout2.dirs,DSpecObs2/max(DSpecObs2),'LineWidth',3, ...
    'Color',[111,160,135,0.7*255]/255);
plot(SMout1.dirs,DSpecObs1/max(DSpecObs1),'LineWidth',3, ...
    'Color',[192,24,135,0.7*255]/255);
grid();
%title('Directional spectrum at X='+string(xpos)+' and Y='+string(ypos))
legend('Input','Default','Corrected','Interpreter','latex','FontSize',15);
ylabel('Directional spread $D(\theta)$','Interpreter','latex','FontSize',15);
xlabel('$\theta$ (deg)','Interpreter','latex','FontSize',15);
xlim([-40+wd*180/pi,40+wd*180/pi]);
ylim([0 1.1]);
print(gcf, '-dpdf', './Figures/dspread-comparison.pdf');


%% Frequency spectrum
figure('PaperPositionMode', 'auto');
alpha=0.005;
Nf=500;
[w,S]=jonswap_spectrum(alpha,Nf,wp,gamma,0);
dw=(w(2)-w(1));
S=(wa*sqrt(8))^2/16*S/sum(S*dw);
FSpecObs1=trapz(SMout1.dirs,SMout1.S');
FSpecObs2=trapz(SMout2.dirs,SMout2.S');
scatter(w,S,'.','k'); %see Suanda 2016;
hold on
plot(SMout2.freqs,FSpecObs2,'LineWidth',3, ...
    'Color',[111,160,135,0.7*255]/255);
plot(SMout1.freqs,FSpecObs1,'LineWidth',3, ...
    'Color',[192,24,135,0.7*255]/255);
grid("minor");
%title('Frequency spectrum at X='+string(xpos)+' and Y='+string(ypos))
legend('Input','Default','Corrected','Interpreter','latex','FontSize',15);
ylabel('Frequency spectrum $D(f)$','Interpreter','latex','FontSize',15);
xlabel('$f$ (Hz)','Interpreter','latex','FontSize',15);
xlim([0.05 0.20]);
print(gcf, '-dpdf', './Figures/freq-comparison-nfft256.pdf');


%%

alpha=0.005;
Nf=500;
[w,S]=jonswap_spectrum(alpha,Nf,wp,gamma,0);
dw=(w(2)-w(1));
S=(wa*sqrt(8))^2/16*S/sum(S*dw);

scatter(w,S,70,'.','k'); %see Suanda 2016;
hold on
plot(freqs1,trapz(dirs2,S2'/19),'LineWidth',3, ...
    'Color',[111,160,135,0.7*255]/255);
plot(freqs2,trapz(dirs1,S1'/19),'LineWidth',3, ...
    'Color',[192,24,135,0.7*255]/255);
grid("minor");
%title('Frequency spectrum at X='+string(xpos)+' and Y='+string(ypos))
legend('Input','Default','Corrected','Interpreter','latex','FontSize',15);
ylabel('Frequency spectrum $D(f)$','Interpreter','latex','FontSize',15);
xlabel('$f$ (Hz)','Interpreter','latex','FontSize',15);
xlim([0.05 0.20]);
