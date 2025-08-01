%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code highly inspired by Rene Schubert Python's code
% (https://github.com/reneschubert/keflux). This code computes the kinetic
% energy flux with the coarse-graining approach.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ====================
%
% --- model params ---
%
fname='/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/rip_avg_out.nc'; 
fname='/Users/simon/Code/CONFIGS/IB09_randomphase_S30/rip_avgtot.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
%fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';
%fname = '/Users/simon/Code/CONFIGS/IB09_S10/rip_avg.nc';
%fname = '/Users/simon/Code/CONFIGS/IB09_DB/rip_avg.nc';


%
%===============================================================

tstr=5;  % Starting time-step
Jstr=1; % Alongshore starting point

nc=netcdf(fname);
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1;  end;

xl=nc{'xl'}(:);
x=squeeze(nc{'x_rho'}(1,:))-xl;
y=squeeze(nc{'y_rho'}(Jstr:end,1));
xplt=(x(2:end)+x(1:end-1))/2;xplt=(xplt(2:end)+xplt(1:end-1))/2;
xplt=(xplt(2:end)+xplt(1:end-1))/2;
yplt=(y(2:end)+y(1:end-1))/2;yplt=(yplt(2:end)+yplt(1:end-1))/2;
yplt=(yplt(2:end)+yplt(1:end-1))/2;
rho=squeeze(nc{'rho'}(:));

%% COMPUTING ENERGY FLUX SCALE EF(L,t,y,x)
gs=1;
scales=[2,3,4,5,10,20,30,40,50,75,100];
%scales=[1,2,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,80,100];
EF=zeros(length(scales),length(tstr:tend),length(y)-3,length(x)-3);
n=1;
for t=tstr:tend
    v=squeeze(nc{'vbar'}(t,Jstr:end,:));
    v=(v(:,2:end)+v(:,1:end-1))/2;
    u=squeeze(nc{'ubar'}(t,Jstr:end,:));
    u=(u(2:end,:)+u(1:end-1,:))/2;
    EF(:,n,:,:)=scalekineticflux_CG(x,y,u,v,scales,gs);
    n=n+1;
end

%% PLOTTING EF(X,Y)
nS=7;
var=1e3*squeeze(sum(EF(nS,:,:,:))/(tend-tstr)); % AVERAGING ON TIME AT SCALE nS
%var=1e3*squeeze(EF(nS,end,:,:));

figure('Units','pixels','Position',[500 500 500 800]);
set(gca,'FontSize',15)
h=pcolor(xplt,yplt,var);  
set(h,'EdgeColor','none');
l=max(max(abs(var)));
cmap=redwhiteblue(-l,l);
colormap(cmap);
caxis([-l l])
colorbar();
xlim([min(x)+scales(nS)/2+1,max(x)-scales(nS)/2-1])
ylim([min(y)+scales(nS)/2+1,max(y)-scales(nS)/2-1])
xlabel('$ x (m) $','Interpreter','latex','FontSize',15);
ylabel('$ y (m) $','Interpreter','latex','FontSize',15);
title("Kinetic energy flux $\Pi$ at "+string(scales(nS))+"m ($mW/m^3$)",'Interpreter','latex','FontSize',15);

%% PLOTTING <EF>(L)
spec=1e3*squeeze(sum(EF(:,:,:,:),2)/(tend-tstr));
casc=mean(mean(spec,3),2);
%cascstd=std(std(spec,3),2);
semilogx(1./scales,smooth(casc,1));
hold on
line([2.e-3 1],[0 0],'color','k')
grid();
ylabel('Spectral Energy Flux m^{2}.s^{-3}');
%axis([6 800 -1.e-4 0.5e-4])

%%
T=(casc(2:end)-casc(1:end-1))./(scales(2:end)-scales(1:end-1))';
semilogx(scales(2:end),T);