%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname = '/Users/simon/Code/IB09/testVADVC2periodok/rip_his.nc';

makepdf   = 0;                % make pdf file
meanlongs = 1;                % allow to reduce incertainty
%
%======================================================================
yindex = 350;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tstr=1;
tend=length(nc{'scrum_time'}(:)); % reads last record
y=squeeze(nc{'y_rho'}(:,1));
if meanlongs
    ystr=1;
    yend=400; %length(y);
    ystep=30; 
else
    ystr=yindex;
    yend=yindex;
    ystep=1;
end

%% COMPUTING ETA AND ITS DERIVATIVE
n=1;
for yindex=ystr:ystep:yend
    disp("Y="+string(y(yindex)));
    for tindex=tstr:tend % ---------------------------------------------
     hr=squeeze(nc{'h'}(yindex,:));
     xr=squeeze(nc{'x_rho'}(yindex,:));
     dx=xr(2)-xr(1);
     Dcrit=nc{'Dcrit'}(:);
    
    % vertical grid
     N=length(nc('s_rho'));
     theta_s=nc.theta_s(:); 
     theta_b=nc.theta_b(:); 
     hc=nc.hc(:); 
     zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
     zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,N,'w',2));
    
     D=hr+zeta;
     ztop=squeeze(zw(end,:,:));
     ztop(D<Dcrit+0.01)=NaN;
     Z(tindex,:)=ztop'; 
     dZdx(tindex,:)=(ztop(2:end)'-ztop(1:end-1)')/dx;
    end
    
    % COMPUTING BREAKER POSITION
    x=squeeze(nc{'x_rho'}(1,:))-300;
    time=squeeze(nc{'scrum_time'}(:));
    [~,x50]=min(abs(x+50));
    for i=1:size(Z,1)
        [~,index]=min(abs(dZdx(i,1:x50)-5000));
        indexes(i)=index;
    end
    if meanlongs
        x_breaker(n,:) = x(indexes);
        n=n+1;
    else
        x_breaker = x(indexes);
    end
end

%% BREAKING POSITION AGAINST TIME
if meanlongs
    x_breakerM=mean(x_breaker,1);
else
    x_breakerM=x_breaker;
end
figure();
plot(time,x_breakerM,'color','b','LineWidth',1);
hold on;
xb_smooth=smooth(x_breakerM(1:end,:),1);
plot(time,xb_smooth,'color','r','LineWidth',2);
grid();
xlabel("$t$ ($s$)",'Interpreter','latex');
ylabel("$x$ ($m$)",'Interpreter','latex');
legend('Raw','Smoothed','Location','best');
ylim([-300,0]);
title('Evolution of the breaker position (maximum wave steepness) for y='+string(y(yindex)))
subtitle('Mean breaker position: '+string(mean(xb_smooth(10:end-10)))+"$\pm$"+string(std(xb_smooth(10:end-10))),'Interpreter','latex');

%% VISUAL VERIFICATION

% Chosen timestep
tex=19;

figure();
title('Timestep='+string(time(tex)));

subplot(2,1,1);
plot((xr(2:end)+xr(1:end-1))*0.5-xr(1,end)+50,dZdx(tex,:),'color','r','LineWidth',2);
grid on
xlim([-150,0]);
ylabel('$\frac{\partial \eta}{\partial x}$','Interpreter','latex');
xlabel('$x$ ($m$)','Interpreter','latex');
title('$x$-gradient of the free surface','Interpreter','latex');

subplot(2,1,2);
plot(xr-xr(1,end)+50,-hr,'color','k','LineWidth',3);
hold on
hn=plot(xr-xr(1,end)+50,Z(tex,:),'color','r','LineWidth',2);
hn=plot(x_breaker(tex)*(xr./xr),linspace(-10,2,size(xr,2)),'--b','LineWidth',2);
grid on
xlim([-150,0]);
title('Alongshore free surface','Interpreter','latex');
ylabel('$\eta$ ($m$)','Interpreter','latex');
xlabel('$x$ ($m$)','Interpreter','latex');
legend('Bathymetry','Free surface','Detected breaker limit','Location','best');

