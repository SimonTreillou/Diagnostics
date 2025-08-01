%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First moment plot
%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname1     = '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/rip_histot.nc';   % croco history file name
fname2     = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_histot.nc';   % croco history file name

makepdf   = 0;                % make pdf file
trac      = "tpas01";         % which tracer temperature or passive tracer
%
%======================================================================

nfiles = 2;
g = 9.81;
figure('Position',[100 100 1000 400])
% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

for i=1:2
    if i==1
        fname=fname1;
    elseif i==2
        fname=fname2;
    end
    nc=netcdf(fname);
    time=nc{'scrum_time'}(:);
    y=squeeze(nc{'y_rho'}(:,1));
    x=squeeze(nc{'x_rho'}(1,:));
    [~,ix]=min(abs(y-60));

    t=squeeze(nc{'tpas01'}(50:end,10,1:end,:));
    Dbar=squeeze(mean(t,1));
    tmp1=trapz(x,Dbar.*x,2);
    tmp2=trapz(x,Dbar,2);
    mu=tmp1./tmp2;
    
    plot(y(ix:end-5),mu(ix:end-5)-300,'LineWidth',3);
    hold on
end

ylabel("$\mu$ (m)",'Interpreter','latex','FontSize',14);
xlabel("$y$ (m)",'Interpreter','latex','FontSize',14);
grid();
legend('$\sigma_{\theta}$=30','$\sigma_{\theta}$=10','Interpreter','latex')
