clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%
fname3     = '/Users/simon/Code/IB09/IB09_fini/rip_avg.nc';    % croco file name
fname1='/Users/simon/Code/IB09/IB09PATPSdx1bathynewang7strat/rip_avg.nc';
fname1='/Users/simon/Code/testAVDV/testVADV/rip_avg.nc';
%fname1='/Users/simon/Code/IB09/IB09_PSbienposang17deb3sansspread2sansang/rip_avg.nc';
fname2='/Users/simon/Code/testAVDV/testVADVsans/rip_avg.nc';
%fname2='/Users/simon/Code/IB09/IB09_PSbienposang17deb3sansspread2/rip_avg.nc';
fname3='/Users/simon/Code/testAVDV/testVADVintweno/rip_avg.nc';
%fname     = '../rip_avg_2D_SC.nc';
%fname     = '../rip_avg_2D_LC.nc';

makepdf   = 0;                       % make pdf file
%
%===============================================================
model3D=1;


for i=1:3
    if i==1
        fname=fname1;
        c="red";
    elseif i==2
        fname=fname2;
        c="blue";
    elseif i==3
        fname=fname3;
        c="green";
    end
    nc=netcdf(fname);
    tstr=15;
    tend=length(nc{'scrum_time'}(:));
    if tend<tstr; tstr=1;  end;
    xl=nc{'xl'}(:);
    x=nc{'x_rho'}(1,:)-xl+50;
    xu=squeeze((x(2:end)+x(1:end-1))*0.5);

    
    ubar=squeeze(nc{'ubar'}(:,:,:));
    vbar=squeeze(nc{'vbar'}(:,:,:));
    pm=squeeze(nc{'pm'}(:));
    pn=squeeze(nc{'pn'}(:));
    
    w=zeros(length(tstr:tend),size(vbar,2),size(ubar,3));
    for t=tstr:tend
        ubart=squeeze(ubar(t,:,:));
        vbart=squeeze(vbar(t,:,:));
        w(t,:,:)=vorticity(ubart,vbart,pm,pn);
    end
    
    Jstr=300;
    mw=squeeze(mean(w(:,Jstr:end,:).^2,2));
    tmw=squeeze(mean(mw));
    stdmw=squeeze(std(mw));
    mwh=tmw+stdmw;
    mwb=tmw-stdmw;

    x2 = [xu, fliplr(xu)];
    inBetween = [mwh, fliplr(mwb)];
    m=fill(x2, inBetween,c,'linestyle','none');
    set(m,'facealpha',.1);
    hold on;
    plot(xu,smooth(tmw),'Color',c,'LineWidth',3);
    xlim([-150,0]);
    xlabel("$x (m)$",'Interpreter','latex');
    ylabel("$\epsilon = <\omega^2>$",'Interpreter','latex');

end