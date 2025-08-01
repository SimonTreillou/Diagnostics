clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%
dirpath   = '/Users/simon/Code/IB09/IB09PATPSdx1bathynewang11stratradiatif/';    % croco file name
%dirpath   = '/Users/simon/Code/IB09/passolo/';    % croco file name
fname     = [dirpath,'rip_avg.nc'];
dname     = [dirpath,'rip_diags_eddy_avg.nc'];
fname1='/Users/simon/Code/testAVDV/testVADV/rip_avg.nc';
fname2='/Users/simon/Code/testAVDV/testVADVsans/rip_avg.nc';
fname3='/Users/simon/Code/testAVDV/testVADVintweno/rip_avg.nc';
dname1='/Users/simon/Code/testAVDV/testVADV/rip_diags_eddy_avg.nc';
dname2='/Users/simon/Code/testAVDV/testVADVsans/rip_diags_eddy_avg.nc';
dname3='/Users/simon/Code/testAVDV/testVADVintweno/rip_diags_eddy_avg.nc';

fname1= '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/rip_histot.nc';   % croco history file name
fname2= '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_histot.nc';   % croco history file name
dname1= '/Users/simon/Code/CONFIGS/IB09_randomphase_S30/rip_diags_eddy_avgtot.nc';   % croco history file name
dname2= '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_diags_eddy_avgtot.nc';   % croco history file name

%fname     = '../rip_avg_2D_SC.nc';
%fname     = '../rip_avg_2D_LC.nc';

makepdf   = 0;                       % make pdf file
%
%===============================================================

for i=1:2
    if i==1
        fname=fname1;
        dname=dname1;
        c="red";
    elseif i==2
        fname=fname2;
        dname=dname2;
        c="blue";
    elseif i==3
        fname=fname3;
        dname=dname3;
        c="green";
    end      
    nc=netcdf(fname);
    ncd=netcdf(dname);
    tstr=5;
    tend=length(nc{'scrum_time'}(:));
    if tend<tstr; tstr=1;  end;
    %tend=88;
    
    h=nc{'h'}(:);
    xl=nc{'xl'}(:);
    x=nc{'x_rho'}(:)-xl+50;
    y=nc{'y_rho'}(:);
    N=length(nc('s_rho'));
    
    u=squeeze(nc{'u'}(tstr:tend,N,:,:));
    v=squeeze(nc{'v'}(tstr:tend,N,:,:));
    


%%
    u=squeeze(mean(nc{'u'}(tstr:tend,:,:,:),2));
    v=squeeze(mean(nc{'v'}(tstr:tend,:,:,:),2));
    for it=1:tend-tstr+1
      ur(it,:,:)=u2rho_2d(squeeze(u(it,:,:)));
      vr(it,:,:)=v2rho_2d(squeeze(v(it,:,:)));
    end
    mu=mean(ur);
    mv=mean(vr);
    
    for it=1:tend-tstr+1
     eke(it,:,:)=0.5*((ur(it,:,:)-mu).^2+(vr(it,:,:)-mv).^2);
    end
    meke=squeeze(mean(eke));

%%
    V=squeeze(mean(mean(nc{'vbar'}(5:end,:,:),1),2));
    stdV=squeeze(std(nc{'vbar'}(5:end,:,:),1));
    stdV=squeeze(std(stdV));

    for it=1:tend-tstr+1
        upuv(it,:,:)=(ur(it,:,:)-mu).*(vr(it,:,:)-mv);
    end

    tt=squeeze(mean(upuv));
    pm=nc{'pm'}(:);
    tx=squeeze((tt(:,2:end)-tt(:,1:end-1))*pm(1,1));
    txr=u2rho_2d(tx);

%%
    uv = ncd{'uv'}(:,10,:,:);
    uvt = squeeze(mean(uv));
    pm = nc{'pm'}(:);
    uvtx = squeeze((uvt(:,2:end)-uvt(:,1:end-1))*pm(1,1));
    uvtxr=u2rho_2d(uvtx);
%plot(x(1,:),smooth(squeeze(mean(uvtxr))));
%%
    subplot(2,1,1)
    Vh = V + stdV';
    Vb = V - stdV';
    x2 = [x(1,:), fliplr(x(1,:))];
    inBetween = [Vh', fliplr(Vb')];
    m=fill(x2, inBetween,c,'linestyle','none');
    set(m,'facealpha',.1)
    hold on
    plot(x(1,:),V,'LineWidth',2,'Color',c);
    ylabel("$V$ ($m.s^{-1}$)",'Interpreter','latex','FontSize',14);
    title("Alongshore velocity",'Interpreter','latex','FontSize',14);
    axis([-300 0 -0.05 0.4])
    grid()
    legend('','INTC2','','Sans modif VADV','','INTWENO');

    subplot(2,1,2)
    txrt=squeeze(-smooth(mean(uvtxr)));
    txrstd=squeeze(std(-uvtxr));
    txrh=txrt'+txrstd;
    txrb=txrt'-txrstd;
    inBetween = [txrh, fliplr(txrb)];
    m=fill(x2, inBetween,c,'linestyle','none');
    %m=fill([x(1,:)';flipud(x(1,:)')], ...
    %    [txrt-txrstd;flipud(txrt+txrstd)],c,'linestyle','none');
    set(m,'facealpha',.1)
    hold on
    plot(x(1,:),txrt,'LineWidth',2,'Color',c);
    axis([-300 0 -0.005 0.005])
    ylabel("$-\frac{\delta}{\delta x}  <u' v'>$ ($m.s^{-2}$)",'Interpreter', ...
        'latex','FontSize',14);
    xlabel("$x$ ($m$)",'Interpreter','latex');
    grid()
    title("Cross-shore eddy advection term ($V$ equation)",'Interpreter','latex', ...
        'FontSize',14);
    legend('','INTC2','','Sans modif VADV','','INTWENO');
    
end

    
