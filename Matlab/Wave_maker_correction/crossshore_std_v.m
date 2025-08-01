clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%
fname1 = '/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/rip_avg_out.nc';
fname2 = '/Users/simon/Code/CONFIGS/IB09_randomphase/rip_avgtot.nc';   % croco history file name
fname1 = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname2 = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';
nbF = 2 % nb of files
makepdf   = 0;                       % make pdf file
%
%===============================================================
model3D=1;

for i=1:nbF
    if i==1
        fname=fname1;
    elseif i==2
        fname=fname2;
    end
    nc=netcdf(fname);
    tstr=5 ;
    tend=length(nc{'scrum_time'}(:));
    if tend<tstr; tstr=1;  end;
    
    xl=nc{'xl'}(:);
    x=nc{'x_rho'}(1,:)-xl;
   
    ubar=squeeze(nc{'ubar'}(tstr:tend,:,:));
    vbar=squeeze(nc{'vbar'}(tstr:tend,:,:));
    pm=squeeze(nc{'pm'}(:));
    pn=squeeze(nc{'pn'}(:));
    for t=1:(tend-tstr)
        w(t,:,:)=vorticity(squeeze(ubar(t,:,:)),squeeze(vbar(t,:,:)),pm,pn);
    end
    
    sigmaV=squeeze(std(w));
    sigmaVint(i,:)=std(sigmaV);
end

x=(x(2:end)+x(1:end-1))/2;
hold on
for i=1:nbF
    plot(x,squeeze(sigmaVint(i,:)));
end
legend('Without correction','With correction');
