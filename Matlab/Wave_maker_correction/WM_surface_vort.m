%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
%fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_his.nc';
makemovie = 1;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
%
%======================================================================

yindex = 101;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record
time=nc{'scrum_time'}(:);


if makemovie
    tstr=1;
    tend=tindex;
    v = VideoWriter("./Figures/Vorticity_Default_S10", 'MPEG-4');
    v.FrameRate=5;  % nombre d'image par seconde
    open(v);
else
    tstr=tindex;
    tend=tstr;
end

hf = figure('position',[1000 500 400 800]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');


rgb = [ ...
    94    79   162
    50   136   189
   102   194   165
   171   221   164
   230   245   152
   255   255   255
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;


for tindex=tstr:tend % ---------------------------------------------
    ubar=squeeze(nc{'ubar'}(tindex,:,:,10));
    vbar=squeeze(nc{'vbar'}(tindex,:,:,10));
    pm=squeeze(nc{'pm'}(:));
    pn=squeeze(nc{'pn'}(:));

    w=vorticity(ubar,vbar,pm,pn);
    xr=squeeze(nc{'x_rho'}(:))-300;
    yr=squeeze(nc{'y_rho'}(:));
    [P,L]=size(yr);
    w(:,2:L-1)=0.5*(w(:,1:L-2)+w(:,2:L-1));
    w(:,1)=w(:,2);w(:,L)=w(:,L-1);

    clf
    hold on

    w(2:P-1,:)=0.5*(w(1:P-2,:)+w(2:P-1,:));
    w(1,:)=w(2,:);w(P,:)=w(P-1,:);
    cmin=-0.1;
    cmax=-cmin;
    cint=(cmax-cmin)/20;
    h=pcolor(xr(:,:),yr(:,:),w(:,:));
    set(h,'EdgeColor','none');
    caxis([cmin cmax]);
    cmap=redwhiteblue(cmin,cmax);
    %colormap(rgb);
    colormap(cmap);
    colorbar();
    ylim([min(yr(:,1)),max(yr(:,1))]);
    xlim([min(xr(1,:)),max(xr(1,:))]);
    set(gcf, 'Color','w');
    hh=floor(time(tindex)/3600);
    mm=round((time(tindex)-hh*3600)/60);
    title("Vertical vorticity at t="+string(hh)+'h'+string(mm)+'min', ...
        'Interpreter','latex','FontSize',20);
    ylabel('$y$ ($m$)','Interpreter','latex','FontSize',16);
    xlabel('$x$ ($m$)','Interpreter','latex','FontSize',16);

    if makemovie
        frame=getframe(gcf);
        %frame.cdata=frame.cdata(1:xframe,1:yframe,:);
        writeVideo(v,frame);
    else
        drawnow
        pause(.2)
    end
end

if makemovie
 close(v)
end
