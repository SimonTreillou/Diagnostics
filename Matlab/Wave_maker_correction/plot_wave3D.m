clear all
close all

%====================================================================
make_movie = 0;   % = 1 to make movie
makepdf    = 0;   % = 1 to save figure 

nmfile = 'Wave3D';    % movie name

DIR    = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_histot.nc';
fname  = [DIR,'rip_his.nc'];  % input file name
fname = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_histot.nc'; 
fname0 = fname; %'../zeta_1sec.nc';

Jstr=500;
Jend=1000;

if make_movie,
 t1=100; t2=140;
else
 t1=input('Time index: '); 
 t2=t1; 
end

%====================================================================

disp('Grid ...')

nc=netcdf(fname);
X=squeeze(nc{'x_rho'}(1,:));
Y=squeeze(nc{'y_rho'}(Jstr:Jend,1));
h=squeeze(nc{'h'}(Jstr:Jend,:));
zeta=squeeze(mean(nc{'zeta'}(:,Jstr:Jend,:)));
Dcrit=nc{'Dcrit'}(:);
theta_s=squeeze(nc{'theta_s'}(:));
theta_b=squeeze(nc{'theta_b'}(:));
hc=squeeze(nc{'hc'}(:));
N=length(nc{'sc_r'}(:));
close(nc)
dx=X(2)-X(1);
dy=Y(2)-Y(1);
[M L]=size(h);
%h=repmat(h(1,:),[M 1]);
X3=repmat(X,[M 1]);X3=repmat(X3,[1 1 N]);X3=permute(X3,[3 1 2]);
Y3=repmat(Y,[1 L]);Y3=repmat(Y3,[1 1 N]);Y3=permute(Y3,[3 1 2]);
X2=squeeze(X3(1,:,:));
Y2=squeeze(Y3(1,:,:));
%h3=permute(repmat(h,[1 1 t0]),[3 1 2]);

x=squeeze(X3(1,:,:));
y=squeeze(Y3(1,:,:));
depth=min(min(-h));
X3min=min(min(min(X3)));
X3max=max(max(max(X3)));
Ylimits=sort([Y(1) Y(M)]);

zeta(h<Dcrit)=zeta(h<Dcrit)-h(h<Dcrit); 
%Z=zeros(N,M,L);
zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2);
Z=squeeze(zr);
dz3=squeeze(Z(2,:,:)-Z(1,:,:));
dz3=repmat(dz3,[1 1 N]);
dz3=permute(dz3,[3 1 2]);  % equispacing z

if make_movie,
 nc=netcdf(fname0);
 t=nc{'scrum_time'}(t1:t2);
 close(nc)
end

disp('PLOT')

% Implementing map
numColors = 256;
oceanR = 15;
oceanG = 127;
oceanB = 230;
map=[linspace(oceanR, 255, numColors)'/255,...
  linspace(oceanG, 255, numColors)'/255,...
  linspace(oceanB, 255, numColors)'/255];

%===============================================================
% Make movie

figure;
%xframe=1600; yframe=1600;
%set(gcf,'position',[100 100 xframe yframe])
%set(gcf,'color','white')

if make_movie
    v = VideoWriter(nmfile,'MPEG-4');
    v.FrameRate=3;  % nombre d'image par seconde
    open(v);
end

for it=t1:t2  % -----------------------------
    
    disp(it);
    
    clf
    hold on
    
    nc=netcdf(fname0);
    zeta0=squeeze(nc{'zeta'}(it,Jstr:Jend,:));
    %k=squeeze(nc{'AKv'}(it,7,y1:y2,x1:x2));
    k=squeeze(nc{'tke'}(it,7,Jstr:Jend,:));
    close(nc)
    zeta0(h<Dcrit)=zeta0(h<Dcrit)-h(h<Dcrit);
    zeta0(h+zeta0<Dcrit*1.001)=NaN;
    ktmp=0.*k;
    ktmp(abs(zeta0)./h>0.1 | h<3.0)=1;
    k=k.*ktmp;
    k(k>1.e-2)=1;
    
    % Plot orientation
    view(-24,42)  % rear view
    
    % sea surface
    ax1=axes;
        view(-24,42)  % rear view

    f=surf(x,y,zeta0);
    f.MeshStyle='none';
    f.FaceColor='interp';
    f.AmbientStrength=0.6;
    
    % camlight headlight
    lighting gouraud
    light('Position',[-1 0 0],'Style','local')
    camlight right
    
    % tracer
    nc=netcdf(fname0);
    tpas=squeeze(nc{'tpas01'}(it,10,Jstr:Jend,:));
    nc.close();
%     z=zeta0;
%     z(tpas<0.5)=NaN;
%     p=surf(x,y,z);
%     p.FaceColor='red';
%     p.EdgeColor='none';
%     p.EdgeAlpha=0.3;
%     p.FaceAlpha=0.2;
%     p.FaceLighting='none';
    
    ax2 = axes;
        view(-24,42)  % rear view

    tpas(tpas<0.5)=NaN;
    f=surf(x,y,zeta0,tpas);
    f.MeshStyle='none';
    f.FaceColor='interp';
    f.FaceAlpha=0.3;
    f.AmbientStrength=0.6;
    colormap(ax2,hot);
        % camlight headlight
    lighting gouraud
    light('Position',[-1 0 0],'Style','local')
    camlight right

    % foam effect
    ax3=axes;
    z=zeta0;
    z(k<1)=NaN;
    p=surf(x,y,z);
    p.FaceColor='w';
    p.EdgeColor='none';
    %p.EdgeAlpha=0.7;
    p.FaceAlpha=0.3;
    p.FaceLighting='none';
    
    % Bathy isosurface or isocontours
    g=surf(x,y,-h);
    g.FaceColor=0.8/255*[230,185,99];
    g.EdgeColor='none';
    %g.EdgeAlpha=0.8;
    g.FaceAlpha=0.8;
    %g.FaceLighting='none';
    
    % map
    colormap(map);
    caxis([-3 4]);
    axis tight;
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    
    set(gca, 'Color','none', 'XColor','w', 'YColor','w', 'ZColor','w')
    set(gcf, 'Color','k')
    
    xlim([X3min X3max]);
    zlim([depth 3]);
    ylim(Ylimits);
    
    %box on
    %set(gca,'boxstyle','full')
    pbaspect([1 2.5 0.4]); % change box aspect ratio
    
    %%%%%%% WRITE VIDEO %%%%%%%%%%%%
    if make_movie
        frame=getframe(gcf);
        %frame.cdata=frame.cdata(1:xframe,1:yframe,:);
        writeVideo(v,frame);
    else
        drawnow
        pause(.2)
    end

end % time index ------------------------------

if make_movie
 close(v)
end

if makepdf,
 set(gcf,'PaperPositionMode','auto');
 %export_fig -transparent wave3D.pdf
 export_fig -transparent wave3D.png
 %print -dpng wave3D.png
end

