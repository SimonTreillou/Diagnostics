clear all
close all

DIR='/Users/pmarches/Roms_tools/Run_rip/';
DIRout='/Users/pmarches/Documents/TEX/BISCAROSSE/FIG_PDF/';
mfile=[DIR,'BISCA_his.nc'];
gfile=[DIR,'BISCA/BISCA_grd.nc'];

tindx=150;

plot_all=0;
waves=1;
%.............................................
nc=netcdf(gfile);
xr   = nc{'x_rho'}(:);
yr   = nc{'y_rho'}(:);
pm   = nc{'pm'}(:);
pn   = nc{'pn'}(:);
h    = nc{'h'}(:);
close(nc)

[Mr,Lr]=size(xr); N=20;
Lp=Lr-1; Mp=Mr-1;
L=Lp-1; M=Mp-1;

xu=rho2u_2d(xr);
yu=rho2u_2d(yr);
pmu=rho2u_2d(pm);
pnu=rho2u_2d(pn);

theta_s = 0.;
theta_b =  0.;
hc  = 1.e16;

grav=9.81;
rho0=1000;

cmin=-1.e-3; cmax=1.e-3; nbcol=40;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];

%.......................................................
 disp('load mean fields ...')

nc=netcdf(mfile);
zeta=squeeze(nc{'zeta'}(tindx,:,:));
ubar=squeeze(nc{'ubar'}(tindx,:,:));
vbar=squeeze(nc{'vbar'}(tindx,:,:));
u=squeeze(nc{'u'}(tindx,:,:,:));
v=squeeze(nc{'v'}(tindx,:,:,:));
bostr=squeeze(nc{'bostr'}(tindx,:,:))./rho0;
if waves
 epb=squeeze(nc{'epb'}(tindx,:,:));
 epr=squeeze(nc{'epr'}(tindx,:,:));
 frq=squeeze(nc{'frq'}(tindx,:,:));
 wkx=squeeze(nc{'wkx'}(tindx,:,:));
 wke=squeeze(nc{'wke'}(tindx,:,:));
end
close(nc)

up=u-tridim(ubar,N);
vp=v-tridim(vbar,N);
up_v=rho2v_3d(u2rho_3d(up));
vp_u=rho2u_3d(v2rho_3d(vp));

ur=u2rho_3d(u);
vr=v2rho_3d(v);
ubr=u2rho_2d(ubar);
vbr=v2rho_2d(vbar);

zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w',2);
dz=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=rho2u_3d(dz);
dzv=rho2v_3d(dz);

ubot=squeeze(ur(1,:,:));
vbot=squeeze(vr(1,:,:));
bustr=bostr./sqrt(ubot.^2+vbot.^2).*ubot;
bustr=bostr./sqrt(ubot.^2+vbot.^2).*vbot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute divergence of Reynolds stresses
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 disp('-duu*dx ...')  % ---------------

uub=squeeze(sum(dzu.*up.*up)./sum(dzu));
duudx=zeros(size(ubar));
duudx(:,2:end-1)=-0.25*(uub(:,3:end)-uub(:,1:end-2)).* ...
                       (pmu(:,3:end)+pmu(:,1:end-2));

dububdx=zeros(size(ubar));
dububdx(:,2:end-1)=-0.25*(ubar(:,3:end).^2-ubar(:,1:end-2).^2).* ...
                         (pmu(:,3:end)+pmu(:,1:end-2));
ubdubdx=zeros(size(ubar));
ubdubdx(:,2:end-1)=-0.25*ubar(:,3:end).* ...
                         (ubar(:,3:end)-ubar(:,1:end-2)).* ...
                         (pmu(:,3:end)+pmu(:,1:end-2));
   
disp(max(max(abs(duudx))))
if plot_all
 figure(1); colormap(map);
 pcolor(xu,yu,duudx);shading flat;colorbar;
 axis([-450 -210 -100 300]); caxis([cmin cmax])
else
 close(1)
end

 disp('-duv*dx ...')  % ---------------

uvb=squeeze(sum(dzu.*up.*vp_u)./sum(dzu));
duvdy=zeros(size(ubar));
duvdy(2:end-1,:)=-0.25*(uvb(3:end,:)-uvb(1:end-2,:)).* ...
                       (pnu(3:end,:)+pnu(1:end-2,:));

vbu=rho2u_2d(vbr);
dubvbdy=zeros(size(ubar));
dubvbdy(2:end-1,:)=-0.25*(ubar(3:end,:).*vbu(3:end,:) ...
                         -ubar(1:end-2,:).*vbu(1:end-2,:)).* ...
                         (pnu(3:end,:)+pnu(1:end-2,:));
vbdubdy=zeros(size(ubar));
vbdubdy(2:end-1,:)=-0.25*vbu(2:end-1,:).* ...
                         (ubar(3:end,:)-ubar(1:end-2,:)).* ...
                         (pnu(3:end,:)+pnu(1:end-2,:));
disp(max(max(abs(duvdy))))
if plot_all
 figure(2); colormap(map);
 pcolor(xu,yu,duvdy);shading flat;colorbar;
 axis([-450 -210 -100 300]); caxis([cmin cmax])
end

% Kinetic energy conversion Kmke   % ---------------
%
 disp('uu*dudx ...')
uududx=zeros(size(xr));
uub_r=u2rho_2d(uub);
uududx(:,2:end-1)=-uub_r(:,2:end-1) ...
                 .*(ubar(:,2:end)-ubar(:,1:end-1)).*pm(:,2:end-1);
 disp('uv*dudy ...')
uvdudy=zeros(size(xr));
uvb_r=u2rho_2d(uvb);
uvdudy(2:end-1,:)=-uvb_r(2:end-1,:) ...
                 .*0.5.*(ubr(3:end,:)-ubr(1:end-2,:)).*pn(2:end-1,:);
 disp('vv*dvdy ...')
vvb=squeeze(sum(dzv.*vp.*vp)./sum(dzv));
vvb_r=v2rho_2d(vvb);
vvdvdy=zeros(size(xr));
vvdvdy(2:end-1,:)=-vvb_r(2:end-1,:) ...
                 .*(vbar(2:end,:)-vbar(1:end-1,:)).*pn(2:end-1,:);
 disp('uv*dvdx ...')
uvdvdx=zeros(size(xr));
uvdvdx(:,2:end-1)=-uvb_r(:,2:end-1) ...
                 .*0.5.*(vbr(:,3:end)-vbr(:,1:end-2)).*pm(:,2:end-1);

kmke=uududx+uvdudy+vvdvdy+uvdvdx;

%
% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% baroclinic conversion: -uududx ...
figure(2)
cmin=-3.e-4; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr,yr,kmke,[cmin:cint:cmax]); hold on
shading flat; colorbar;
contour(xr,yr,h,[-10:0.5:2],'color','r')
quiver(xr,yr,ubr,vbr,1.5);
axis([-450 -210 -100 300]);
caxis([cmin cmax])
title('Conversion to baroclinic KE','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');

% baroclinic advection: -duudx ...
figure(3)
cmin=-3.e-3; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xu,yu,duudx+duvdy,[cmin:cint:cmax]); hold on
shading flat; colorbar;
contour(xr,yr,h,[-10:0.5:2],'color','r')
quiver(xr,yr,ubr,vbr,1.5);
axis([-450 -210 -100 300]);
caxis([cmin cmax])
title('Baroclinic advection','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');
%print -dpdf file.pdf
%eval(['!pdfcrop file.pdf ',DIRout,'duudx_vert.pdf'])

% Reynolds Stress: uu
figure(4)
cmin=-7.e-2; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xu,yu,uub,[cmin:cint:cmax]); hold on
shading flat; colorbar;
contour(xr,yr,h,[-10:1:2],'color','r')
quiver(xr,yr,ubr,vbr,1.5);
axis([-450 -210 -100 300]);
caxis([cmin cmax])
title('Reynolds stress','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');
%print -dpdf file.pdf
%eval(['!pdfcrop file.pdf ',DIRout,'uu_vert.pdf'])

% Breaking
if waves
kw=sqrt(wkx.^2+wke.^2);
kw=wkx;
brk=kw./frq.*(0.5*epb+epr)./(h+zeta+1.5239977172591646);
figure(5)
cmin=-3.e-3; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr,yr,brk,[cmin:cint:cmax]); hold on
shading flat; colorbar;
contour(xr,yr,h,[-10:1:2],'color','r')
quiver(xr,yr,ubr,vbr,2);
axis([-450 -210 -100 300]);
caxis([cmin cmax])
title('Breaking accelaration','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');
%print -dpdf file.pdf
%eval(['!pdfcrop file.pdf ',DIRout,'uu_vert.pdf'])
end

% Bottom drag
if waves
bostr=-bustr./(h+zeta+1.5239977172591646);

figure(6)
cmin=-3.e-3; cmax=3.e-3; nbcol=20;
cint=(cmax-cmin)/nbcol;
map=colormap(jet(nbcol));
map(nbcol/2  ,:)=[1 1 1];
map(nbcol/2+1,:)=[1 1 1];
colormap(map);
contourf(xr,yr,bostr,[cmin:cint:cmax]); hold on
shading flat; colorbar;
contour(xr,yr,h,[-10:1:2],'color','r')
quiver(xr,yr,ubr,vbr,2);
axis([-450 -210 -100 300]);
caxis([cmin cmax])
title('Bottom friction','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');
%print -dpdf file.pdf
%eval(['!pdfcrop file.pdf ',DIRout,'uu_vert.pdf'])
end

% Barotropic advection: -dububdx ...
figure(7)
cmin=-3.e-3; cmax=3.e-3; nbcol=20;
colormap(map);
contourf(xu,yu,ubdubdx+vbdubdy,[cmin:cint:cmax]); hold on
shading flat; colorbar;
contour(xr,yr,h,[-10:1:2],'color','r')
quiver(xr,yr,ubr,vbr,2);
axis([-450 -210 -100 300]);
caxis([cmin cmax])
title('Barotropic advection -dububdx')
%print -dpdf file.pdf
%eval(['!pdfcrop file.pdf ',DIR,'dububdx.pdf'])

