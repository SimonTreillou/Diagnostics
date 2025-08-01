clear all
close all

DIR='/Users/simon/Code/CONFIGS/testIB09GPP/';
DIRout='./';

mfile=[DIR,'rip_avg.nc'];
gfile=mfile; %[DIR,'rip_grd.nc'];

printpdf=0;

tstr=1; tend=20;

%.............................................
nc=netcdf(gfile);
xr   = nc{'x_rho'}(:);
yr   = nc{'y_rho'}(:);
pm   = nc{'pm'}(:);
pn   = nc{'pn'}(:);
h    = nc{'h'}(:);
close(nc)

[Mr,Lr]=size(xr);
Lp=Lr-1; Mp=Mr-1;
L=Lp-1; M=Mp-1;

N=10;
theta_s = 0.;
theta_b =  0.;
hc  = 1.e16;

grav=9.81;
rho0=1000;

%.......................................................
 disp('load mean fields ...')

nc=netcdf(mfile);
mz=squeeze(mean(nc{'zeta'}(tstr:tend,:,:)));
zw=squeeze(zlevs(h,mz,theta_s,theta_b,hc,N,'w',2));
zr=squeeze(zlevs(h,mz,theta_s,theta_b,hc,N,'r',2));

mus=squeeze(mean(nc{'u'}(tstr:tend,N,:,:)));
mvs=squeeze(mean(nc{'v'}(tstr:tend,N,:,:)));
mws=squeeze(mean(nc{'w'}(tstr:tend,N,:,:)));
mur=u2rho_2d(mus);
mvr=v2rho_2d(mvs);
x1=min(min(xr)); x2=max(max(xr));
y1=min(min(yr)); y2=max(max(yr));
[xr2,yr2]=meshgrid([x1:10:x2],[y1:10:y2]);
mur2=griddata(xr,yr,mur,xr2,yr2);
mvr2=griddata(xr,yr,mvr,xr2,yr2);
mwr2=griddata(xr,yr,mws, xr2,yr2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute Reynolds stresses and conversion terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

muududx=zeros(size(h));
mvvdvdy=zeros(size(h));
mvvdvdy=zeros(size(h));
muvdvdx=zeros(size(h));
muwdudz=zeros(size(h));
mvwdvdz=zeros(size(h));
H=zeros(size(h));

for k=1:N;

 disp([' --- Level ',num2str(k)])

mu=squeeze(mean(nc{'u'}(tstr:tend,k,:,:)));
mv=squeeze(mean(nc{'v'}(tstr:tend,k,:,:)));
mw=squeeze(mean(nc{'w'}(tstr:tend,k,:,:)));
mur=u2rho_2d(mu);
mvr=v2rho_2d(mv);

%................................................
 disp('uu*dudx ...')

du=zeros(size(xr));
uu=zeros(size(xr));
uur=zeros(tend-tstr+1,Mr,Lr);
%uu=u2rho_2d(squeeze(mean((nc{'u'}(tstr:tend,k,:,:)).^2)));
for it=tstr:tend;
  uur(it,:,:)=u2rho_2d((squeeze(nc{'u'}(it,k,:,:))-mu).^2);
end
uu=squeeze(mean(uur,1));
du(2:Mp,2:Lp)=mu(2:Mp,2:Lp)-mu(2:Mp,1:Lp-1);
uududx(k,:,:)=uu.*du.*pm;

clear du;
disp(max(max(max(abs(uududx(k,:,:))))))
%figure; pcolor(uududx(20:end-10,15:end));shading flat;colorbar;
%return
%................................................
 disp('uv*dudy ...')

du=zeros(size(xr));
uv=zeros(size(xr));
mu_v=zeros(size(mv));
uvr=zeros(tend-tstr+1,Mr,Lr);
for it=tstr:tend;
  uvr(it,:,:)=u2rho_2d(squeeze(nc{'u'}(it,k,:,:))-mu).* ...
              v2rho_2d(squeeze(nc{'v'}(it,k,:,:))-mv);
end
uv=squeeze(mean(uvr,1));
mu_v(1:Mp,2:Lp)=0.25*(mu(2:Mr,2:Lp)+mu(1:Mr-1,2:Lp)+ ...
                      mu(2:Mr,1:Lp-1)+mu(1:Mr-1,1:Lp-1));
du(2:Mp,2:Lp)=mu_v(2:Mp,2:Lp)-mu_v(1:Mp-1,2:Lp);
uvdudy(k,:,:)=uv.*du.*pn;

clear du mu_v;
disp(max(max(max(abs(uvdudy(k,:,:))))))
%figure; pcolor(uvdudy(20:end-10,15:end));shading flat;colorbar;
%return
%....................................................
 disp('vv*dvdy ...')

dv=zeros(size(xr));
vv=zeros(size(xr));
vvr=zeros(tend-tstr+1,Mr,Lr);
%vv=v2rho_2d(squeeze(mean((nc{'v'}(tstr:tend,k,:,:)).^2)));
for it=tstr:tend;
  vvr(it,:,:)=v2rho_2d((squeeze(nc{'v'}(it,k,:,:))-mv).^2);
end
vv=squeeze(mean(vvr,1));
dv(2:Mp,2:Lp)=mv(2:Mp,2:Lp)-mv(1:Mp-1,2:Lp);
vvdvdy(k,:,:)=vv.*dv.*pn;

clear dv;
disp(max(max(max(abs(vvdvdy(k,:,:))))))
%figure; pcolor(vvdvdy(20:end-10,15:end));shading flat;colorbar;
%return

%........................................................
 disp('uv*dvdx ...')

dv=zeros(size(xr));
mv_u=zeros(size(mu));
mv_u(2:Mp,1:Lp)=0.25*(mv(2:Mp,2:Lr)+mv(1:Mp-1,2:Lr)+ ...
                      mv(2:Mp,1:Lr-1)+mv(1:Mp-1,1:Lr-1));
dv(2:Mp,2:Lp)=mv_u(2:Mp,2:Lp)-mv_u(2:Mp,1:Lp-1);
uvdvdx(k,:,:)=uv.*dv.*pm;

clear uv du mu_v mv_u uvr;
disp(max(max(max(abs(uvdvdx(k,:,:))))))
%figure; pcolor(uvdvdx(20:end-10,15:end));shading flat;colorbar;
%return

%........................................................
 disp('uw*dudz ...')

uwr=zeros(tend-tstr+1,Mr,Lr);
for it=tstr:tend;
  uwr(it,:,:)=u2rho_2d(squeeze(nc{'u'}(it,k,:,:))-mu).* ...
                      (squeeze(nc{'w'}(it,k,:,:))-mw);
end
uw=squeeze(mean(uwr,1));
kp=min(k+1,N);
km=max(k-1,1);
mu_kp=squeeze(mean(nc{'u'}(tstr:tend,kp,:,:)));
mu_km=squeeze(mean(nc{'u'}(tstr:tend,km,:,:)));
mur_km=u2rho_2d(mu_km);
du=mur-mur_km;
dz=squeeze(zr(kp,:,:)-zr(km,:,:));
uwdudz(k,:,:)=uw.*du./dz;

clear uw du mur_km mur_kp uwr;
disp(max(max(max(abs(uwdudz(k,:,:))))))
%figure; pcolor(uwdudz(20:end-10,15:end));shading flat;colorbar;
%return

%........................................................
 disp('vw*dvdz ...')

vwr=zeros(tend-tstr+1,Mr,Lr);
for it=tstr:tend;
  vwr(it,:,:)=v2rho_2d(squeeze(nc{'v'}(it,k,:,:))-mv).* ...
                      (squeeze(nc{'w'}(it,k,:,:))-mw);
end
vw=squeeze(mean(vwr,1));
kp=min(k+1,N);
km=max(k-1,1);
mv_kp=squeeze(mean(nc{'v'}(tstr:tend,kp,:,:)));
mv_km=squeeze(mean(nc{'v'}(tstr:tend,km,:,:)));
mvr_km=v2rho_2d(mv_km);
dv=mvr-mvr_km;
vwdvdz(k,:,:)=vw.*dv./dz;

clear vw dv mvr_km mvr_kp vwr;
disp(max(max(max(abs(vwdvdz(k,:,:))))))
%figure; pcolor(vwdvdz(20:end-10,15:end));shading flat;colorbar;
%return
%........................................................

dz=squeeze(zw(k+1,:,:)-zw(k,:,:));
muududx=muududx+squeeze(uududx(k,:,:)).*dz;
mvvdvdy=mvvdvdy+squeeze(vvdvdy(k,:,:)).*dz;
mvvdvdy=mvvdvdy+squeeze(vvdvdy(k,:,:)).*dz;
muvdvdx=muvdvdx+squeeze(uvdvdx(k,:,:)).*dz;
muwdudz=muwdudz+squeeze(uwdudz(k,:,:)).*dz;
mvwdvdz=mvwdvdz+squeeze(vwdvdz(k,:,:)).*dz;
H=H+dz;

end % k loop

muududx=muududx./H;
mvvdvdy=mvvdvdy./H;
mvvdvdy=mvvdvdy./H;
muvdvdx=muvdvdx./H;
muwdudz=muwdudz./H;
mvwdvdz=mvwdvdz./H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT KMKE
%
 disp('Km->Ke ...');

kmke=-(muududx+mvvdvdy+muwdudz+mvvdvdy+muvdvdx+mvwdvdz);
%kmke=-(muududx+mvvdvdy+mvvdvdy+muvdvdx);

%vort=psi2rho(vorticity(mu,mv,pm,pn));
%vort=psi2rho(shear(mu,mv,pm,pn));

maxkmke=max(max(abs(kmke)))
kmke=1.e4*kmke;

figure('position',[0 0 400 600]);
%
cmin=-80; cmax=-cmin; nbcol=20;
cint=(cmax-cmin)/nbcol;
AdvancedColormap('vbgswoyrm',nbcol)
contourf(xr,yr,kmke,[cmin:cint:cmax]); hold on
contour(xr,yr,kmke,[cmin:cint:cmax]);
shading flat; colorbar;
%hh=contour(xr,yr,vort,5,'b'); hold on
%quiver(xr2,yr2,mur2,mvr2,1.5,'color','k');
%axis([-450 -190 -25 220]);
caxis([cmin cmax])
title('KmKe [cm^2/s^3]','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');
%
if printpdf
 export_fig -transparent file.pdf
 eval(['!mv file.pdf ',DIRout,'kmke.pdf'])
end

%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EKE
%
figure('position',[0 0 400 600]);
cmin=0; cmax=1500; nbcol=20;
cint=(cmax-cmin)/nbcol;
AdvancedColormap('worm',nbcol);
eke=0.5e4*(uu+vv);
contourf(xr,yr,eke,[0:cint:cmax]); hold on
contour(xr,yr,eke,[0:cint:cmax]);
shading flat; colorbar;
%hh=contour(xr,yr,vort,10); hold on
quiver(xr2,yr2,mur2,mvr2,1.5,'color','k');
%axis([-450 -190 -25 225]);
caxis([cmin cmax])
title('EKE [cm^2/s^2]','FontSize',15)
xlabel('Cross-shore distance [m]','FontSize',15)
ylabel('Alongshore distance [m]','FontSize',15)
set(gca,'FontSize',15)
set(gcf,'PaperPositionMode','auto');
%
if printpdf
 export_fig -transparent file.pdf
 eval(['!mv file.pdf ',DIRout,'eke.pdf'])
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mub=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
mvb=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
mub=u2rho_2d(mub);
mvb=v2rho_2d(mvb);

ub1=u2rho_2d(squeeze(nc{'ubar'}(201,:,:)))-mub;
vb1=v2rho_2d(squeeze(nc{'vbar'}(201,:,:)))-mvb;
ub2=u2rho_2d(squeeze(nc{'ubar'}(205,:,:)))-mub;
vb2=v2rho_2d(squeeze(nc{'vbar'}(205,:,:)))-mvb;
ub3=u2rho_2d(squeeze(nc{'ubar'}(209,:,:)))-mub;
vb3=v2rho_2d(squeeze(nc{'vbar'}(209,:,:)))-mvb;
ub4=u2rho_2d(squeeze(nc{'ubar'}(213,:,:)))-mub;
vb4=v2rho_2d(squeeze(nc{'vbar'}(213,:,:)))-mvb;

figure
subplot(2,2,1)
quiver(xr,yr,ub1,vb1,5); hold on
quiver(xr,yr,mur,mvr,'r'); hold off
axis([-500 -190 -100 350]);
axis([-800 -100 -250 500]);

subplot(2,2,2)
quiver(xr,yr,ub2,vb2,5); hold on
quiver(xr,yr,mur,mvr,'r'); hold off
axis([-500 -190 -100 350]);
axis([-800 -100 -250 500]);

subplot(2,2,3)
quiver(xr,yr,ub3,vb3,5);hold on
quiver(xr,yr,mur,mvr,'r'); hold off
axis([-500 -190 -100 350]);
axis([-800 -100 -250 500]);

subplot(2,2,4)
quiver(xr,yr,ub4,vb4,5);hold on
quiver(xr,yr,mur,mvr,'r'); hold off
axis([-500 -190 -100 350]);
axis([-800 -100 -250 500]);


close(nc)




