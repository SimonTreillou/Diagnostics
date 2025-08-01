clear all
close all

DIR='/Users/simon/Code/CONFIGS/testIB09GPP/';
DIRout='../';
mfile=[DIR,'rip_avg.nc'];


load kmke

nc=netcdf(mfile);
u=nc{'u'}(:);
v=nc{'v'}(:);
zeta=nc{'zeta'}(:);
xr=nc{'x_rho'}(:);
yr=nc{'y_rho'}(:);
h=nc{'h'}(:);

x1=400; x2=240;
y1= 120; y2= 180;

xsec=linspace(x1,x2);
ysec=linspace(y1,y2);

xmin=x1-10; xmax=x2+10;
ymin=y1-10; ymax=y2+10;
sub=xr>xmin & xr<xmax & yr>ymin & yr<ymax;
ival=sum(sub,1);
jval=sum(sub,2);
imin=min(find(ival~=0));
imax=max(find(ival~=0));
jmin=min(find(jval~=0));
jmax=max(find(jval~=0));
xr=xr(jmin:jmax,imin:imax);
yr=yr(jmin:jmax,imin:imax);
h=h(jmin:jmax,imin:imax);
sub=sub(jmin:jmax,imin:imax);
zeta=zeta(jmin:jmax,imin:imax);

N=20;
theta_s = 0.;
theta_b =  0.;
hc  = 1.e16;

us=u2rho_2d(squeeze(u(end,:,:)));
us=us(jmin:jmax,imin:imax);
contourf(xr,yr,us,N); shading flat; hold on
line(xsec,ysec,'color','r')
axis([-500 -200 50 200]);

L=length(xsec);
gamma=atan((y2-y1)/(x2-x1));
usec=zeros(N,L);
for k=1:N
  ur=u2rho_2d(squeeze(u(k,:,:)));
  ur=ur(jmin:jmax,imin:imax);
  usec(k,:)=griddata(xr,yr,ur,xsec,ysec);
%  vr=v2rho_2d(squeeze(v(k,:,:)));
%  vr=vr(jmin:jmax,imin:imax);
%  vsec(k,:)=griddata(xr,yr,vr,xsec,ysec);
%  alpha=atan(vsec./usec)-gamma;
%  usec=sqrt(usec.^2+vsec.^2).*cos(alpha);
end
hsec=griddata(xr,yr,h,xsec,ysec);
zetasec=griddata(xr,yr,zeta,xsec,ysec);
zsec=squeeze(zlevs(hsec,zetasec,theta_s,theta_b,hc,N,'r',2));
L=length(hsec);
xsec=repmat(xsec,N,1);

mask=ones(size(hsec));
mask((hsec+zsec(N,:))<0.2)=NaN;
mask=repmat(mask,N,1);
usec=usec.*mask;

figure('position',[500 500 700 400])
contourf(xsec,zsec,100*usec,[-100:5:100]); hold on
[C,h1]=contour(xsec,zsec,100*usec,[-50:5:5],'k');
contour(xsec,zsec,100*usec,0,'k','linewidth',2);
clabel(C,h1,'LabelSpacing',400,'fontsize',14)
map=colormap(jet(40));
%map(1:64-31,:)=[];
colormap(map);
caxis([-100 100])
%colorbar
line(xsec,zsec(1,:),'color','k','linewidth',3)
xlabel('Cross-shore distance (m)','fontsize',15)
ylabel('Depth (m)','fontsize',15)
%title('Cross-shore current profile','fontsize',15)
set(gca,'FontSize',15)
%
axes('position',[0.62 0.2 0.25 0.3])
[C,h2]=contour(XR,YR,H,20);
set(h2,'color','k'); hold on
line(xsec,ysec,'color','r')
axis([-600 -150 -100 300])
%
set(gca,'FontSize',13)
set(gcf, 'PaperPositionMode','auto')
print -dpdf file.pdf
eval('!pdfcrop file.pdf uprof.pdf ')
eval('!mv  uprof.pdf /Users/pmarches/Documents/TEX/BISCAROSSE/FIG_PDF/.')


close(nc)


