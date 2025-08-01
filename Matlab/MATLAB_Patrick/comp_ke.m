%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = '/Users/simon/Code/CONFIGS/testIB09GPP/rip_avg.nc';  % roms file name
%fname     = 'rip_his_3D_9h.nc';
%fname     = 'rip_his_2D_9h_Cd004.nc';
%fname     = 'rip_his_mu003.nc';

model3D   = 1;

makepdf   = 0;             % make pdf file
Dcrit=0.2;
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);

dx=3;

tstr=30;
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1; end;

h=nc{'h'}(:,:);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:,:)-xl+93;
y=nc{'y_rho'}(:,:);
pm=nc{'pm'}(:,:);
pn=nc{'pn'}(:,:);
N=length(nc('s_rho'));

if model3D,
  ui=squeeze(nc{'u'}(1,N,:,:));
  vi=squeeze(nc{'v'}(1,N,:,:));
  u=squeeze(nc{'u'}(tstr:tend,N,:,:));
  v=squeeze(nc{'v'}(tstr:tend,N,:,:));
  ub=squeeze(nc{'u'}(tstr:tend,1,:,:));
  vb=squeeze(nc{'v'}(tstr:tend,1,:,:));
  ubar=squeeze(nc{'ubar'}(tstr:tend,:,:));
  vbar=squeeze(nc{'vbar'}(tstr:tend,:,:));
else
  ui=squeeze(nc{'ubar'}(1,:,:));
  vi=squeeze(nc{'vbar'}(1,:,:));
  mu=squeeze(mean(nc{'ubar'}(tstr:tend,:,:)));
  mv=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
  u=nc{'ubar'}(tstr:tend,:,:);
  v=nc{'vbar'}(tstr:tend,:,:);
end

eb =squeeze(nc{'epb'}(tend,:,:));
er =squeeze(nc{'epr'}(tend,:,:));
wke=squeeze(nc{'wke'}(tend,:,:));
frq=squeeze(nc{'frq'}(tend,:,:));
ar=0.5;
eb=(1-ar)*eb+er;

ui=u2rho_2d(ui);
vi=v2rho_2d(vi);
u=u2rho_3d(u);
v=v2rho_3d(v);
ubar=u2rho_3d(ubar);
vbar=v2rho_3d(vbar);
mu=squeeze(mean(u));
mv=squeeze(mean(v));
if model3D,
 ub=u2rho_3d(ub);
 vb=v2rho_3d(vb);
 mub=squeeze(mean(ub));
 mvb=squeeze(mean(vb));
 mubar=squeeze(mean(ubar));
 mvbar=squeeze(mean(vbar));
end
mx=mean(x);
mh=mean(h);
lt=size(u,1);
[M L]=size(mu);

ke=u.^2+v.^2;
kem=mean(ke);
mkei=squeeze(ui.^2+vi.^2);
mke=mu.^2+mv.^2;

up=u-reshape(repmat(mu,lt,1),lt,M,L);
vp=v-reshape(repmat(mv,lt,1),lt,M,L);;
meke=squeeze(mean(up.^2+vp.^2));

akem=mean(mean(kem));
akei=mean(mean(mkei));
amke=mean(mean(mke));
ameke=mean(mean(meke));
propKE=100*ameke/amke;

disp(['Mean KE         : ',num2str(akem)])
disp(['Initial KE      : ',num2str(akei)])
disp(['KE of mean flow : ',num2str(amke)])
disp(['Eddy KE         : ',num2str(ameke)])
disp(['% Eddy KE       : ',num2str(propKE)])

aui=sqrt(akei);
amu=sqrt(amke);
aup=sqrt(ameke);
proprms=100*aup/amu;

disp(' ...  ')
disp(['mean flow       : ',num2str(amu)])
disp(['Initial rms     : ',num2str(aui)])
disp(['Eddy rms        : ',num2str(aup)])
disp(['% Eddy rms      : ',num2str(proprms)])

clear u ke mke ui vi 

%================================================================
% -------- Bstr+Advx+Break=0 ---------
%

% Eddy mixing ---
Rstr=0.*mv;
advx=0.*Rstr;
dvdx=0.*mv;
grad=0.*mv;
diss=0.*mv;
cff=0.5/dx;
for k=1:N
  disp([num2str(k)])
  u=squeeze(nc{'u'}(tstr:tend,k,:,:));
  v=squeeze(nc{'v'}(tstr:tend,k,:,:));
  u=u2rho_3d(u);
  v=v2rho_3d(v);
  mu=squeeze(mean(u));
  mv=squeeze(mean(v));
  up=u-reshape(repmat(mu,lt,1),lt,M,L);
  vp=v-reshape(repmat(mv,lt,1),lt,M,L);;
  Rstr=squeeze(mean(up.*vp));
  advx(:,2:end-1)=advx(:,2:end-1) ...
                  - 0.5/dx*(Rstr(:,3:end)-Rstr(:,1:end-2));

  clear Rstr u mu mv up vp
  for it=1:size(v,1);
   dvdx(:,2:end-1)=squeeze(cff*(v(it,:,3:end)-v(it,:,1:end-2)));
   nusmag=0.03*dx^2*abs(dvdx);
   grad=nusmag.*dvdx;                   
   diss(:,2:end-1)=diss(:,2:end-1)+cff*(grad(:,3:end)-grad(:,1:end-2));
  end
end
advx=advx./N;
diss=diss./(N*(tend-tstr+1));
madvx=mean(advx);
mdiss=mean(diss);
%
% Bottom friction ---
cd=3.e-3; 
if model3D,
 bvstr=-cd*squeeze(mean(sqrt(ub.^2+vb.^2).*vb))./h;
else
 bvstr=-cd*squeeze(mean(sqrt(u.^2+v.^2).*v))./h;
end
cff=4./h; cff(cff>1)=1;
bvstr=cff.*bvstr;
mbvstr=mean(bvstr);
%
% Wave breaking ---
brk=eb.*wke./frq./h;
mbrk=mean(brk);
%
% subgrid-scale mixing
%if model3D
%  mv0=mvbar;
%else
%  mv0=mv;
%end
%dvdx=0.*mv;
%grad=0.*mv;
%diss=0.*mv;
%dvdx(:,2:end-1)=0.5/dx*(mv0(:,3:end)-mv0(:,1:end-2));
%nusmag=0.03*dx^2*abs(dvdx);
%grad(:,2:end-1)=nusmag(:,2:end-1).*(0.5/dx).*(mv0(:,3:end)-mv0(:,1:end-2));
%diss(:,2:end-1)=0.5/dx*(grad(:,3:end)-grad(:,1:end-2));
%mdiss=mean(diss);
%
% Residual ---
res=madvx+mbvstr+mbrk+mdiss;
%
figure
plot(mx,madvx,mx,mbvstr,mx,mbrk,mx,mdiss,mx,res);
legend('Eddy mixing','Bottom friction','Wave breaking','Subgrid mixing','Residual', ...
       'location','northwest')
axis([-150 -20 -5e-3 3.5e-3])
grid on
export_fig -transparent Advx.pdf

return
%===================================================================
%---------------------
cd=3.e-3; mu0=3.e-3;
if model3D,
 bvstr=squeeze(mean(sqrt(ub.^2+vb.^2).*vb));
 mbvstr=sqrt(mub.^2+mvb.^2).*mv;
 mbvstr2=-mu0*mean(mvb);
else
 bvstr=squeeze(mean(sqrt(u.^2+v.^2).*v));
 mbvstr=sqrt(mu.^2+mv.^2).*mv;
 mbvstr2=-mu0*mean(mv);
end
bvstr=-cd*mean(bvstr);
mbvstr=-cd*mean(mbvstr);
%figure; plot(mx,bvstr,mx,mbvstr,mx,mbvstr2);
%legend('quad drag','mean quad drag','lin drag','location','northwest')
%export_fig -transparent bvstr.pdf

% -------- KX ---------
dx=3;
Rstr=squeeze(mean(up.*vp));
%Rstr=filter_pp(x,y,Rstr,50.);
dvdx=0.*Rstr;
dvdx(:,2:end-1)=0.5/dx*(mv(:,3:end)-mv(:,1:end-2));
%dvdx=filter_pp(x,y,dvdx,50.);
mdvdx=mean(dvdx);
mRstr=mean(Rstr);
mx=mean(x);
Kx=-mRstr./mdvdx;
%
Kx=-mean(Rstr./dvdx);
%Kx(abs(mdvdx)<0.005 & abs(mRstr)>0.0001)=NaN;
%Kx=naninterp(Kx);
% --- Plot
figure
mdv=mdvdx/max(abs(mdvdx));
mr=mRstr/max(abs(mRstr));
plot(mx,smooth(Kx,3),'r',mx,-mean(mv),'k',mx,100*mRstr,'b',mx,(1/13)^2*(mx.^2).*abs(mdvdx),'g');
legend('Eddy viscosity','Mean longshore current','location','northwest')
axis([-140 -20 -2 2])
grid on
export_fig -transparent Kx.pdf

% ---------- Lx: mixing length ------------------
%uscale=squeeze(mean(std(up)))';
%Lx=Kx./uscale.;
%
%Lx=sqrt(Kx./mdvdx);
%
Lx=sqrt(abs(-mRstr))./abs(mdvdx);
%
k=mean(0.5*(up.^2+vp.^2));
mk=squeeze(mean(k))';
Lx=Kx./sqrt(mk);
figure; plot(mx,Lx); axis([-150 0 0 10])
%
mh=mean(h);
ic=max(find(abs(mh)==min(abs(mh)))); 
xc=0; %mx(ic);
C=1/13;
Lx=Lx./(C*abs(mx-xc));
figure
plot(mx,smooth(Lx,3),mx,mean(h),mx,-mean(mv));
legend(['Mixing Length / x*',num2str(C)],'location','northwest')
axis([-150 0 0 10])
grid on
export_fig -transparent Lx.pdf




return


