%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname = '/Users/simon/Code/IB09/testVADVC2periodok/rip_his.nc';

makepdf   = 0;                % make pdf file
meanlongs = 1;                % allow to reduce incertainty
%
%======================================================================
yindex = 350;
g = 9.81;
tempo=0;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record
%tindex=5; 

if tempo,
 tstr=1;
 tend=tindex;
else,
 tstr=tindex;
 tend=tindex;
end



for tindex=tstr:tend % ---------------------------------------------
%
% horizontal grid
 hr=squeeze(nc{'h'}(yindex,:));
 L=length(hr);
 xr=squeeze(nc{'x_rho'}(yindex,:));
 yr=squeeze(nc{'y_rho'}(yindex,:));

 dx=xr(2)-xr(1);
 dy=dx;
 xmin=min(xr);
 xmax=max(xr);
 Dcrit=nc{'Dcrit'}(:);
%
% vertical grid
 N=length(nc('s_rho'));
 N=1;
 theta_s=nc.theta_s(:); 
 theta_b=nc.theta_b(:); 
 hc=nc.hc(:); 
 zeta=squeeze(nc{'zeta'}(tindex,yindex,:));
 zr=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,10,'r',2));
 zw=squeeze(zlevs(hr,zeta,theta_s,theta_b,hc,10,'w',2));
 dzr=zr(2:end,:)-zr(1:end-1,:);         % ---> zw(2:N,:)
 zru=0.5*(zr(:,1:end-1)+zr(:,2:end));
 dzru=zru(2:end,:)-zru(1:end-1,:);      % ---> zwu(2:N,:)
 dzw=zw(2:end,:)-zw(1:end-1,:);         % ---> zr
 zwu=0.5*(zw(:,1:end-1)+zw(:,2:end));
 dzwu=zwu(2:end,:)-zwu(1:end-1,:);      % ---> zru
%
 xr2d=repmat(xr,[N 1]);
 D=hr+zeta;
 D2d=repmat(D,[N 1]);
%
 %xmin= min(xr); xmax=90;
 zmin=-max(hr); zmax=0.2*max(hr);


 % ---------------------------------------------------------------------
 % --- read/compute numerical model fields (index 1) ---
 % ---------------------------------------------------------------------
 time=nc{'scrum_time'}(tindex);

 % ... zonal velocity ...                         ---> xu,zru
 u=squeeze(nc{'u'}(tindex,N,yindex,:));
 v=squeeze(nc{'v'}(tindex,N,yindex,:));
 v=(v(2:end)+v(1:end-1))/2;

 % ... vertical velocity ...                      ---> xr,zw
 w=squeeze(nc{'w'}(tindex,:,yindex,:));
 w=(w(2:end,:)+w(1:end-1,:))/2;                   %--> w-->r
 w=squeeze(w(N,:));
 w=(w(2:end)+w(1:end-1))/2;

 % ... grad V ...
 dzru2=(dzru(N,3:end)+dzru(N,2:end-1));

 % ... advection ...
 ududx = (u(3:end)+u(2:end-1))/2 .* (u(3:end)-u(1:end-2))/(2*dx);
 vdudy = (v(3:end)+v(2:end-1))/2 .* (u(3:end)-u(1:end-2))/(2*dx);
 wdudz = (w(3:end)+w(2:end-1))/2 .* (u(3:end)-u(1:end-2))./dzru2;
    
 udvdx = (u(3:end)+u(2:end-1))/2 .* (v(3:end)-v(1:end-2))/(2*dx);
 vdvdy = (v(3:end)+v(2:end-1))/2 .* (v(3:end)-v(1:end-2))/(2*dx);
 wdvdz = (w(3:end)+w(2:end-1))/2 .* (v(3:end)-v(1:end-2))./dzru2;

 udwdx = (u(3:end)+u(2:end-1))/2 .* (w(3:end)-w(1:end-2))/(2*dx);
 vdwdy = (v(3:end)+v(2:end-1))/2 .* (w(3:end)-w(1:end-2))/(2*dx);
 wdwdz = (w(3:end)+w(2:end-1))/2 .* (w(3:end)-w(1:end-2))./dzru2;

 ADV=[ududx+vdudy+wdudz; udvdx+vdvdy+wdvdz; udwdx+vdwdy+wdwdz];

 % ... diffusion ...
 d2udx2=(u(3:end)-2*u(2:end-1)+u(1:end-2))/dx^2;
 d2udy2=(u(3:end)-2*u(2:end-1)+u(1:end-2))/dx^2;
 d2udw2=(u(3:end)-2*u(2:end-1)+u(1:end-2))./(dzru2/2).^2;

 d2vdx2=(v(3:end)-2*v(2:end-1)+v(1:end-2))/dx^2;
 d2vdy2=(v(3:end)-2*v(2:end-1)+v(1:end-2))/dx^2;
 d2vdw2=(v(3:end)-2*v(2:end-1)+v(1:end-2))./(dzru2/2).^2;

 d2wdx2=(w(3:end)-2*w(2:end-1)+w(1:end-2))/dx^2;
 d2wdy2=(w(3:end)-2*w(2:end-1)+w(1:end-2))/dx^2;
 d2wdw2=(w(3:end)-2*w(2:end-1)+w(1:end-2))./(dzru2/2).^2;

 DIF=[d2udx2+d2udy2+d2udw2; d2vdx2+d2vdy2+d2vdw2; d2wdx2+d2wdy2+d2wdw2];

 % ... viscosity ... 
 t=squeeze(nc{'AKv'}(tindex,N,yindex,:));

 % ... density ...
 rho=squeeze(nc{'rho'}(tindex,N,yindex,:));

 % ... Reynolds number ...
 Re=rho*norm(ADV)./(t.*norm(DIF));
 Re=sum(ADV.^2/2,1)./sum(DIF.^2/2,1);

 plot(xr(1:end-3)-300,Re);
 hold on

end
