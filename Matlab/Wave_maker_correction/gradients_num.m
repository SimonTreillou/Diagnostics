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
g = 9.81;
tempo=0;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record
time=nc{'scrum_time'}(:);
yr=squeeze(nc{'y_rho'}(:,1));






for yindex=1:size(yr,1) % ---------------------------------------------
%
% horizontal grid
 hr=squeeze(nc{'h'}(yindex,:));
 L=length(hr);
 xr=squeeze(nc{'x_rho'}(yindex,:));

 dx=xr(2)-xr(1);
 dy=dx;
 xmin=min(xr);
 xmax=max(xr);
 Dcrit=nc{'Dcrit'}(:);
%
% vertical grid
 N=length(nc('s_rho'));
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
 
 % ... zonal velocity ...                         ---> xu,zru
 v1=squeeze(nc{'vbar'}(tindex,yindex,:));
 if yindex<size(yr,1)
    v2=squeeze(nc{'vbar'}(tindex,yindex+1,:));
 else
    v2=squeeze(nc{'vbar'}(tindex,1,:));
 end
 dv=v2-v1;
 dvs(yindex)=norm(dv);
end
