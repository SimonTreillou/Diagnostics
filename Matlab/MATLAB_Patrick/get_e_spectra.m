function [k,E]=get_e_spectra(fname,tindex,vlevel,i,npts);
%
% Get the energy spectrum for a constant latitude.
%
rmpts=1;
nc=netcdf(fname);
pn=squeeze(nc{'pn'}(1+rmpts:end-rmpts,i));
  
close(nc)
if vlevel==0
  u=u2rho_2d(get_hslice(fname,fname,'ubar',tindex,vlevel,'u'));
  v=v2rho_2d(get_hslice(fname,fname,'vbar',tindex,vlevel,'v'));
else
  u=u2rho_2d(get_hslice(fname,fname,'u',tindex,vlevel,'u'));
  v=v2rho_2d(get_hslice(fname,fname,'v',tindex,vlevel,'v'));
end
u=squeeze(u(1+rmpts:end-rmpts,i));
v=squeeze(v(1+rmpts:end-rmpts,i));

%h = 1/4*ones(2);
%u= filter2(h,u);
%v= filter2(h,v);

%
% Get distances
%
dx=1./pn;
x=0.*dx;
dx_v=0.5*(dx(1:end-1)+dx(2:end));
for i=1:length(dx_v)
  x(i+1)=dx_v(1)+x(i);
end
%
% correct series for 1D nonperiodic data
% (limited area grid) to avoid aliasing ...
%
detrending=0;
windowing=0;
N=length(u);
%
u(isnan(u))=0;
v(isnan(v))=0;
if windowing,
  u=u-mean(u);
  v=v-mean(v);
  f=window(@tukeywin,N);
  u=u.*f;
  v=v.*f;
elseif detrending
 u=u-(u(end)-u(1))*(x-x(1))./(x(end)-x(1));
 v=v-(v(end)-v(1))*(x-x(1))./(x(end)-x(1));
end
%
% Get the  energy spectrum 
%
FU=fft(u,npts);
FV=fft(v,npts);
E=real(FU.*conj(FU) +  FV.*conj(FV))./npts;
k=2*pi*mean(pn)*(0:npts/2)/npts;
imin=max(find(k<2/max(x)));
disp(imin)
imax=length(k);
E=E(imin:imax);
k=k(imin:imax);
%
% 
%
return
