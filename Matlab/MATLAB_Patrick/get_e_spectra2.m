function [k,E]=get_e_spectra(fname,tindex,vlevel,imin,imax,npts);
%
%=================================================================
% Get the flux energy spectrum
%=================================================================
%

%
% Get grid
%
nc=netcdf(fname);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
close(nc)

%
% Get u,v
%
if vlevel==0
  u=get_hslice(fname,fname,'ubar',tindex,vlevel,'u');
  v=get_hslice(fname,fname,'vbar',tindex,vlevel,'v');
else
  u=get_hslice(fname,fname,'u',tindex,vlevel,'u');
  v=get_hslice(fname,fname,'v',tindex,vlevel,'v');
end

u=u2rho_2d(u);
v=v2rho_2d(v);

%
% Initalize
%
[k,E]=get_e_spectra(fname,1,vlevel,imin,npts);
E=0*E;
k=0*k;
n=0;


for i=imin:imax; %  ------------- i loop ----------------

ui=squeeze(u(:,i));
vi=squeeze(v(:,i));
%
%
% Get distances
%
pni=squeeze(pn(:,i));
dx=1./pni;
x=0.*dx;
dx_v=0.5*(dx(1:end-1)+dx(2:end));
for j=1:length(dx_v)
  x(j+1)=dx_v(1)+x(j);
end
%
% correct series for 1D nonperiodic data
% (limited area grid) to avoid aliasing ...
%
detrending=0;
windowing=0;
N=length(ui);
%
if windowing,
  ui=ui-mean(ui);
  vi=vi-mean(vi);
  f=hanning(N);
  ui=ui.*f;
  vi=vi.*f;
elseif detrending
 ui=ui-(ui(end)-ui(1))*(x-x(1))./(x(end)-x(1));
 vi=vi-(vi(end)-vi(1))*(x-x(1))./(x(end)-x(1));
end
%
% Get the  energy spectrum 
%
FU=fft(ui,npts);
FV=fft(vi,npts);
E1=real(FU.*conj(FU) +  FV.*conj(FV))./npts;
%E1=real(FV.*conj(FV))./npts;
k=2*pi*mean(pni)*(0:npts/2)/npts;
jmin=max(find(k<2/max(x)));
jmax=length(k);
E1=2*E1(jmin:jmax);
k=k(jmin:jmax);

%
% Spatial integration
%
n=n+1;
E=E+E1;

end %  -------------end i loop ----------------

E=E/n;

% 
%
return
