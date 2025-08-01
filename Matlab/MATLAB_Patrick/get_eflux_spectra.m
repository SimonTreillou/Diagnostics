function [k,Eflux]=get_e_spectra(fname,tindex,vlevel,imin,imax,npts);
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

%
% Get advu,advv
%
advC2=0;
[Mu,Lu]=size(u); L=Lu; M=Mu-1; Lm=L-1; Mm=M-1; Lm2=L-2; Mm2=M-2;
advu=zeros(size(u)); advv=zeros(size(v));
if advC2,
    advu(2:Mm2,2:Lm2)=-0.25*u(2:Mm2,2:Lm2) ...
                    .*(u(2:Mm2,3:Lm)-u(2:Mm2,1:Lm2-1)) ...
                    .*(pm(2:Mm2,2:Lm2)+pm(2:Mm2,3:Lm)) ...
                    -0.0625*(v(2:Mm2  , 2:Lm2 )+v(2:Mm2  ,1:Lm2-1)+ ...
                             v(1:Mm2-1,1:Lm2-1)+v(1:Mm2-1,2:Lm2  )) ...
                    .*(u(3:Mm,2:Lm2)-u(1:Mm2-1,2:Lm2)) ...
                    .*(pn(2:Mm2,2:Lm2)+pn(2:Mm2,3:Lm));
    advv(2:Mm2,2:Lm2)=-0.25*v(2:Mm2,2:Lm2) ...
                    .*(v(3:Mm,2:Lm2)-v(1:Mm2-1,2:Lm2)) ...
                    .*(pn(2:Mm2,2:Lm2)+pn(3:Mm,2:Lm2)); % ...
            %        -0.0625*(u(2:Mm2  ,2:Lm2)+u(2:Mm2  ,3:Lm )+ ...
            %                 u(1:Mm2-1,3:Lm )+u(1:Mm2-1,2:Lm2)) ...
            %        .*(v(2:Mm2,3:Lm)-v(2:Mm2,1:Lm2-1)) ...
            %        .*(pm(2:Mm2,2:Lm2)+pm(3:Mm,2:Lm2));
else,
    advu(2:Mm2,2:Lm2)=-0.25*u(2:Mm2,2:Lm2) ...
                    .*( (u(2:Mm2,3:Lm)-u(2:Mm2,1:Lm2-1)) ...
                        -1/6*(u(2:Mm2,4:L)-3.*u(2:Mm2,3:Lm)+ ...       % 4th order term
                              3.*u(2:Mm2,2:Lm2)-u(2:Mm2,1:Lm2-1)) ...  %
                      ) ...
                    .*(pm(2:Mm2,2:Lm2)+pm(2:Mm2,3:Lm)) ...
                    -0.0625*(v(2:Mm2  ,2:Lm2  )+v(2:Mm2  ,1:Lm2-1)+ ...
                             v(1:Mm2-1,1:Lm2-1)+v(1:Mm2-1,2:Lm2  )) ...
                    .*( (u(3:Mm,2:Lm2)-u(1:Mm2-1,2:Lm2)) ...
                        -1/6*(u(4:M,2:Lm2)-3.*u(3:Mm,2:Lm2)+ ...       % 4th order term
                              3.*u(2:Mm2,2:Lm2)-u(1:Mm2-1,2:Lm2)) ...  %
                      ) ...
                     .*(pn(2:Mm2,2:Lm2)+pn(2:Mm2,3:Lm));

    advv(2:Mm2,2:Lm2)=-0.25*v(2:Mm2,2:Lm2) ...
                    .*( (v(3:Mm,2:Lm2)-v(1:Mm2-1,2:Lm2)) ...
                        -1/6*(v(4:M,2:Lm2)-3.*v(3:Mm,2:Lm2)+ ...       % 4th order term
                              3.*v(2:Mm2,2:Lm2)-v(1:Mm2-1,2:Lm2)) ...  %
                      ) ...
                     .*(pn(2:Mm2,2:Lm2)+pn(3:Mm,2:Lm2)); % ...
         %           -0.0625*(u(2:Mm2  ,2:Lm2)+u(2:Mm2  ,3:Lm )+ ...
         %                    u(1:Mm2-1,3:Lm )+u(1:Mm2-1,2:Lm2)) ...
         %           .*( (v(2:Mm2,3:Lm)-v(2:Mm2,1:Lm2-1)) ...
         %               -1/6*(v(2:Mm2,4:L)-3.*v(2:Mm2,3:Lm)+ ...       % 4th order term
         %                     3.*v(2:Mm2,2:Lm2)-v(2:Mm2,1:Lm2-1)) ...  %
         %             ) ...
         %            .*(pm(2:Mm2,2:Lm2)+pm(3:Mm,2:Lm2));
end

u=u2rho_2d(u);
v=v2rho_2d(v);
advu=u2rho_2d(advu);
advv=v2rho_2d(advv);

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
advui=squeeze(advu(:,i));
advvi=squeeze(advv(:,i));
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
windowing=1;
N=length(ui);
%
if windowing,
  ui=ui-mean(ui);
  vi=vi-mean(vi);
  advui=advui-mean(advui);
  advvi=advvi-mean(advvi);
  f=tukeywin(N,0.9); %window(N,'flattop'); %window(@tukeywin,N);
  ui=ui.*f;
  vi=vi.*f;
  advui=advui.*f;
  advvi=advvi.*f;
elseif detrending
 ui=ui-(ui(end)-ui(1))*(x-x(1))./(x(end)-x(1));
 vi=vi-(vi(end)-vi(1))*(x-x(1))./(x(end)-x(1));
 advui=advui-(advui(end)-advui(1))*(x-x(1))./(x(end)-x(1));
 advvi=advvi-(advvi(end)-advvi(1))*(x-x(1))./(x(end)-x(1));
end
%
% Get the  energy spectrum 
%
FU=fft(ui,npts);
FV=fft(vi,npts);
FADVU=fft(advui,npts);
FADVV=fft(advvi,npts);
%E1=real(FADVU.*conj(FU) +  FADVV.*conj(FV))./npts;
E1=real(FADVV.*conj(FV))./npts;
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

Eflux=E/n;

% 
%
return
