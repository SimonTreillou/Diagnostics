% program spec_flux
%
clear all
close all
%=====================================
%
indir='/Users/simon/Code/IB09/passolo/';
%indir='/Users/simon/Code/CONFIGS/CALMIP/IB09_dx1_ang7_spread10_init016/';
fname='rip_his.nc';
fname1='rip_his.nc';
fname2='rip_his_2D_9h_Cd0035.nc';
%
imin0=250; imax0=300;
%imin0=114; imax0=114;
%
makepdf=1;
compute=1;
wavenumber=1;
%
%=====================================
disp('')
disp(['Get energy spectrum'])

npts=1024;

figure('Units','pixels','Position',[500 500 600 500]);

set(gca,'FontSize',15)

for i=1:1;

if i==1,
 fname=fname1;
 imin=imin0; imax=imax0;
 vlevel=10;
 col='k';
 linew=2;
elseif i==2,
 fname=fname2;
 imin=imin0; imax=imax0;
 vlevel=0;
 col=0.8*[0 1 0];
 linew=3;
else
 fname=fname1;
 imin=imin0; imax=imax0;
 vlevel=0;
 col='g';
 linew=1;
end

outname=['rip_Eflux',num2str(i),'.mat'];

if compute,

fname=[indir,fname]
nc=netcdf(fname);
x=nc{'xi_rho'}(:);
h=nc{'h'}(:);
close(nc)
[M,L]=size(x);
% Initalize
%
[k,E]=get_e_spectra(fname,1,vlevel,imin,npts);
E=0*E;
k=0*k;
n=0;
%
% Perform the FFT
%
disp(['Opening : ',fname])
nc=netcdf(fname);
ntime=length(nc('time'));
close(nc)
for t=20:ntime; %1:ntime
  disp(['Time index: ',num2str(t)])
  n=n+1;
  [k,E1]=get_eflux_spectra(fname,t,vlevel,imin,imax,npts);
  E=E+E1;
end
E=E/n;

%
% Spectral integration
%
Eflux=cumsum(E1);
E1=Eflux(end)-Eflux;

%
% Save
%
%save(outname,'E','k')

else % do not compute

%load(outname,'E','k')

end

%
% Plot
%
if wavenumber
fname=[indir,fname]
 semilogx(k,smooth(E,10),'color',col,'linewidth',linew)
else
 lambda=2*pi./k;
 semilogx(lambda,smooth(E,10),'color',col,'linewidth',linew)
end
hold on

end % file name index i ---

legend('ROMS 2D', ...
       'ROMS 3D') %, ...
%       'ROMS 2D \mu = 0.006 m/s'); %, ...
    %   'ROMS 3D depth-mean V')

ylabel('Spectral Energy Flux m^{2}.s^{-3}')
line([2.e-3 1],[0 0],'color','k')
if wavenumber
 axis([4.e-3 1 -0.9e-4 0.6e-4])
 xlabel('Wavenumber m^{-1}')
 annotation('textarrow',[0.5 0.5],[0.2 0.4], ...
           'String','Injection','fontsize',15,'color','r');
 annotation('textarrow',[0.45 0.3],[0.5 0.5], ...
             'String','','fontsize',15,'color','r');
 annotation('textarrow',[0.52 0.68],[0.73 0.73], ...
            'String','','fontsize',15,'color','r');
else
 axis([6 800 -1.e-4 0.5e-4])
 xlabel('Wavelength m}')
 annotation('textarrow',[0.4 0.4],[0.2 0.4], ...
           'String','Injection','fontsize',15,'color','r');
end
hold off

if makepdf,
 set(gcf,'PaperPositionMode','auto');
 export_fig -transparent spec_flux.pdf
end


