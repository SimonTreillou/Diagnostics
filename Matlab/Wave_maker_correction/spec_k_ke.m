% program spectrum_latband_quick
%
clear all
close all
%=====================================
%
indir='/Users/simon/Code/IB09/IB09PATdxbathynewang11strat10radiatif_newWM2023/';
fname1='rip_avg_out.nc';
fname2='rip_avg.nc';%'rip_avg_2D_SC.nc';
fname3='rip_avg_2D_LC.nc';
fname1 = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname1 = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

%
nbfiles=1;
%
imin0=1; imax0=imin0;
%
makepdf=0;
compute=1;
logplot=1;
%
%=====================================
disp('')
disp(['Get energy spectrum'])

npts=1024;

figure('Units','pixels','Position',[500 500 600 500]);

set(gca,'FontSize',15)

for i=1:nbfiles;

if i==1,
 fname=fname1;
 imin=imin0; imax=imax0;
 vlevel=10;
 col='k';
 linew=2;
elseif i==2,
 fname=fname2;
 imin=imin0; imax=imax0;
 vlevel=10;
 col='b';
 linew=1;
elseif i==3,
 fname=fname3;
 %imin0=110; imax0=imin0;
 %imin=imin0+96; imax=imax0+96;
 imin=imin0; imax=imax0;
 vlevel=10;
 col='m--';
 linew=1;
else
 fname=fname1;
 imin=imin0; imax=imax0;
 vlevel=0;
 col='g';
 linew=1;
end

vlevel=1;

outname=['rip_E',num2str(i),'.mat'];

if compute,

%fname=[indir,fname]
nc=netcdf(fname);
x=nc{'xi_rho'}(:);
h=nc{'h'}(1,:);
close(nc)
[M,L]=size(x);
%
% Initialize
%
[k,E]=get_e_spectra(fname,1,vlevel,imin,npts);
E=0*E;
k=0*k;
n=0;

% Perform the FFT
%
disp(['Opening : ',fname])
nc=netcdf(fname);
ntime=length(nc('time'));
close(nc)
for t=15:ntime; %1:ntime
 disp(['Time index: ',num2str(t)])
 for vlevel=1:10;
  n=n+1;
  [k,E1]=get_e_spectra2(fname,t,vlevel,imin,imax,npts);
  E=E+E1;
 end
end
dk=k(2)-k(1);
E=E/n/dk;

%
% Save
%
save(outname,'E','k')

else % do not compute

load(outname,'E','k')

end

%if i==3
 for ismoo=1:20
  E(k>0.2)=smooth(E(k>0.2),10);
 end
%end

%
% Plot
%
if logplot,
  loglog(k,E,col,'linewidth',linew)
else
  semilogx(k,k'.^2.*smooth(E,1),col,'linewidth',linew)
end
hold on

end % file name index i ---

legend( ...
       'CROCO 3D SC', ...
       'CROCO 2D SC', ...
       'CROCO 2D LC');

axis([-Inf 3 1.e-3 Inf])

morestuff=1;
if morestuff,
if logplot,
 grey=0.7*[1 1 1];
 k53=(k.^(-5/3)); 
 k3=(k.^(-3));
 loglog(k,0.4*k53)
 text(0.006,5.e2,'k^{-5/3}','color',grey,'fontsize',15)
 loglog(k,0.8e-2*k3,'color',grey)
 text(0.4,0.2,'k^{-3}','color',grey,'fontsize',15)
 %axis([2.e-3 1 5.e-5 1.e4])
 %title('KE spectrum')
 %annotation('textarrow',[0.53 0.53],[0.4 0.65], ...
 %           'String','Injection','fontsize',15,'color','r')
else
 %axis([2.e-3 1 0 0.20])
 %title('Compensated KE spectrum')
 annotation('textarrow',[0.53 0.53],[0.2 0.3], ...
            'String','Injection','fontsize',15,'color','r')
 %grid on
end
end

hold off
xlabel('Wavenumber m^{-1}')
ylabel('Spectral density m^{3}.s^{-2}')

if makepdf,
 set(gcf,'PaperPositionMode','auto');
 if logplot
  export_fig -transparent spec_k_ke.pdf
 else
  export_fig -transparent spec_k_compke.pdf
 end
end


