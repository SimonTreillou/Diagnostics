%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname1     = '/Users/simon/Code/CONFIGS/testIB09GPP/rip_avg.nc';       % roms file name
fname1     = 'rip_his_3D_9h.nc';
fname2     = 'rip_his_2D_9h.nc'; 
fname3     = 'rip_his_mu006.nc';

makepdf   = 1;                    % make pdf file
smoo      = 1;
freq      = 0;
%
%======================================================================

figure('Units','pixels','Position',[500 500 600 500]);
set(gca,'FontSize',15)

%=================================================================

for ifile=1:1  %------------------------------- loop on ROMS files

if ifile==1,
 fname=fname1;
 col='k';
 model3D   = 1;
 smoo=1;
 linew=2;
else
 fname=fname3;
 col=0.7*[1 1 1];
 model3D   = 0;
 smoo=1;
 linew=2;
end

%=================================================================
 load('DATA/V_GPP2.mat')
 %specRad(specT<3)=NaN;
 if freq,
  specF=1./specT/60.;
  loglog(specF,smooth(specRad,1),'color',0.7*[1 1 1],'linewidth',3);
 else
  %semilogy(specT,smooth(specRad,1),'color',0.75*[1 1 1],'linewidth',2);
  loglog(specT,smooth(specRad,1),'color',0.8*[1 1 1],'linewidth',3);
 end
%=================================================================

hold on

if ifile==1,
 nc=netcdf(fname);

 tstr=60;
 tend=length(nc{'scrum_time'}(:));
 if tend<tstr; tstr=1;  end;

 eb=squeeze(mean(squeeze(nc{'epb'}(tend,:,:)),1));  % !<-- Locate breaking point
 er=squeeze(mean(squeeze(nc{'epr'}(tend,:,:)),1));
 ar=0.5; eb=(1-ar)*eb+er;
 istr0=find(eb==max(eb));
 disp(['Breaking point: ',num2str(istr0)])

 istr=112; %istr0;
 iend=istr+6;

 jstr=1; jend=256;

 h=nc{'h'}(jstr:jend,istr:iend);
 xl=nc{'xl'}(:);
 x=nc{'x_rho'}(jstr:jend,istr:iend)-xl+93;
 y=nc{'y_rho'}(jstr:jend,istr:iend);
 eb=nc{'epb'}(tstr,jstr:jend,istr:iend);
 N=length(nc('s_rho'));

 if model3D,
  u=squeeze(nc{'u'}(tstr:tend,N,jstr:jend,istr:iend));
  v=squeeze(nc{'v'}(tstr:tend,N,jstr:jend,istr:iend));
 else
  u=nc{'ubar'}(tstr:tend,jstr:jend,istr:iend);
  v=nc{'vbar'}(tstr:tend,jstr:jend,istr:iend);
 end
 
 [N M L]=size(v);
 vm=squeeze(v(:,jend-jstr+1,iend-istr+1));
 dtm=60;
 tm=[0:dtm:length(vm)*dtm-dtm]';
 vv=vm; t=tm; dt=dtm;

else

 load('DATA/GPP2014_currents.mat')
 x0=75;
 xV=-X+x0;
 %
 load('DATA/GPP2014_currents_2.mat')
 v=CRad(1:360,10,:);
 [N M L]=size(v);
 t=86400*(datef(1:360)-datef(1));
 dt=t(2)-t(1);

end

n = 2^(ceil(log2(N))); 
Xm=zeros(n,1);


for j=1:M; for i=1:L; % ---------- j i

%
% Get detrended time series
%
vv=squeeze(v(:,j,i));
vv=vv-mean(vv);
x=detrend(vv,'linear');
%
% Windowing
%
x=x.*tukeywin(length(x),0.05);

%
% --- Compute Freq spectrum ---
%
%% Time specifications:
Fs = 1/dt;                    % samples per second
dt = 1/Fs;                    % seconds per sample
N = size(x,1);                % signal size
n = 2^(ceil(log2(N)));        % fft size %1024

%% Fourier Transform:
X = fftshift(fft(x,n));
Xm=Xm+abs(X)/N;

end; end; % ------------------- j i


Xm=Xm./(M*L);
if smoo,
  Xm=smooth(Xm,5);
end

%% Frequency specifications:
dF = Fs/n;                      % hertz
f = -Fs/2:dF:Fs/2-dF;           % hertz
p = (1./f)/60;                  % mn

if ifile==1
  save spec_f_ke.mat p Xm
end

%% Plot spectrum: -------------------

if freq,
 loglog(f,Xm,col,'linewidth',linew);
else
 %semilogy(p,Xm,col,'linewidth',linew); hold on
 loglog(p,Xm,col,'linewidth',linew);
end

end % ---------------------------- end loop ifile


% FINALIZE PLOT

annotation('textarrow',[11.5 11.5]./20,[0.2 0.4], ...
           'String','Injection','fontsize',15,'color','r')

hold off
%grid on
legend('VIDEO','ROMS')

if freq,
 axis([5.e-4 Inf 0 0.05])
 xlabel('Frequency (Hz)');
else
 axis([2.5 15 0.003 0.025])
 xlabel('Period (min)');
 x=[3:1:10 15];
 set(gca,'XTick',x)
 set(gca,'XTickLabel',sprintf('%3i|',x))
end
ylabel('Spectral density m^{3}.s^{-2}')
%title('Longshore current spectrum');

set(gcf,'PaperPositionMode','auto');
export_fig -transparent spec_f_comp.pdf


