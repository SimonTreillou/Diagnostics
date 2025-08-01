%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ====================
%
% --- model params ---
%
fname1     = '/Users/simon/Code/supp/rip_avg.nc';    % croco file name
%fname1     = '/Users/simon/Code/CONFIGS/testIB09GPP/rip_avg.nc';    % croco file name
fname2     = '../rip_avg_2D_SC.nc';
fname3     = '../rip_avg_2D_LC.nc';

nbfiles   = 1;

makepdf   = 0;                       % make pdf file
smoo      = 1;
freq      = 1;

istr0     = 160; 
%
%===============================================================

figure('Units','pixels','Position',[500 500 600 500]);
set(gca,'FontSize',15)

for ifile=1:nbfiles  %--------------------- loop on ROMS files

if ifile==1,
 fname=fname1;
 col='k';
 istr=istr0;
 %model3D=1;
 %smoo=1;
 %nsmoo=8;
 linew=2;
elseif ifile==2,
 fname=fname2;
 col='b';
 istr=istr0;
 %model3D=1;
 %smoo=1;
 %nsmoo=8;
 linew=2;
else
 fname=fname3;
 col='m--';
 istr=istr0; %+96;
 %model3D=0;
 %smoo=1;
 %nsmoo=6;
 linew=1;
end

smoo=1;
nsmoo=10;
model3D=1;

nc=netcdf(fname);
tstr=10;
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1;  end;
tend=21;

xl=nc{'xl'}(:);
x=nc{'x_rho'}(1,:)-xl;
%istr=160; %min(find(x>-100.));
istr=1;
h=nc{'h'}(:);
[M L]=size(h);
iend=istr+40;
jstr=2; jend=M-1;
clear x h;

h=nc{'h'}(jstr:jend,istr:iend);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(jstr:jend,istr:iend)-xl; %+93;
y=nc{'y_rho'}(jstr:jend,istr:iend);
N=length(nc('s_rho'));

if model3D,
 u=squeeze(nc{'u'}(tstr:tend,N,jstr:jend,istr:iend));
 v=squeeze(nc{'v'}(tstr:tend,N,jstr:jend,istr:iend));
else
 u=nc{'ubar'}(tstr:tend,jstr:jend,istr:iend);
 v=nc{'vbar'}(tstr:tend,jstr:jend,istr:iend);
end

v=v;

[N M L]=size(v);
vm=squeeze(v(:,jend-jstr+1,iend-istr+1));
%dtm=440*0.025;   % his 
dtm=880*0.025; % avg
tm=[0:dtm:length(vm)*dtm-dtm]';

vv=vm; t=tm; dt=dtm;

n = 2^(ceil(log2(N))); 
Xm=zeros(n,1);


for j=1:M; for i=1:L; % ---------- j i

%
% Get detrended and windowed time series
%
vv=squeeze(v(:,j,i));
vv=vv-mean(vv);
x=detrend(vv,'linear');
%x=x.*hann(length(x),'periodic'); 
%x=x.*tukeywin(length(x),0.25);

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
  Xm=smooth(Xm,nsmoo);
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
 loglog(f,Xm,col,'linewidth',linew); hold on
else
 %semilogy(p,Xm,col,'linewidth',linew); hold on
 loglog(p,Xm,col,'linewidth',linew); hold on
end

end % ---------------------------- end loop ifile

hold off
grid on
legend('CROCO 3D SC','CROCO 2D SC','CROCO 2D LC', ...
       'location','northwest');

if freq,
 axis([5.e-4 Inf 0 Inf])
 xlabel('Frequency (Hz)');
else
 %axis([0 20 0.002 0.04])
 xlabel('Period (min)');
 %x=[3:1:10 15];
 %set(gca,'XTick',x)
 %set(gca,'XTickLabel',sprintf('%3i|',x))
end
ylabel('Spectral density m^{2}.s^{-2}.Hz^{-1}')
%title('Longshore current spectrum');

if makepdf,
 set(gcf,'PaperPositionMode','auto');
 export_fig -transparent spec_f_ke.pdf
end

