%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
%fname1     = 'rip_his_3D_noCEW.nc';
fname1     = '../rip_his.nc';
fname2     = 'rip_his_2D_9h.nc'; 
fname3     = 'rip_his_mu006.nc';  % roms file name

makepdf   = 1;                    % make pdf file
smoo      = 1;
freq      = 0;

istr0     = 115; 
%
%======================================================================

figure('Units','pixels','Position',[500 500 600 500]);

%=================================================================

for ifile=1:1  %------------------------------- loop on ROMS files

if ifile==1,
 fname=fname1;
 col='k';
 istr=istr0;
 model3D   = 1;
 smoo=1;
 linew=2;
elseif ifile==2,
 fname=fname2;
 col='b';
 istr=istr0;
 model3D   = 0;
 smoo=1;
 linew=1;
else
 fname=fname3;
 col='m--';
 istr=istr0+96;
 model3D   = 0;
 smoo=0;
 linew=1;
end

nc=netcdf(fname);
tstr=10;
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1;  end;

iend=istr;
jstr=1; jend=256;

h=nc{'h'}(jstr:jend,istr:iend);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(jstr:jend,istr:iend)-xl+93;
y=nc{'y_rho'}(jstr:jend,istr:iend);
N=length(nc('s_rho'));

if model3D
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
x=x.*tukeywin(length(x),0.25);

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
if smoo
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
 loglog(f,Xm,col); hold on
else
 semilogy(p,Xm,col); hold on
end

end % ---------------------------- end loop ifile

%=================================================================

% VIDEO DATA

% --- LOAD DATA -----------------------------------
datafile='spectre_video_Mar11_2.fig';
openfig(datafile,'reuse','invisible')

h = gcf; %current figure handle
axesObjs = get(h, 'Children');        % axes handles
dataObjs = get(axesObjs, 'Children'); % handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');     % type of low-level graphics object 
periode = get(dataObjs, 'XData');     % data from low-level grahics objects
spectre = get(dataObjs, 'YData');

load spectre_video_Mar11.mat

load spec_video.mat
per_video=periode;
spe_video=spectre;

load spec_ADV.mat
per_ADV=periode;
spe_ADV=spectre;

load video_V_Mar11.mat      % <--

%spe_video=smooth(spe_video,5)./2.5;
%semilogy(per_video,spe_video);

%spe_ADV=smooth(spe_ADV,5)./2.5;
%semilogy(per_ADV,spe_ADV,'color',0.7*[1 1 1],'linewidth',2);


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

%end; end; % ------------------- j i


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
 loglog(f,Xm,col); hold on
else
 semilogy(p,Xm,col); hold on
end

%=================================================================

% FINALIZE PLOT

annotation('textarrow',[8 8]./20,[0.2 0.4], ...
           'String','Injection','fontsize',15,'color','r')

hold off
%grid on
legend('ROMS','ADV','VIDEO')

if freq,
 axis([5.e-4 Inf 0 0.05])
 xlabel('Frequency (Hz)');
else
 axis([0 20 0.002 0.02])
 xlabel('Period (min)');
end
ylabel('Spectral density m^{3}.s^{-2}')
%title('Longshore current spectrum');

set(gcf,'PaperPositionMode','auto');
export_fig -transparent spec_f_comp.pdf


