close all
clear all
%
%======================================================================
%  COMPUTE LINEAR SUMMATION OF WAVE SPECTRUM
%======================================================================
%
% Vyzikas et al., 2018: The evolution of free and bound waves during 
% dispersive focusing in a numerical and physical flume. Coastal 
% Engineering, 132, 95–109.
%
scale=0.154/0.05;  % scaling for flume experiment
h=1;               % Tank depth
g=9.8;

fname='flume_his_A00.nc';
%======================================================================
%
% Read CROCO waves     ------------------
%
nc=netcdf(fname);
z=100*squeeze(nc{'zeta'}(1:10:end,2,:));
x=squeeze(nc{'x_rho'}(2,:));
t=nc{'scrum_time'}(1:10:end);
close(nc)
[LT LX]=size(z);
maxz=max(max(z));
[it ix]=find(z==maxz);
maxt=t(it);
maxx=x(ix);
maxz=z(it,ix);
disp(['max level ',num2str(maxz),' at time ',num2str(maxt), ...
                           ' and position ',num2str(maxx)]);
%
%
% Read linear waves     -------------
%
load datwaves.txt
amp=datwaves(:,1);
frq=datwaves(:,2);
pha=datwaves(:,3);
kw =datwaves(:,4);
%
t0=64;        % focal time
x0=14.1;      % and position
zl=zeros(length(t),length(x));
tl=t-(maxt-t0);
xl=x-(maxx-x0);
for i=1:length(t);
  zl(i,:)=100*scale*sum(amp.*cos(kw*(xl-x0) -frq*(tl(i)-t0) -pha));
end
%
%  MAKE PLOT
%
movObj = QTWriter('rogue.mov');
hf = figure('position',[500 500 700 150]);
axis tight; set(hf,'DoubleBuffer','on');

for i=50:it;

 hold on
 plot(x-2,z(i,:),'color','k','linewidth',2);
 plot(x-2,zl(i,:),'color',[0.7 0.7 0.7],'linewidth',2);
 grid on
 axis([0 28 -20 +20])
 %xlabel('Distance [m]')
 %ylabel('Water level [cm]')
 %title('Rogue wave')
 set(gca,'fontsize',14)
 hold off

 movObj.FrameRate =  10;
 writeMovie(movObj,getframe(hf));
 clf('reset')

end

movObj.Loop = 'loop'; % Set looping flag
close(movObj);        % Finish writing movie and close file


