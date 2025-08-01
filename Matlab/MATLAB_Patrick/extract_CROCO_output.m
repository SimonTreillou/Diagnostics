close all
clear all

ifile='flume_his_A00.nc';
ofile='croco_A00.txt';

% 0.0   time
% 0.00  1
% 0.78  2
% 1.49  3
% 1.63  4 -
% 2.51  5
% 4.17  6
% 5.17  7 -
% 6.17  8
% 7.1   9
% 7.8   10
% 8.44  11
% 9.4   12 -
% 10.5  13
% 11.5  14 -
% 12.5  15
% 13.9  16 -
% 14.1  17 -
% 14.5  18
% 15.5  19
% 16.5  20
% 17.5  21
% 18.5  22
% 19.5  23
% 20.0  24
 
fs=100;

%xo=[0.   0.78 1.49 1.63 2.51 4.17 5.17 6.17 7.10 7.80 ...
%    8.44 9.40 10.5 11.5 12.5 13.9 14.1 14.5 15.5 16.5 ...
%    17.5 18.5 19.5 20.0];
xo=[0. 1.63 5.17 9.4 11.5 13.9 14.1];
to=[40.:0.01:70]';

nc=netcdf(ifile);

zi=squeeze(nc{'zeta'}(:,2,:));
xi=squeeze(nc{'x_rho'}(2,:));
ti=nc{'scrum_time'}(:);
close(nc)
[LT LX]=size(zi);
for k=1:LT
 zo(k,:)=interp1(xi,zi(k,:),xo);
end
dataout=[to zo];
dataout(1,:)=[];

save(ofile,'dataout','-ascii')

plot(ti,zo(:,7));
grid on
axis([40 70 -0.1 0.2])
