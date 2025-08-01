%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SHOREFACE test case
% 
%  Further Information:  
%  http://www.romsagrif.org/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Ref: Penven, P., L. Debreu, P. Marchesiello and J.C. McWilliams,
%       Application of the ROMS embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%  Patrick Marchesiello, IRD 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname      = '../../CONFIGS/testIB09GPP2/rip_his.nc';  % roms file name
%fname1     = 'rip_his_Cd005_wind.nc';
fname1     = 'rip_his_mu003_wind.nc';
fname2     = 'rip_his_mu001.nc';
fname3     = 'rip_his_mu006.nc';
fname5     = fname;

makepdf   = 0;             % make pdf file
Dcrit=0.2;

%
%======================================================================

nc=netcdf(fname);

tstr=30;
tend=length(nc{'scrum_time'}(:));

h=nc{'h'}(:);
x=nc{'x_rho'}(:)-675;
y=nc{'y_rho'}(:);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
N=length(nc('s_rho'));

close(nc)

% ------ DATA -------
% xZ/Z: topo
% xB/B: breaking
% xV/V: longshore current
%load('DATA/GPP2014_topo_currents_breaking.mat')
x0=75;
xZ=-xZ+x0; Z=Z-0.8;
xV=-xV+x0;
xB=-xB+x0;

clear xV V
%load('DATA/GPP2014_currents.mat')
V=mean(MCRad(:,4:6),2);
xV=-X+50;

clear xB B
%load('DATA/BreakingGPP2014_20141219_2.mat')
B=b/max(b);
xB=xb+50;

figure('Units','pixels','Position',[500 500 700 400]);

plot(xV,V,'r--','linewidth',2); hold on


for i=1:3; %----------------------------------------------

if i==1,
 fname=fname1;
 col='b';
 lw=2;
elseif i==2,
 fname=fname2;
 col='b--';
 lw=1;
elseif i==3,
 fname=fname3;
 col='b+';
 lw=1;
elseif i==4,
 fname=fname4;
 col='b*';
 lw=1;
elseif i==5,
 fname=fname5;
 col='b:';
 lw=1;
end

nc=netcdf(fname);
v=squeeze(mean(nc{'vbar'}(tstr:tend,:,:)));
mv=squeeze(mean(v));
mx=squeeze(mean(x));
close(nc)

plot(mx,-mv,col,'linewidth',lw); hold on


end % ------------------------------------------- fname loop

hold off

grid('on')
axis([-200 -20 0 1.3])
ylabel('Longshore Current [m/s]');
legend('Video',...
       '\mu 0.003', ...
       '\mu 0.001', ...
       '\mu 0.006', ...
       'Location','northwest')

export_fig -transparent V_model_test.pdf

return



