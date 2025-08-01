%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make plot from the results of the SHOREFACE test case
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
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
%       Application of the CROCO embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%  Patrick Marchesiello, IRD 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%close all
%================== User defined parameters ===========================
%
% --- model params ---
%
fname1 = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname2 = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

varname   = 'Akv';              % var name [ 'u' 'w' 'Akv' 'Temp' ]

makemovie = 0;                % make movie using QTWriter
makepdf   = 0;                % make pdf file
%
%======================================================================

xindex=160;
g = 9.81;

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

fname=fname2;
nc=netcdf(fname);
tindex=38;
length(nc{'scrum_time'}(:)); % reads last record
%tindex=5; 




x=squeeze(nc{'x_rho'}(1,:));
y=squeeze(nc{'y_rho'}(:,1));
time=squeeze(nc{'scrum_time'}(:));
tke=squeeze(nc{'tke'}(:,10,:,xindex));

% Plotting
hf = figure('position',[1000 500 800 400]);
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren', 'LineWidth', 3);
s=contourf(y,time,tke);
set(s,'EdgeColor',[1,1,1]);
ylabel('Time ($s$)','Interpreter','latex','FontSize',15);
xlabel('$y$ ($m$)','Interpreter','latex','FontSize',15);
colorbar();
title('TKE at $x$='+string(x(xindex)-300)+'m','Interpreter','latex' ...
    ,'FontSize',15);
ax=gca; ax.YAxis.Exponent = 3;

return


[w,S]=jonswap_spectrum(0.1,300,2.0,3.3,0); 
dw=w(2)-w(1);
plot(x,y)
hold on
S1=S*dw/sum(S*dw);
plot(w,S1*0.25/sqrt(8))


[~,ix1]=min(abs(x-0.03));
[~,ix2]=min(abs(x-1.2));
4*sqrt(trapz(x(ix1:ix2),y(ix1:ix2)))
4*sqrt(trapz(w,S/sum(S*dw))/16*0.25/sqrt(8))
