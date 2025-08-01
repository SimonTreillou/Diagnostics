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
warning off
%================== User defined parameters ===========================
%
% --- model params ---
%
fname     = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx1_tr1000_w011_sansinit_suite3/rip_avg.nc';  % croco file name
%fname     = '/Users/simon/Code/CONFIGS/CALMIP/IB09_dx1_tr1000_w011_sansinit_DB10.2_suite/rip_avg.nc';  % croco file name
%fname     = '../rip_avg_2D_SC.nc';
%fname     = '../rip_avg_2D_LC.nc';

model3D   = 1;
title0    = 'CROCO';

tpas      = 0;
perturb   = -100;             % plot longshore perturbation fields
makemovie = 0;             % make movie using QTWriter
makepdf   = 0;             % make pdf file

Dcrit     = 0.2;

%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

nc=netcdf(fname);
tindex=length(nc{'scrum_time'}(:)); % reads last record
t0=nc{'scrum_time'}(1)/60;

if makemovie,
 tstr= 1;
 tend=tindex;
 movObj = QTWriter('rip_vort.mov');
else
 tstr=input('Give record number:')+1;
 tstr=min(tindex,tstr);
 tend=tstr;
end

hf = figure('position',[500 500 400 700]);
set(gca,'FontSize',15)
axis tight; set(hf,'DoubleBuffer','on');
set(gca,'nextplot','replacechildren');
box on

for tindex=tstr:tend % ---------------------------------------------

 if tindex==tstr
  h=nc{'h'}(:);
  xl=nc{'xl'}(:);
  x=nc{'x_rho'}(:)-xl+50; %+93;
  y=nc{'y_rho'}(:);
  pm=nc{'pm'}(:);
  pn=nc{'pn'}(:);
  N=length(nc('s_rho'));
 end

 time=nc{'scrum_time'}(tindex)/60;
 zeta=nc{'zeta'}(tindex,:,:);
 if model3D,
  u=nc{'u'}(tindex,N,:,:);
  v=nc{'v'}(tindex,N,:,:);
  if tpas,
    t=nc{'tpas01'}(tindex,N,:,:);
  end
 else
  u=nc{'ubar'}(tindex,:,:);
  v=nc{'vbar'}(tindex,:,:);
 end
 ur=u2rho_2d(u);
 vr=v2rho_2d(v);

% Mask
 mask=ones(size(zeta));
 mask((h+zeta)<=Dcrit)=0;

% Vorticity
 vort=mask.*psi2rho(vorticity(u,v,pm,pn));

% Make lower resolution grid for vectors
 x1=min(min(x)); x2=max(max(x));
 y1=min(min(y)); y2=max(max(y));
 incvec=10; % resolution in m
 [x2,y2]=meshgrid([x1:incvec:x2],[y1:incvec:y2]);
 ur2=griddata(x,y,ur,x2,y2);
 vr2=griddata(x,y,vr,x2,y2);

 %============================================================
 % --- plot ---
 %=============================================================
 if tpas
  cmin=0; cmax=0.8; nbcol=20;
  AdvancedColormap('woyorm',20)
 else
  cmin=-0.1; cmax=-cmin; nbcol=20;
  AdvancedColormap('vbgswwyorm',20)
 end
 cint=(cmax-cmin)/nbcol;
%
 %set(gcf,'color','w');

 if tpas
  [C,hh]=contourf(x,y,t,[cmin:cint:cmax]); hold on
 else
  if perturb
   vort=vort-repmat(mean(vort),size(vort,1),1);
  end
  %[C,hh]=contourf(x,y,vort,[cmin:cint:cmax]); 
  pcolor(x,y,vort)
 end
 shading flat; 
 colorbar;

 I = ~(sqrt(ur2.^2+vr2.^2)<0.1);
 if perturb
  ur2=ur2-repmat(mean(ur2),size(ur2,1),1);
  vr2=vr2-repmat(mean(vr2),size(vr2,1),1);
 end
 %quiver(x2(I),y2(I),ur2(I),vr2(I)); hold off

 axis([-200 10 0 max(y(:,1))])
 caxis([cmin cmax])
 thour=floor(time/60);
 tmin=floor(time)-60*thour;
 clock=[num2str(thour),' h ',num2str(tmin),' min '];
 %title(['Time: ',clock])
 %title(title0)
 disp(['Plot for time: ',clock,' (Record: ',num2str(tindex),')'])
 hold off

 if makemovie,  
  % Write each frame to the file
  movObj.FrameRate =  5;
  writeMovie(movObj,getframe(hf));
  clf('reset')
 end

%----------------------------------

end % time loop

close(nc);

if makemovie,  
    movObj.Loop = 'loop'; % Set looping flag
    close(movObj);        % Finish writing movie and close file
end

if makepdf, 
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'color','w');
    if tpas,
      export_fig -transparent model_tpas.pdf
    else
      export_fig -transparent model_vort.pdf
    end
end

return





