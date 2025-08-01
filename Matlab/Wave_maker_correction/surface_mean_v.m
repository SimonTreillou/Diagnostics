clear all
%close all
%================== User defined parameters ====================
%
% --- model params ---
%

fname = '/Users/simon/Code/CONFIGS/IB09_21AVR_S10/rip_avgtot.nc'; 
fname = '/Users/simon/Code/CONFIGS/WM-Default-S10/rip_avg.nc';
fname = '/Users/simon/Code/CONFIGS/WM-Corrected-S10/rip_avg.nc';

makepdf   = 0;                       % make pdf file
%
%===============================================================
model3D=1;

nc=netcdf(fname);
tstr=1 ;
tend=length(nc{'scrum_time'}(:));
if tend<tstr; tstr=1;  end;
%tend=88;

h=nc{'h'}(:);
xl=nc{'xl'}(:);
x=nc{'x_rho'}(:)-xl+50;
y=nc{'y_rho'}(:);
N=length(nc('s_rho'));



v=nc{'vbar'}(tstr:tend,:,:);
for it=1:tend-tstr+1
  vr(it,:,:)=v2rho_2d(squeeze(v(it,:,:)));
end

mv=squeeze(mean(vr));

%% Plotting
figure('Position',  [100, 1000, 400, 800]);
Jstr=1;
pcolor(x(Jstr:end,:),y(Jstr:end,:)-y(Jstr),mv(Jstr:end,:))
shading flat;
caxis;
colorbar;
l=max(-min(min(mv(Jstr:end,:))),max(max(mv(Jstr:end,:))));
l=0.25;
cmin=-l;
cmax=l;
cmap=redwhiteblue(cmin,cmax);
colormap(cmap);
caxis([cmin,cmax])
xlim([-300 0]);
%ylim([0 800]);
%set(gcf, 'Position',  [100, 1000, 400, 800]);
%title('Mean longshore velocity')
xlabel('$x$ (m)','Interpreter','latex','FontSize',15)
ylabel('$y$ (m)','Interpreter','latex','FontSize',15)
print(gcf, '-dpdf', './Figures/Mean-V-Corrected.pdf');



