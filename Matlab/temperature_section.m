%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Compute Fig. X from Treillou & Marchesiello (draft)
%   -> velocity average
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

% Parameters
FILE = '/scratch/users/treillou/IB09_IS_STRAT06';
%FILE = '/scratch/users/treillou/IB09_2024_3';
makepdf=0;

%% Compute
fn=[FILE '/rip_avg.nc'];

% Extract data
time = ncread(fn,'scrum_time');
x = findx_IB09(fn);
y = ncread(fn,'y_rho');y = y(1,:);
h = ncread(fn,'h'); hr = h(:,1); h = hr(1) * linspace(-1,0,10);
ixt = find(abs(time-5500) == min(abs(time-5500)));
ixe = find(abs(time-6100) == min(abs(time-6100)));

tpas01 = ncread(fn,'temp');
tpas01 = tpas01(:,:,:,ixt:ixe);%.*hr/10;
tpas01 = squeeze(mean(tpas01,4));
zeta = ncread(fn,'zeta'); zeta = squeeze(zeta(:,:,ixe));
ttt = squeeze(mean(tpas01, 2));
zr = squeeze(mean(zeta, 2));


ht = zeros(size(ttt));
xt = zeros(size(ttt));
for j = 1:size(ht, 1)
    ht(j,:) = linspace(-hr(j), (zr(j)), 10);
    xt(j,:) = x(j) * ones(10,1);
end

%%
    
figure('Position', [100, 100, 800, 400], 'PaperOrientation', 'landscape');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

plot(x, -hr, 'k', 'LineWidth', 3); hold on;
contourf(xt, ht, ttt, 20, 'LineStyle', 'none');
plot(x, zr, 'r', 'LineWidth', 2);
    
xlabel('$x$ ($m$)','Interpreter','latex','FontSize',15);
ylabel('$z$ ($m$)','Interpreter','latex','FontSize',15);
xlim([-330, 0]);
plot(-80*ones(10),linspace(-1.489,0.,10),'LineWidth',2,'Color','black','LineStyle',':')
set(gca,'linewidth',2);
caxis([16,20]);
title('Temperature for idealized case ','Interpreter','latex','FontSize',15);
%title('Temperature for IB09','Interpreter','latex','FontSize',15);
colorbar()
