clear all
close all
clf
xl =650;
x=0:5:xl;

C1=[9.509, 0.002914, 0.014];
C2=[0.5618,0.02377,-2.684];
C3=[0.1274,0.05868,1.32];
% 
% D1=[7.469, 0.004057, -0.204];
% D2=[1.59, 0.004936, -4.481];
% D3=[0.5793, 0.01982, 3.158];
% D4=[0.09405, 0.03475, -0.8076];
% D5=[0.0926, 0.04956, 0.1163];


E1=[7.528,0.003551,-0.01278];
E2=[0.6487,0.01703,-3.534];
E3=[0.1431,0.02953,4.913];


load Data/HR16bathy.mat
% tofitx = HR16bathy(:,1);
% tofity = HR16bathy(:,2);
% [~,indexmax] = min(abs(tofity));
% [~,indexmin] = min(abs(tofitx+300));
% tofitx = tofitx(indexmin:indexmax);
% tofity = tofity(indexmin:indexmax);
% tofity(1:3) = -7;
% 
% 
% D1=[7.341, 0.008515, 0.3361];
% D2=[2.594, 0.01699, -1.947];
% D3=[0.6649, 0.02734, -3.056];
% D4=[0.1026, 0.06184, -4.476];
% 
for i=1:length(x)
    if x(i)>320
        y(i)=7-0.05;
        %y2(i)=y(i)+0.05;
    elseif x(i)<50
        y(i)=-(50-x(i))*0.02;
        %y2(i)=y(i);
    else
        y(i)=-(sinf(C1,50-x(i)) + sinf(C2,50-x(i)) + sinf(C3,50-x(i)));
        %y2(i)=-(sinf(D1,50-x(i)) + sinf(D2,50-x(i)) + sinf(D3,50-x(i)) ...
        %    + sinf(D4,50-x(i)));
    end
end


for i=1:length(x)
    if x(i)>410
        y2(i)=7;
    elseif x(i)<50
        y2(i)=-(50-x(i))*0.02;
    else
        y2(i)=-(sinf(E1,50-x(i)) + sinf(E2,50-x(i)) + sinf(E3,50-x(i)));
    end
end


plot(x-xl+50,-flip(y),'LineWidth',3);
hold on
%plot(x-xl+50,-flip(y2),'LineWidth',3);
scatter(HR16bathy(:,1),HR16bathy(:,2),'MarkerFaceColor','blue','MarkerFaceAlpha',0.5);
%scatter(tofitx,tofity);
xlim([-300,50]);
plot(x-xl+50,(x-xl+50)*0.02)
plot(x-xl+50,(x-xl+50)*0.027)
legend('BathySimon','Obs','Location','best');



%%
% fname='/Users/simon/Code/IB09/IB09PATPSdx1bathynewang7strat/rip_his.nc';
% nc=netcdf(fname);
% 
% hr=squeeze(nc{'h'}(1,:));
% L=length(hr);
% xr=squeeze(nc{'x_rho'}(1,:));
% 
% plot(xr-300,-hr);
% hold on
% scatter(HR16bathy(:,1),HR16bathy(:,2))
% legend('refined bathymetry','obs.')