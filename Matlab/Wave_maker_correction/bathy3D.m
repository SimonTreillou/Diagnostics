load Data/HR16_fig1_bathy

seuils = [-500,-230, -190, -158,-131.7,-110,-73,-29,20];
[Lx,Ly] = size(HR16fig1bathy);
h = zeros([Lx,1]);
bathy3d = [HR16fig1bathy h];

for i=1:(length(seuils)-1)
    bathy3d((bathy3d(:,1)<seuils(i+1))&(bathy3d(:,1)>seuils(i)),3) = 8 - i; 
end

x = bathy3d(:,1);
y = bathy3d(:,2);
z = bathy3d(:,3);
[X,Y,Z] = meshgrid(x,y,z);



[XX,YY] = meshgrid(x,y);
[XXq,YYq] = meshgrid(0:1:350,0:1:1600);
%[XX,YY] = meshgrid(linspace(min(X),max(X),100),linspace(min(Y),max(Y),72));
%first grid 
ZZ = griddata(x,y,z,XX,YY);

% new interpolated grid 
lonx=min_x:20:max_x;
laty=min_y:20:max_y;
[Xi,Yi]=meshgrid(lonx,laty);
bathymetry=interp2(XX,YY,ZZ);
%surf(XX,YY,ZZ);
%ZZ(ZZ>0)=nan; 
figure(1)
%pcolor(XX,YY,ZZ),shading interp ,colorbar;
pcolor(Xi,Yi,bathymetry),shading flat ,colorbar;  % if error transpose bathymetry and see
saveas(gcf,'bat_fine','eps');



[y0,x0,z0] = ndgrid(d,d,d);
XI = [x0(:) y0(:) z0(:)];
YI = griddatan(X,Y,XI);

