function [Z]=filter_my(lon,lat,z,R);

pio2=pi/2;
cff1=pio2/R;
R2=R^2;

%[x,y]=mercator(lon,lat,mean(mean(lon)),mean(mean(lat)));
x=lon; y=lat;

[L M]=size(x);
x1=reshape(x,L*M,1);
y1=reshape(y,L*M,1);
z1=reshape(z,L*M,1);

for i=1:L*M; 

%  disp([i,L*M])
  
  cosX=zeros(L*M,1);
  x0=x1-x1(i); y0=y1-y1(i);
  R0=x0.^2+y0.^2;
  D=find(R0<R2 & z1~=0);
  cosX(D)=cos(cff1*sqrt(R0(D)));
  cff=1/sum(cosX(D));
  Z1(i)=cff*sum(z1(D).*cosX(D));

end;
Z=reshape(Z1,L,M);

return


