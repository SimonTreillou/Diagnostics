fname='rip_his.nc'
dir='../'
nc=netcdf([dir,fname])

h=nc{'h'}(:);
x=nc{'x_rho'}(:,:);
xi=nc{'xi_rho'}(:,:);
x=x(1,:);
xi2= 650-xi;

%%
clear y; 


for i=1:length(xi)
    xo = xi2(i);
    if xo>=(410)
        y(i)=7;
    elseif xo<100
        y(i)=-(xo-550)*0.02;
    elseif (xo>=100)&&(xo<410)
        y(i)=-(sinf(C1,(xo-550)) + sinf(C2,(xo-550)) + sinf(C3,(xo-550)));
    end
end


plot(xi2,y);
%%
xit=linspace(240,550,300)
xito=xit-550
t = sinf(C1,(xito)) + sinf(C2,(xito)) + sinf(C3,(xito))
plot(xit,t)

%%
xiR = linspace(1,240,240);
xiR2 = 240 - xiR + 85;

%%
clear y; 
xs = 85.;
db=50;
alpha=0.025;
lambda=10;

for i=1:length(xiR2)
    xx = xiR2(i);
    y(i) = -4.6 - 1.5*exp(-6*(((xx-xs-db)/db)^2)) +3.1*(1+tanh(0.025*(xx-xs))) +0.014*(xx+log(cosh(alpha*(xx-xs))/cosh(alpha*xs))/alpha);
    y(i) = min(y(i), 11.3-0.88e-5*(abs(xx-700))^2.17);
end

plot(xiR2,y)
