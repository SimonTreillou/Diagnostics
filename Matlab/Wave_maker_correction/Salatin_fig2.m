alpha=0.005;
Nf=1550;
Ndir=1550;
T=6.7;
wd=0*(pi/180);

gamma=2.0;

wds=30*(pi/180);
Ndir=320;

[w,S]=jonswap_spectrum(alpha,Nf,T,gamma);
Hs=sqrt(trapz(w,S)*16)
plot(w,S);
hold on

%% 

clear wa_bry_d;
clear wd_bry;
dmin=wd-90*(pi/180);
dmax=wd+90*(pi/180);
dd=(dmax-dmin)/Ndir;

cff4=0;
for i=1:Ndir
    wdbryi=dmin+i*dd;
    wd_bry(i)=wdbryi;
    cff3=exp(-((wdbryi-wd)^2/(max(1.5*wds,1.e-12)^2))); %*1/(sqrt(2*pi)*wds);
    wa_bry_d(i)=cff3;
    %wa_bry_d2(i)=exp(-((wdbryi-wd)/max(wds,1.e-12))^2);
    cff4=cff4+cff3;
end

for i=1:Ndir
    wa_bry_d(i)=sqrt(wa_bry_d(i)/cff4);
end
plot(wd_bry*180/(pi),wa_bry_d)%*180/(pi))
plot(wd_bry,wa_bry_d)%*180/(pi))

CS=cumtrapz(wa_bry_d)/trapz(wa_bry_d);
disp(trapz(wd_bry,wa_bry_d));
hold on

%% Plot 3D spectra 
figure();
tes=S.*wa_bry_d';
[X,Y]=meshgrid(w,wd_bry);
h=surf(X,Y,tes);
set(h,'EdgeColor','none');
caxis([0, 0.8]);
colorbar();
xlim([0.05,0.25]);
ylim([-1.5,1.5]);