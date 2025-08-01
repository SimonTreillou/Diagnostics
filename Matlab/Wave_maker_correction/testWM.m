%% CORRECTED
wd=0;
gamma=3.3;
wp=2.1;
wa=0.099;
wds=0.1;
wf=2*pi/wp;
x0=0;
y0=0.;
time0=0.;
wd =wd * pi/180;
wds=wds * pi/180;
Ndir=31;
Nfrq=50*Ndir;
depth=1.07;
g=9.81;

fmin=0.2*wf;
fmax=5.0*wf;
df=(fmax-fmin)/Nfrq;
cff2=0.;
for iw=1:Nfrq
    wf_bry(iw)=fmin+iw*df;
    wk_bry(iw)=w_to_k(wf_bry(iw),depth);
end

sumspec=0.0;
for iw=1:Nfrq
    sigma=0.5*( 0.09*(1.+sign(wf_bry(iw)-wf))+ ...
        0.07*(1.-sign(wf_bry(iw)-wf)));
    cff0=exp(-0.5*((wf_bry(iw)-wf)/(sigma*wf))^2);
    cff1=0.3119*(wf^4)/(wf_bry(iw)^5) ...
         *exp(-1.25*(wf/wf_bry(iw))^4)*gamma^cff0;
    cff2=16*cff1*df;
    wa_bry(iw)=cff2;
    sumspec=sumspec+cff2;
end

[~,displacetheta]=min(abs(wf_bry-wf));

cff2=mod(displacetheta(1),Ndir);
cff4=0.0;

for jw=1:Nfrq
    wd_bry(jw)=mod(jw-cff2,Ndir);
    if (wd_bry(jw) <= 0.0) 
        wd_bry(jw)=wd_bry(jw)+Ndir;
    end
    wd_bry(jw)=(-1.0)^real(jw)*(-pi*0.5 + ...
               pi*(floor(real(wd_bry(jw))/2.0 - ...
                0.5))/(real(Ndir)-1.0));
    wd_bry(jw)=wd_bry(jw)+wd;
    if (wd_bry(jw) >= 0.5*pi) 
        wd_bry(jw)=0.5*pi;
    end
    if (wd_bry(jw) <= -0.5*pi) 
        wd_bry(jw)=-0.5*pi;
    end
    cff3=exp(-((wd_bry(jw)-wd)/max(1.5*wds,1e-12))^2);
    wa_bry_d(jw)=cff3;
    cff4=cff4+wa_bry_d(jw)*wa_bry(jw);
end

cff1=sumspec/dot(wa_bry_d,wa_bry);
S1=wa_bry_d'*cff1 * wa_bry/sumspec;
disp("Corr "+num2str(sum(trapz(wd_bry,S))))

for jw=1:Nfrq
    wa_bry_d(jw)=sqrt(wa_bry_d(jw)*cff1);
    wa_bry(jw)=wa*sqrt(wa_bry(jw)/sumspec);
end

S=wa_bry_d'*wa_bry;
disp("Corr "+num2str(sum(trapz(wd_bry,S))))

%%
el=30;
cff3=2*pi/el;
for jw=1:Nfrq
    cff1=wd_bry(jw);
    cff2=wk_bry(jw);
    if (cff3 < cff2)
        cff4=cff2*sin(cff1)/cff3;
        cff5=999.0;
        iw=1;
        if (abs(cff4) >= 1.0)
            while ((abs(int16(cff4)-cff5) > 3e-3)) & (iw < 10000)
                cff1=cff1-sign(cff1)*0.00005;
                cff5=cff2*sin(cff1)/cff3;
                iw=iw+1;
            end
        else
            cff1=0.0;
        end
        wd_bry2(jw)=cff1;
    else
        wd_bry2(jw)=0.0;
    end
end


for iw=1:Nfrq
    wpha_bry(iw)=rand*2.*pi;
    wkx_bry(iw)=wk_bry(iw)*cos(wd_bry(iw));
    wky_bry(iw)=wk_bry(iw)*sin(wd_bry(iw));
end

plot(wd_bry,wa_bry_d);
hold on
plot(wd_bry2,wa_bry_d);
hold off


S=wa_bry_d'*wa_bry;

Sw=trapz(wf_bry,S');
Sd=trapz(wd_bry2,Sw);
plot(Sw);


%% DEFAULT
wd=0;
gamma=3.3;
wp=2.1;
wa=0.099;
wds=0.1;
wf=2*pi/wp;
x0=0;
y0=0.;
time0=0.;
wd =wd * pi/180;
wds=wds * pi/180;
Nfrq=50;
Ndir=31;
depth=1.07;
g=9.81;
       

fmin=0.2*wf;
fmax=5.*wf;
df=(fmax-fmin)/Nfrq;
cff2=0.;
for iw=1:Nfrq
    wf_bry(iw)=fmin+iw*df;
    wk_bry(iw)=w_to_k(wf_bry(iw),depth);
end

cff3=0;
for iw=1:Nfrq
    sigma=0.5*( 0.09*(1.+sign(wf_bry(iw)-wf))+ ...
        0.07*(1.-sign(wf_bry(iw)-wf)));
    cff0=exp(-0.5*((wf_bry(iw)-wf)/(sigma*wf))^2);
    cff1=0.3119*(wf^4)/(wf_bry(iw)^5) ...
         *exp(-1.25*(wf/wf_bry(iw))^4)*gamma^cff0;
    cff2=16*cff1*df;
    wa_bry(iw)=cff2;
    cff3=cff3+cff2;
end
for iw=1:Nfrq
    wa_bry(iw)=wa*sqrt(wa_bry(iw)/cff3);
end

dmin=wd-30*pi/180;
dmax=wd+30*pi/180;
dd=(dmax-dmin)/Ndir;
cff4=0.;
clear wa_bry_d
clear wd_bry
for jw=1:Ndir
    wd_bry(jw)=dmin+jw*dd;
    cff3=exp(-((wd_bry(jw)-wd)/(1.5*max(wds,1.-12)))^2);
    wa_bry_d(jw)=cff3;
    cff4=cff4+cff3;
end

disp(sum(wa_bry_d/cff4));
disp(sum(sqrt(wa_bry_d/cff4)))

for jw=1:Ndir
    wa_bry_d(jw)=sqrt(wa_bry_d(jw)/cff4);
end

S=wa_bry_d'*wa_bry;
disp("Def "+num2str(trapz(wf_bry,trapz(wd_bry,S))))
