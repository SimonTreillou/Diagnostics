Nf=[50,70,100,200,300,400,500,600,700,1000,1500,2000,3000,4000,5000,7000,10000];
i=1;
wa=0.23;
Tp=13;
gamma=20.;
for n=Nf
    alpha=0.3119;
    [w,S]=jonswap_spectrum(alpha,n,Tp,gamma,1);
    plot(w,S);
    disp("Int:"+trapz(w,S));
    hold on
    S=wa*S;
    ints(i)=trapz(w,S);
    i=i+1;
end
figure(); semilogx(Nf,ints); 
xlabel('Nfrq'); ylabel('$\int S(\omega)d\omega$','Interpreter','latex');
grid(); title('Spectrum energy as a function of freq. resolution')
disp("Attendu: "+string(wa*2*sqrt(2)));


%%

Nf=[50,70,100,200,300,400,500,600,700,1000,1500,2000,3000,4000,5000,7000,10000];
i=1;
wa=0.23;
wp=13;
gamma=20.;
for n=Nf
    clear wf_bry
    clear wa_bry
    wf=2*pi/wp;
    
    fmin=0.2*wf;
    fmax=5.0*wf;
    df=(fmax-fmin)/n;
    
    cff2=0.;
    for iw=1:n
        wf_bry(iw)=fmin+iw*df;
    end
    
    cff3=0;
    for iw=1:n
        sigma=0.5*( 0.09*(1.+sign(wf_bry(iw)-wf))+0.07*(1.-sign(wf_bry(iw)-wf)));
        cff0=exp(-0.5*((wf_bry(iw)-wf)/(sigma*wf))^2);
        cff1=0.3119*(wf^4)/(wf_bry(iw)^5)*exp(-1.25*(wf/wf_bry(iw))^4)*gamma^cff0;
        cff2=16.*cff1*df;  % integral must be 1
        wa_bry(iw)=cff2;
        cff3=cff3+cff2;
    end
    
    for iw=1:n
        wa_bry(iw)=wa*sqrt(wa_bry(iw)/cff3); % normalize
    end
    ints(i)=trapz(wf_bry,wa_bry);
    i=i+1;
end
figure(); semilogx(Nf,ints); 
xlabel('Nfrq'); ylabel('$\int S(\omega)d\omega$','Interpreter','latex');
grid(); title('Spectrum energy as a function of freq. resolution')
disp("Attendu: "+string(wa*2*sqrt(2)));

%%
clear wf_bry
clear wa_bry
Nfrq=100;
wp=13;
wa=0.23;
wf=2*pi/wp;

fmin=0.2*wf;
fmax=5.0*wf;
df=(fmax-fmin)/Nfrq;

cff2=0.;
for iw=1:Nfrq
    wf_bry(iw)=fmin+iw*df;
end

cff3=0;
for iw=1:Nfrq
    sigma=0.5*( 0.09*(1.+sign(wf_bry(iw)-wf))+0.07*(1.-sign(wf_bry(iw)-wf)));
    cff0=exp(-0.5*((wf_bry(iw)-wf)/(sigma*wf))^2);
    cff1=0.3119*(wf^4)/(wf_bry(iw)^5)*exp(-1.25*(wf/wf_bry(iw))^4)*gamma^cff0;
    cff2=16.*cff1*df;  % integral must be 1
    wa_bry(iw)=cff2;
    cff3=cff3+cff2;
end

for iw=1:Nfrq
    wa_bry(iw)=wa*sqrt(wa_bry(iw)/cff3); % normalize
end
plot(wf_bry,wa_bry);
disp(4*sqrt(trapz(wf_bry,wa_bry)));

