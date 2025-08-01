
hname1="BAKER_G2B_default";
hname2="BAKER_G2B";
nbfiles=2;

wd=15;%0;
wds=10;%30;
gamma=20.0;%2.0;
wa=0.23;%0.17;
wp=13.0;%6.7;
load /Users/simon/Code/BAKER/Dspread_G2B
load /Users/simon/Code/BAKER/Fspec_G2B

SMout1=0;SMout2=0;
for n=1:nbfiles
    if n==1
        hname=hname1;
    elseif n==2
        hname=hname2;
    end
    dirpath = '/Users/simon/Code/CONFIGS/'+hname+'/';
    stname  = strcat(dirpath,'stations.nc');
    nc=netcdf(stname,'r');

    xpos=nc{'Xgrid'}(:);
    ypos=nc{'Ygrid'}(:);
    staT=find(((xpos==40)));
    for ss=1:length(staT)
        sta=staT(ss);
        xpos=nc{'Xgrid'}(:);
        xpos=xpos(sta);
        ypos=nc{'Ygrid'}(:);
        ypos=ypos(sta);
        tend=length(nc{'scrum_time'}(:));
        zeta1=nc{'zeta'}(100:tend,:);
        zeta1=zeta1(:,sta);
        u1=squeeze(nc{'u'}(100:tend,:,10));
        u1=u1(:,sta);
        v1=squeeze(nc{'v'}(100:tend,:,10));
        v1=v1(:,sta);
        time=nc{'scrum_time'}(100:tend);
        depth=1.07;
        disp(time(end));

        
        %
        % DATA
        %
        ID.data=[zeta1 u1 v1];
        xtmp=repmat((xpos),size(ID.data,2)/size(xpos,1),1)/10;
        ytmp=repmat((ypos),size(ID.data,2)/size(ypos,1),1)/10;
        ztmp=zeros(size(ID.data,2),1);
        ID.layout=[xtmp'; ytmp'; ztmp'];
        ztmp={'elev'};utmp={'velx'};vtmp={'vely'};
        for i=1:size(zeta1,2)-1
            ztmp{end+1}='elev';utmp{end+1}='velx';vtmp{end+1}='vely';
        end
        ID.datatypes=[ztmp utmp vtmp];
        ID.depth=depth;
        ID.depth=1.07;
        ID.fs=1/(time(2)-time(1));
        
        check_data(ID,1)
        
        %
        % SPECTRAL MATRIX STRUCTURE
        %
        nf=1000;
        nd=1000;
        SM.freqs=linspace(0.001,1.,nf);
        SM.dirs=linspace(-180,180,nd);
        SM.S=zeros(nf,nd);
        SM.axisdir=0.0;
        SM.funit='Hz';
        SM.dunit='cart';
        
        check_data(SM,2)
        
        %
        % ESTIMATION PARAMETER STRUCTURE
        %
        EP.method='EMEP';
        EP.dres=180;
        EP.iter=100;
        EP.smooth='ON';
        EP.nfft=64;
        %EP.nfft=128;
    
        check_data(EP,3)
        
        %
        % COMPUTING SPECTRUM
        %
        [SMout,EPout]=dirspec(ID,SM,EP,{'MESSAGE',0,'PLOTTYPE',0});
        disp(i)
        if n==1
            SMout1=SMout.S+SMout1;
            EPout1=EPout;
        elseif n==2
            SMout2=SMout.S+SMout2;
            EPout2=EPout;
        end
    end
end

%% PLOTTING FREQUENCY SPECTRUM

plot(Fspec_G2B(:,1),Fspec_G2B(:,2),'LineStyle',':','Color','k','LineWidth',2); %see Suanda 2016;
hold on
%title('Frequency spectrum at X='+string(xpos)+' and Y='+string(ypos))
ylabel('Frequency spectrum $D(f)$','Interpreter','latex','FontSize',15);
xlabel('$f$ (Hz)','Interpreter','latex','FontSize',15);

for i=1:nbfiles
    if i==1
        S=SMout1;
        c=[111,160,135,0.7*255]/255;
    elseif i==2
        S=SMout2/9;
        c=[192,24,135,0.7*255]/255;
    end
    FSpecObs=trapz(SMout.dirs,S');
    plot(SMout.freqs,smooth(FSpecObs,1),'LineWidth',3,'Color',c);
end
legend('Data','Default','Corrected','Interpreter','latex','FontSize',15);
grid("minor");
print(gcf, '-dpdf', './Figures/Baker-freq-comparison.pdf');

%% PLOTTING DIRECTIONAL SPECTRUM

plot(Dspread_G2B(:,1),Dspread_G2B(:,2)/max(Dspread_G2B(:,2)),'LineStyle',':','Color','k','LineWidth',2);
hold on
legend('Input','Default','Corrected','Interpreter','latex','FontSize',15);
ylabel('Directional spread $D(\theta)$','Interpreter','latex','FontSize',15);
xlabel('$\theta$ (deg)','Interpreter','latex','FontSize',15);

for i=1:nbfiles
    if i==1
        S=SMout1;
        c=[111,160,135,0.7*255]/255;
    elseif i==2
        S=SMout2;
        c=[192,24,135,0.7*255]/255;
    end
    [~,ix1]=min(abs(SMout.freqs<0.3));
    [~,ix2]=min(abs(SMout.freqs<0.8));
    dw=SMout.freqs(2)-SMout.freqs(1);
    DSpecObs=sum(S(ix1:ix2,:))*dw/9;

    plot(SMout.dirs,smooth(DSpecObs,1)/max(DSpecObs),'LineWidth',3,'Color',c);
end
legend('Data','Default','Corrected','Interpreter','latex','FontSize',15);
grid("minor");
print(gcf, '-dpdf', './Figures/Baker-dir-comparison.pdf');

