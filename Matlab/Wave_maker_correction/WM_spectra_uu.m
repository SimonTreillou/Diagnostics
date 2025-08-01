hname1="WM-Corrected-S10";
hname2="WM-Default-S10";

hname1 = "BAKER_G2B";
hname2 = "BAKER_G2B_default";

nbfiles=2;

wd=15;%0;
wds=10;%30;
gamma=20.0;%2.0;
wa=0.23;%0.17;
wp=13.0;%6.7;
load /Users/simon/Code/BAKER/Dspread_G2B
load /Users/simon/Code/BAKER/Fspec_G2B

    
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
    sta=find(((xpos==120)));
    xpos=nc{'Xgrid'}(:);
    xpos=xpos(sta);
    ypos=nc{'Ygrid'}(:);
    ypos=ypos(sta);
    tend=length(nc{'scrum_time'}(:));
    zeta1=nc{'zeta'}(100:tend,:);
    zeta1=zeta1(:,sta);
    u=squeeze(mean(nc{'u'}(100:tend,:,:),3));
    u=u(:,sta);
    v=squeeze(nc{'v'}(100:tend,:,10));
    v=v(:,sta);
    time=nc{'scrum_time'}(100:tend);
    depth=55;
    disp(time(end));
    dt=time(2)-time(1);
    disp(dt);

    S=0;
    for j=1:size(u,2)
        u1=squeeze(u(:,j));
        [S1,f]=mycspd(u1,u1,2^7,1/dt);
        S=S+S1;
    end
    if n==1
        u11=u1;
    elseif n==2
        u22=u1;
    end
    S=S/size(u1,2);

    loglog(f,smooth(S,2));
    hold on


end
    grid("minor")