%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Compute Fig. X from Treillou & Marchesiello (draft)
%   -> vorticity movie
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%================== User defined parameters ===========================
%
% --- model params ---
%
% Simulations
dirpath = '/scratch/users/treillou/';
fnames{1}='IB09_2024_2';
fnames{2}='IB09_2024_3';
fnames{3}='IB09_2024_4';
fnames{4}='IB09_2024_5';
fnames{5}='IB09_2024_5_repeat';
fnames{6}='IB09_2024_6';
nbfiles=length(fnames);
makepdf   = 0;             % make pdf file
Dcrit=0.05;
%
%======================================================================

% ---------------------------------------------------------------------
% --- get grid from numerical model ---
% ---------------------------------------------------------------------

%%
% mean x= two surfzone widths

for i=1:nbfiles
    fname=[dirpath,fnames{i},'/rip_avg.nc'];
    if i==1
        time=ncread(fname, 'scrum_time');
        x=findx_IB09(fname);
        [~,ixIS]=min(abs(x+160)); %limit of the SZ-IS
        h=ncread(fname,'h'); h=h(1:ixIS,1)/10;
        lvls=linspace(1,10,10);
        hIS=h(ixIS);
        
        t=ncread(fname,'tpas01');
        t=squeeze(t(ixIS,:,:,:));
        
        for l=lvls
            if (l+1)*hIS>1
                break
            end
        end
        L=l;
        
        t_1m=t(:,L,:) + (1-hIS*L)*(t(:,L+1,:)-t(:,L,:))/(hIS*(L+1)-hIS*L);
        res{i}=zeros([10,1]);
        for l=lvls
            t_l=squeeze(t(:,l,:));
            res{i}(l)=nanmean(t_l(t_l>2));
        end


    else
        time=cat(1,time,ncread(fname,'scrum_time'));
        t=ncread(fname,'tpas01');
        t=squeeze(t(ixIS,:,:,:));
        
        for l=lvls
            if (l+1)*hIS>1
                break
            end
        end
        L=l;
        
        t_1m=t(:,L,:) + (1-hIS*L)*(t(:,L+1,:)-t(:,L,:))/(hIS*(L+1)-hIS*L);
        res{i}=zeros([10,1]);
        for l=lvls
            t_l=squeeze(t(:,l,:));
            res{i}(l)=nanmean(t_l(t_l>2));
        end
    end
end

%%

C_obs=[2.29,2.968,3.935,4.903,5.645];
C_obsp=[3.839,4.903,6.194,7.29,8.097];
C_obsm=[0.742,1.161,1.806,2.452,3.065];
h_obs=[-3.017,-2.51,-2.012,-1.506,-1.008];
   

%%
RES=0;
for i=1:nbfiles
    RES=res{i}/nbfiles+RES;
end
plot(RES,lvls*hIS-hIS*10,'LineWidth',3,'Color','r')
hold on
plot(C_obs,h_obs,'d-','LineWidth',3,'Color',[0.4,0.4,0.4])
plot(C_obsm,h_obs,'--','LineWidth',3,'Color',[0.6,0.6,0.6])
plot(C_obsp,h_obs,'--','LineWidth',3,'Color',[0.6,0.6,0.6])
xlabel('$\bar{D}$ $(ppb)$','Interpreter','latex','FontSize',15)
ylabel('$z$ $(m)$','Interpreter','latex','FontSize',15)
xlim([0,10])
grid()
ylim([-3.25,-0.75])
legend('CROCO','Obs.','Interpreter','latex','FontSize',15)
set(gca,'linewidth',2)
set(gca,'layer','top')
