wtime=[34*60+56,19*60,11*60+49,5*60+51,3*60+26];
cpus=[1,2,4,8,16];
th=[1,2,4,8,16];
obs=wtime(1)./wtime;

loglog(cpus,th)
hold on
loglog(cpus,obs)