import numpy as np

def csf(sc, theta_s,theta_b):
    if (theta_s > 0 ):
        csrf=(1-np.cosh(sc*theta_s))/(np.cosh(theta_s)-1);
    else:
        csrf=-sc**2;
    if (theta_b > 0):
        h = (np.exp(theta_b*csrf)-1)/(1-np.exp(-theta_b));
    else:
        h  = csrf;
    return h

def zlevs(h,zeta,theta_s,theta_b,hc,N,ttype,vtransform):
    [M,L]=np.shape(h)
    sc_r=np.zeros([N,1]);
    Cs_r=np.zeros([N,1]);
    sc_w=np.zeros([N+1,1]);
    Cs_w=np.zeros([N+1,1]);

    if (vtransform == 2):
        ds=1./N
        if ttype=='w':
            sc_w[0] = -1.0
            sc_w[-1] =  0
            Cs_w[0] = -1.0
            Cs_w[-1] =  0
            
            sc_w[1:-1] = ds*(np.arange(1,N)-N)
            Cs_w=csf(sc_w, theta_s,theta_b);
            N=N+1;
        else:
            sc= ds*(np.arange(1,N+1)-N-0.5)   
            Cs_r=csf(sc, theta_s,theta_b);
            sc_r=sc;

    else:
        cff1=1./np.sinh(theta_s);
        cff2=0.5/np.tanh(0.5*theta_s);
        if ttype=='w':
            sc=(np.arange(0,N+1)-N)/N;
            N=N+1;
        else:
            sc=(np.arange(1,N+1)-N-0.5)/N;
        Cs=(1.-theta_b)*cff1*np.sinh(theta_s*sc)+theta_b*(cff2*np.tanh(theta_s*(sc+0.5))-0.5)


    Dcrit=0.01; 
    zeta[zeta<(Dcrit-h)]=Dcrit-h[zeta<(Dcrit-h)];
    

    z=np.zeros([N,M,L])
    
    if (vtransform == 2):
        if ttype=='w':
            cff1=Cs_w;
            cff2=sc_w+1;
            sc=sc_w;
        else:
            cff1=Cs_r;
            cff2=sc_r+1;
            sc=sc_r;
        h2=(np.abs(h)+hc);
        cff=hc*sc;
        h2inv=1/h2;
        for k in range(N):
            z0=cff[k]+cff1[k]*np.abs(h);
            z[k,:,:]=z0*h/(h2) + zeta*(1.+z0*h2inv);  
        
    else:
        h[h==0]=1.e-2;
        hinv=1./h;
        cff1=Cs;
        cff2=sc+1;
        cff=hc*(sc-Cs);
        cff2=sc+1;
        for k in range(N):
            z0=cff[k]+cff1[k]*h;
            z[k,:,:]=z0+zeta*(1.+z0*hinv);


    return z