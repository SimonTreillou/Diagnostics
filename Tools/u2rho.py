import numpy as np
import pandas as pd
from netCDF4 import Dataset
import datetime

def u2rho(var_u):
    if len(var_u.shape)==4:
        [T,N,Mp,L]=var_u.shape
        Lp=L+1
        Lm=L-1
        var_rho=np.zeros((T,N,Mp,Lp))
        var_rho[:,:,:,1:L-1]=0.5*(var_u[:,:,:,0:Lm-1]+var_u[:,:,:,1:L-1])
        var_rho[:,:,:,0]=var_rho[:,:,:,1]
        var_rho[:,:,:,L]=var_rho[:,:,:,L-1]
    elif len(var_u.shape)==3:
        [N,Mp,L]=var_u.shape
        Lp=L+1
        Lm=L-1
        var_rho=np.zeros((N,Mp,Lp))
        var_rho[:,:,1:L-1]=0.5*(var_u[:,:,0:Lm-1]+var_u[:,:,1:L-1])
        var_rho[:,:,0]=var_rho[:,:,1]
        var_rho[:,:,L]=var_rho[:,:,L-1]
    elif len(var_u.shape)==2:
        [Mp,L]=var_u.shape
        Lp=L+1
        Lm=L-1
        var_rho=np.zeros((Mp,Lp))
        var_rho[:,1:L-1]=0.5*(var_u[:,0:Lm-1]+var_u[:,1:L-1])
        var_rho[:,0]=var_rho[:,1]
        var_rho[:,L]=var_rho[:,L-1]
    return var_rho