import numpy as np
import pandas as pd
from netCDF4 import Dataset
import datetime

def v2rho(var_v):
    if len(var_v.shape)==4:
        [T,N,M,Lp]=var_v.shape
        Mp=M+1
        Mm=M-1
        var_rho=np.zeros((T,N,Mp,Lp))
        var_rho[:,:,1:M-1,:]=0.5*(var_v[:,:,0:Mm-1,:]+var_v[:,:,1:M-1,:])
        var_rho[:,:,0,:]=var_rho[:,:,1,:]
        var_rho[:,:,M,:]=var_rho[:,:,M-1,:]
    elif len(var_v.shape)==3:
        [N,M,Lp]=var_v.shape
        Mp=M+1
        Mm=M-1
        var_rho=np.zeros((N,Mp,Lp))
        var_rho[:,1:M-1,:]=0.5*(var_v[:,0:Mm-1,:]+var_v[:,1:M-1,:])
        var_rho[:,0,:]=var_rho[:,1,:]
        var_rho[:,M,:]=var_rho[:,M-1,:]
    elif len(var_v.shape)==2:
        [M,Lp]=var_v.shape
        Mp=M+1
        Mm=M-1
        var_rho=np.zeros((Mp,Lp))
        var_rho[1:M-1,:]=0.5*(var_v[0:Mm-1,:]+var_v[1:M-1,:])
        var_rho[0,:]=var_rho[1,:]
        var_rho[M,:]=var_rho[M-1,:]
    return var_rho