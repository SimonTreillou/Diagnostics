import numpy as np
import scipy.io
import re
from netCDF4 import Dataset
import os

def swash_to_netcdfdomain(source, dir_path, var):
    filepath = os.path.join(source, dir_path, f"{var}.mat")
    data = scipy.io.loadmat(filepath)

    if var == 'domain':
        swashname = 'Botlev'
    else:
        raise ValueError("Only 'zeta' is supported currently.")

    Xp=data['Xp'][:]
    h=data['Botlev'][:]
    print(h.shape)
    xSize, ySize = h.shape
    
    # Write to NetCDF
    ncFileName = os.path.join(source, dir_path, f"swash_{var}.nc")
    with Dataset(ncFileName, 'w', format='NETCDF4') as ncfile:
        ncfile.createDimension('x', xSize)
        ncfile.createDimension('y', ySize)
        
        x_nc = ncfile.createVariable('Xp', 'f8', ('x','y'))
        h_nc = ncfile.createVariable('Botlev', 'f8', ('x','y'))

        x_nc[:] = Xp
        h_nc[:] = h

    print(f"NetCDF file created: {ncFileName}")