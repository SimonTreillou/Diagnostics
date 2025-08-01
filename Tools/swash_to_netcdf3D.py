import numpy as np
import scipy.io
import re
from netCDF4 import Dataset
import os

def swash_to_netcdf3D(source, dir_path, var):
    filepath = os.path.join(source, dir_path, f"{var}.mat")
    data = scipy.io.loadmat(filepath)

    if var == 'u':
        swashname = 'Vksi'
    elif var == 'v':
        swashname = 'Veta'
    else:
        raise ValueError("Unsupported variable. Use 'u' or 'v'.")

    times = []
    layers = set()
    parsedData = {}

    time_map = {}

    # Parse all variable names
    for varName in data:
        match = re.match(rf"{swashname}_k(\d+)_([0-9]{{6}})_([0-9]{{3}})$", varName)
        if match:
            k = int(match.group(1))
            hhmmss = match.group(2)
            subsec = match.group(3)

            h = int(hhmmss[0:2])
            m = int(hhmmss[2:4])
            s = int(hhmmss[4:6])
            micro = int(subsec)

            t_seconds = h * 3600 + m * 60 + s + micro / 1000.0
            times.append(t_seconds)

            time_map.setdefault(t_seconds, []).append((k, varName))
            layers.add(k)

    timesSorted = sorted(set(times))
    layersSorted = sorted(layers)

    time_index = {t: i for i, t in enumerate(timesSorted)}
    layer_index = {k: i for i, k in enumerate(layersSorted)}

    # Get dimensions from a sample variable
    # Get dimensions from a sample 2D variable
    sample_data = None
    for v in data.values():
        if isinstance(v, np.ndarray) and v.ndim == 2 and np.issubdtype(v.dtype, np.number):
            sample_data = v
            break

    if sample_data is None:
        raise ValueError("No valid 2D numerical variable found in .mat file.")

    ySize, xSize = sample_data.shape
    
    # Initialize data array with NaNs
    output_data = np.full((xSize, ySize, len(layersSorted), len(timesSorted)), np.nan)

    # Fill the data array
    for t in timesSorted:
        for k, varName in time_map[t]:
            try:
                value = data[varName]
                i = layer_index[k]
                j = time_index[t]
                output_data[:, :, i, j] = value.T
            except KeyError:
                continue
    
    if hasattr(data, 'Xp'):
        Xp=data['Xp'][:]
    if hasattr(data, 'Botlev'):
        h=data['Botlev'][:]

    # Write to NetCDF
    ncFileName = os.path.join(source, dir_path, f"swash_{var}.nc")
    with Dataset(ncFileName, 'w', format='NETCDF4') as ncfile:
        ncfile.createDimension('x', xSize)
        ncfile.createDimension('y', ySize)
        ncfile.createDimension('layer', len(layersSorted))
        ncfile.createDimension('time', len(timesSorted))

        var_nc = ncfile.createVariable(var, 'f8', ('x', 'y', 'layer', 'time'))
        x_nc = ncfile.createVariable('Xp', 'f8', ('x','y'))
        h_nc = ncfile.createVariable('Botlev', 'f8', ('x','y'))
        time_nc = ncfile.createVariable('time', 'f8', ('time',))
        layer_nc = ncfile.createVariable('layer', 'i4', ('layer',))

        var_nc[:, :, :, :] = output_data
        if hasattr(data, 'Xp'):
            x_nc[:] = Xp
        else:
            x_nc[:] = np.tile(np.linspace(0,xSize,xSize),[1,ySize])
        if hasattr(data, 'Botlev'):
            h_nc[:] = h
        else:
            h_nc[:] = np.zeros_like(x_nc)
        time_nc[:] = timesSorted
        layer_nc[:] = layersSorted

    print(f"NetCDF file created: {ncFileName}")