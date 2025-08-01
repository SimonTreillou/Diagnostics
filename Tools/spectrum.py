from scipy.signal import welch, hanning, csd
import numpy as np

def spectrum(var, dt,N=256,along=False):
    nfft=N
    noverlap=N/2
    win = hanning(nfft, True)
    if along:
        S=0
        for i in range(var.shape[1]):
            f, S_tmp = welch(var[:,i], 1/dt, window=win, noverlap=noverlap, nfft=nfft, return_onesided=True)
            S=S+S_tmp
        S=S/var.shape[1]
    else:
        f, S = welch(var, 1/dt, window=win, noverlap=noverlap, nfft=nfft, return_onesided=True)
    return f,S

def moment(f,S,k):
    m = np.trapz(f**k*S)
    return m

def compute_cospectrum_quadspectrum(x, z, fs, nperseg=None):
    """
    Compute co-spectrum and quad-spectrum between signals x and z.

    Parameters:
    - x, z: time series arrays
    - fs: sampling frequency (Hz)
    - nperseg: length of each segment for Welch method (default: auto)

    Returns:
    - f: array of frequency bins
    - Cxz: co-spectrum (real part of cross-spectral density)
    - Qxz: quad-spectrum (imaginary part of cross-spectral density)
    """
    f, Pxy = csd(x, z, fs=fs, nperseg=nperseg)
    Cxz = np.real(Pxy)
    Qxz = np.imag(Pxy)
    return f, Cxz, Qxz
