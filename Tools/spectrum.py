from scipy.signal import welch, csd
from scipy.signal.windows import hann
import numpy as np
import numpy as np
from scipy import signal
from scipy.stats import chi2
import matplotlib.pyplot as plt

def spectrum(var, dt,N=256,along=False):
    nfft=N
    noverlap=N/2
    win = hann(nfft, True)
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


def welch_spectrum_CI(x, fs, nperseg=256, noverlap=None, alpha=0.05):
    """
    Compute power spectrum with DOF and confidence intervals using Welch's method.

    Parameters
    ----------
    x : array-like
        Time series data.
    fs : float
        Sampling frequency in Hz.
    nperseg : int, optional
        Length of each Welch segment. Default 256.
    noverlap : int, optional
        Number of overlapping points. Default is nperseg//2.
    alpha : float, optional
        Significance level for CI (0.05 = 95% CI).

    Returns
    -------
    f : ndarray
        Frequency array.
    Pxx : ndarray
        Power spectral density.
    dof : int
        Degrees of freedom for the PSD estimate.
    ci_lower : ndarray
        Lower bound of confidence interval.
    ci_upper : ndarray
        Upper bound of confidence interval.
    """

    if noverlap is None:
        noverlap = nperseg // 2  # default like scipy

    # --- Welch PSD ---
    f, Pxx = signal.welch(x, fs=fs, nperseg=nperseg, noverlap=noverlap)

    # --- Compute DOF ---
    step = nperseg - noverlap
    n_segments = (len(x) - noverlap) // step
    dof = 2 * n_segments

    # --- Confidence intervals ---
    ci_lower = dof / chi2.ppf(1 - alpha/2, dof)
    ci_upper = dof / chi2.ppf(alpha/2, dof)

    return f, Pxx, dof, ci_lower, ci_upper


def plot_CIbar_loglog(f,S,ci_upper,color='k',x_pos=None,y_pos=None,lw=3):
    if x_pos==None:
        x_pos = f[-1] * 0.8  # near top-right, 80% of max freq
    if y_pos==None:
        y_center = max(S) * 0.8  # position around the top of PSD curve
    else:
        y_center=y_pos
    # CI width (half height in log10 units)
    ci_half = ci_upper * S[0]  # same everywhere

    # Draw the bar (vertical line)
    plt.vlines(x=x_pos, ymin=y_center - ci_half, ymax=y_center + ci_half,
               colors=color, linewidth=lw)

    # Add small horizontal “caps”
    plt.hlines([y_center - ci_half, y_center + ci_half],
               x_pos * 0.95, x_pos * 1.05, colors=color, linewidth=lw)