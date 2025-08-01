import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/users/treillou/Diagnostics/Tools')
import spectrum 

# dispersion_relation function for scalar and array inputs
def dispersion_relation(omega, h=10):
    # other method:
    #    khd = h*wf^2/g;  % dispersion relation
    #    kh  = sqrt(   khd*khd + khd/(1.0 + khd*(0.6666666666 ...
    #          +khd*(0.3555555555 + khd*(0.1608465608 ...
    #          +khd*(0.0632098765 + khd*(0.0217540484 ...
    #                             + khd*0.0065407983)))))) ); % Hunt 1979
    #    wki=kh/h;     % wavenumber
    
    g = 9.81  # gravity [m/s^2]
    const = (omega ** 2) * h / g

    # If const is a scalar, convert it to a 1-element array to handle indexing properly
    if np.isscalar(const):
        const = np.array([const])
        scalar_input = True
    else:
        scalar_input = False

    # Initialize wavenumber (kh)
    kh = np.full_like(const, np.nan)

    # Handle special cases
    kh[const == 0] = 0  # Zero const returns zero wavenumber
    positive_const = const > 0

    # Initial guess for Newton-Raphson iteration
    kh[positive_const] = np.sqrt(const[positive_const])

    # Newton-Raphson iteration to solve kh * tanh(kh) = const
    tolerance = 1e-6
    max_iter = 100
    for _ in range(max_iter):
        f = kh[positive_const] * np.tanh(kh[positive_const]) - const[positive_const]
        fprime = kh[positive_const] / np.cosh(kh[positive_const]) ** 2 + np.tanh(kh[positive_const])
        kh[positive_const] -= f / fprime

        # Check for convergence
        if np.max(np.abs(f)) < tolerance:
            break

    # Compute final wavenumber
    k = kh / h

    # Return scalar if input was scalar
    if scalar_input:
        return k[0]
    return k

def jonswap_spectrum(Tp, gamma, Hs, fmin=0.02, fmax=1.0, df=0.001):
    """
    Generate a JONSWAP spectrum Snn(f).
    
    Parameters:
    - Tp : Peak period (s)
    - gamma : Peak enhancement factor
    - Hs : Significant wave height (m)
    - fmin : Minimum frequency (Hz)
    - fmax : Maximum frequency (Hz)
    - df : Frequency resolution (Hz)
    
    Returns:
    - f : Frequency array (Hz)
    - Snn : Spectral density array (m^2/Hz)
    """
    # Constants
    g = 9.81  # gravitational acceleration (m/s²)
    fp = 1.0 / Tp  # peak frequency (Hz)
    alpha = 0.076 * (g**2 / (2 * np.pi)**4) * fp**(-4) * Hs**2  # PM spectrum alpha

    f = np.arange(fmin, fmax + df, df)  # frequency array
    sigma = np.where(f <= fp, 0.07, 0.09)  # spreading parameter

    r = np.exp(- ( (f - fp)**2 ) / (2 * sigma**2 * fp**2))  # peak enhancement shape

    # Pierson-Moskowitz base spectrum
    S_PM = alpha * g**2 * f**(-5) * np.exp(-1.25 * (fp / f)**4)

    # JONSWAP spectrum
    Snn = S_PM * gamma**r

    return f, Snn

def group_velocity(omega,h):
    k = dispersion_relation(omega, h)
    cp = phase_speed(omega,h)
    cg = 0.5 * cp * (1+ (2*k*h)/(np.sinh(2*k*h)) )
    return cg

def phase_speed(omega,h):
    g=9.81
    k = dispersion_relation(omega, h)
    cp = np.sqrt(g/k * np.tanh(k*h))
    return cp


def backrefract_angle(f,thetaF,hF,hWM,fmin=0.06,fmax=0.18):
    imi = np.argmin(np.abs(f-fmin))
    ima = np.argmin(np.abs(f-fmax))
    
    cpF  = phase_speed(2*np.pi*f,hF)
    cpWM = phase_speed(2*np.pi*f,hWM) 
    
    thetaWM = np.arcsin((np.sin(thetaF*np.pi/180) * cpWM / cpF)) * 180/np.pi
    return thetaWM

def backrefract_Snn(f,Snn1,h1,h0,theta1,theta0,fmin=0.06,fmax=0.18):
    imi = np.argmin(np.abs(f-fmin))
    ima = np.argmin(np.abs(f-fmax))
    
    cg1 = group_velocity(2*np.pi*f[imi:ima],h1)
    cg0 = group_velocity(2*np.pi*f[imi:ima],h0)
    Snn0 = (cg0*np.cos(theta0[imi:ima]))/(cg1*np.cos(theta1[imi:ima])) * Snn1[imi:ima]
    return f[imi:ima],Snn0

def backrefract_dirspread(f,Snn1,h1,h0,theta1,theta0,sigma1,fmin=0.06,fmax=0.18):
    imi = np.argmin(np.abs(f-fmin))
    ima = np.argmin(np.abs(f-fmax))
    
    cp1 = phase_speed(2*np.pi*f,h1)
    cp0 = phase_speed(2*np.pi*f,h0)
    
    sigma0 = cp0/cp1 * np.cos(theta1)/np.cos(theta0) * sigma1
    
    return sigma0

def test_oblique_waves(h,Tp,theta):
    g=9.81
    omega=2*np.pi/Tp
    wd=theta*np.pi/180
    
    k = dispersion_relation(omega, h)
    L = 2*np.pi/k    
    Ly = L/np.abs(np.sin(wd))    
    
    a1 = wd*180/np.pi
    a2 = np.arcsin(np.sin(wd)/2)*180/np.pi
    a3 = np.arcsin(np.sin(wd)/3)*180/np.pi
    
    print("Compatible Y-axis length:")
    print("  "+str(np.round(Ly))+"  "+str(np.round(Ly*2))+"  "+str(np.round(Ly*3)))
    print("Compatible incidence angles:")  
    print("  "+str(np.round(a1))+"  "+str(np.round(a2))+"  "+str(np.round(a3)))
    
    
def compute_a1_a2_b1_b2(Ezz, Exx, Eyy, Cxy, Qxz, Qyz):
    """
    Compute directional Fourier components a1, a2, b1, b2.

    Parameters:
    - Ezz: array or value of surface elevation variance density spectrum
    - Exx: array or value of x-displacement variance density spectrum
    - Eyy: array or value of y-displacement variance density spectrum
    - Cxy: array or value of x–y co-spectrum
    - Qxz: array or value of x–z quad-spectrum
    - Qyz: array or value of y–z quad-spectrum

    Returns:
    - a1, a2, b1, b2: arrays or values of directional Fourier components
    """
    denom_a1b1 = np.sqrt(Ezz * (Exx + Eyy))
    denom_a2b2 = Exx + Eyy

    a1 = Qxz / denom_a1b1
    b1 = Qyz / denom_a1b1
    a2 = (Exx - Eyy) / denom_a2b2
    b2 = 2 * Cxy / denom_a2b2

    return a1, a2, b1, b2

def find_theta_tseries(zeta,u,v,fs,fmin=0.05,fmax=0.2,N=256):
    f,Ezz,Qzz=spectrum.compute_cospectrum_quadspectrum(zeta,zeta, fs,N)
    f,Exx,Qzz=spectrum.compute_cospectrum_quadspectrum(u,u, fs,N)
    f,Eyy,Qzz=spectrum.compute_cospectrum_quadspectrum(v,v, fs,N)
    f,Cxy,Qxy=spectrum.compute_cospectrum_quadspectrum(u,v, fs,N)
    f,Cxz,Qxz=spectrum.compute_cospectrum_quadspectrum(u,zeta, fs,N)
    f,Cyz,Qyz=spectrum.compute_cospectrum_quadspectrum(v,zeta, fs,N)
    a1, a2, b1, b2=compute_a1_a2_b1_b2(Ezz, Exx, Eyy, Cxy, Qxz, Qyz)
    theta2=180/np.pi * np.arctan(b2/a2)*0.5
    sigma2=180/np.pi * 0.5* np.sqrt( 1 - a2*np.cos(2*theta2) -b2*np.sin(2*theta2))
    ifmin = np.argmin(np.abs(f-fmin))
    ifmax = np.argmin(np.abs(f-fmax))
    df=f[1]-f[0]
    a1b = np.trapz(a1[ifmin:ifmax]*Ezz[ifmin:ifmax],dx=df) / np.trapz(Ezz[ifmin:ifmax],dx=df)
    b1b = np.trapz(b1[ifmin:ifmax]*Ezz[ifmin:ifmax],dx=df) / np.trapz(Ezz[ifmin:ifmax],dx=df)
    a2b = np.trapz(a2[ifmin:ifmax]*Ezz[ifmin:ifmax],dx=df) / np.trapz(Ezz[ifmin:ifmax],dx=df)
    b2b = np.trapz(b2[ifmin:ifmax]*Ezz[ifmin:ifmax],dx=df) / np.trapz(Ezz[ifmin:ifmax],dx=df)
    theta2b=180/np.pi * np.arctan(b2b/a2b)*0.5
    sigma2b=180/np.pi * np.sqrt(0.5* ( 1 - a2b*np.cos(2*theta2b*np.pi/180) -b2b*np.sin(2*theta2b*np.pi/180)))
    return f,theta2,sigma2,theta2b,sigma2b 
