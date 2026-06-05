import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uncertainties as un
import scipy
from scipy.optimize import curve_fit
from scipy import constants
from uncertainties import ufloat as un
import matplotlib as mpl
#import uncertainties.umath as umath
from pylab import cm, figure, text, scatter, show
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import scipy.signal as sig
import scipy.constants as const


csfont = {'fontname':'Times New Roman'}
hfont = {'fontname':'Times New Roman'}


#------------------------------------------------------------------------------------------------
#####################
#####Fitten##########
#####################

#Erstellen der Datensets für den Fitplot
def daten_fitten(funct, x_data, y_data, genauigkeit, p, vl=0, vr=0, fl=0, fr=0, bounds=[-np.inf, np.inf], intervall=[0, 0], sigma=None, absolute_sigma=False): #Step4/5

    #für chi_squared
    if sigma == None:
        y_err = np.sqrt(np.abs(y_data))
    else:
        y_err = sigma

    #Fitintervall kürzen
    if (fl != 0 or fr != 0):
        n = len(x_data)
        start_idx = int(fl * n) if fl > 0 else 0
        end_idx = n - int(fr * n) if fr > 0 else n
        if start_idx < end_idx:
            x_data = x_data[start_idx:end_idx]
            y_data = y_data[start_idx:end_idx]
            if sigma != None:
                sigma = sigma[start_idx:end_idx]

    #fiten
    pars, cov = curve_fit(f=funct, xdata=x_data, ydata=y_data, p0=p, bounds=bounds, maxfev=500000, sigma=sigma, absolute_sigma=absolute_sigma)
    # Get the standard deviations of the parameters (square roots of the diagonal of the covariance)
    stdevs = np.sqrt(np.diag(cov))
    #print('Fitparameter')
    #print(pars)
    #print(stdevs)
    #datensets
    if intervall == [0, 0]:
        a = np.abs(max(x_data)-min(x_data))
    else:
        a = intervall[1]-intervall[0]
    xfit = np.linspace(min(x_data)-a*vl, max(x_data)+a*vr, genauigkeit)

    yfit = funct(xfit, *pars)
    ychi2 = funct(np.array(x_data), *pars)

    #chi_sq = chi_squared(funct, pars, ychi2, y_data, y_err)


    return pars, stdevs, xfit, yfit



#BErechnung von Chi_Sqared
def chi_squared(funct, pars, y_fitted, y_data, y_err):
    # Berechnung der gefitteten Werte
    chi_sq = 0
    n = len(y_data)

    for i in range(n):
        if(y_err[i] != 0):
            chi_sq += (y_data[i] - y_fitted[i])**2 / y_err[i]**2

    print('Chi^2: ', chi_sq)
    print('red Chi^2: ', chi_sq/(len(y_fitted)-len(pars)))
    return chi_sq

def lnOfList(A):
    n = len(A)
    C = [0]*n
    for i in range(n):
        C[i] = np.log(A[i])
    return C

def linear_funct(x, m, c):
    return -m*x + c

def mass_gap_finite_size(N, m, A, xi):
    N = np.asarray(N, dtype=float)
    return m + A / np.sqrt(N) * np.exp(-N / xi)

def fit_mass_gap_convergence(N_data, gap_data, fl=0, fr=0):
    N_data = np.asarray(N_data, dtype=float)
    gap_data = np.asarray(gap_data, dtype=float)

    if len(N_data) < 3:
        raise ValueError("At least three N values are needed to fit m, A, and xi")

    m0 = max(float(gap_data[-1]), 1e-12)
    A0 = max(float(gap_data[0] - m0), 0.36)
    xi0 = max((float(max(N_data)) - float(min(N_data))) / 2.0, 1.0)
    #Fitintervall kürzen
    if (fl != 0 or fr != 0):
        n = len(N_data)
        start_idx = int(fl * n) if fl > 0 else 0
        end_idx = n - int(fr * n) if fr > 0 else n
        if start_idx < end_idx:
            x_data = N_data[start_idx:end_idx]
            y_data = gap_data[start_idx:end_idx]

    bounds = ([0.0, 0, 1e-12], [np.inf, np.inf, np.inf])
    pars, cov = curve_fit(
        f=mass_gap_finite_size,
        xdata=N_data,
        ydata=gap_data,
        p0=[m0, A0, xi0],
        bounds=bounds,
        maxfev=500000,
    )
    stdevs = np.sqrt(np.diag(cov))

    xfit = np.linspace(min(N_data), max(N_data), 200)
    yfit = mass_gap_finite_size(xfit, *pars)

    return pars, stdevs, xfit, yfit
