import sys
import numpy as np
from matplotlib import pyplot as plt
import filters as fi

# m^2 * kg^1 * s^(-2) * K^(-1)
kB = 1.38064852e-23

# m^2 * kg^1 * s^(-1)
h = 6.62607004e-34

# m^1 * s^(-1)
c = 299792458

x = 2 * h * c**2

# Units in meters.
solar_r = 6.9599e+8
parsecs = 3.087e+16

"""
Data from
-http://www.astrophysicsspectator.com/topics/observation/MagnitudesAndColors.html
-https://en.wikipedia.org/wiki/Photometric_system

The 0 Magnitude Flux values (in MKS unit system) used below have been derived by
running the program and comparing the result magnitude values with the Flash applet
at http://astro.unl.edu/classaction/animations/light/bbexplorer.html
Using these values results in much lesser difference over the wavelength spectrum
between the magnitude calculated by the two programs than it does when the values
from the above mentioned source are used. Those values are mentioned in the program
as inline comments next to the used values for quick reference and testing purposes.
Perhaps, I am making an error in converting the the source values to the units used
in this program for calculation. It might also be possible that the error in magnitude
is due to the fact that Vega's brightness seen from Earth is variable, or the source
of the values is obselete.
"""

U = 0.0243640  # 0.03980
B = 0.0167443  # 0.06950
V = 0.0110930  # 0.03630
R = 0.0069576  # 0.01196


def integrate(x):
    """
    Return the approximate value of flux of an object in a filter.

    Parameters
    ----------
    x: dict
       Dictionary containing response function values for a filter.

    Attributes
    ----------
    total: float
           Sum of all the response function values.

    Returns
    -------
    flux_: float
           Value of flux in the filter.
    """
    flux_ = 0
    total = 0
    for nm, val in x.items():
        flux_ += flux_wl[nm] * val
        total += val
    flux_ /= total
    return flux_


def flux(T, wl):
    """
    Return the spectral radiance of blackbody.

    Parameters
    ----------
    T: int
       Temperature in kelvin.
    wl: float
       Wavelength in metres.

    Returns
    -------
    B: float
       Spectral radiance in W ^ 1 * sr ^ (-1) * m ^ (-3)

    """
    try:
        a = (h * c) / (wl * kB * T)
        b = np.exp(a) - 1
        d = wl**5
        B = x / (b * d)
        return B
    except ZeroDivisionError:
        return 0


if __name__ == '__main__':
    name = sys.argv[1]
    T = int(name)
    flux_wl = {}  # Stores flux (SI unit), wavelength (nanometres) pairs.

    # Min and max wavelength values in nanometres.
    wl_min = 0
    wl_max = 1100
    wl_nm = wl_min

    # Calculate flux for wavelengths and store in dict.
    while wl_nm <= wl_max:
        wl = wl_nm * 1e-9
        flux_wl[wl_nm] = flux(T, wl)
        wl_nm += 1

    R_ = 1 * solar_r
    r = 10 * parsecs
    ra = (R_ / r)**2

    # Calculate magnitude for filters.
    flux_ = integrate(fi.V)
    mv_V = round(-2.5 * np.log10(flux_ * ra / V), 2)

    flux_ = integrate(fi.U)
    mv_U = round(-2.5 * np.log10(flux_ * ra / U), 2)

    flux_ = integrate(fi.B)
    mv_B = round(-2.5 * np.log10(flux_ * ra / B), 2)

    flux_ = integrate(fi.R)
    mv_R = round(-2.5 * np.log10(flux_ * ra / R), 2)

    u = str(mv_U)
    b = str(mv_B)
    v = str(mv_V)
    r = str(mv_R)
    print("U: " + u + "\nB: " + b + "\nV: " +
          v + "\nR: " + r)

    # Plot flux vs. wavelength graph along with the calculated values and
    # save to PNG file.
    plt.plot(list(flux_wl.keys()), list(flux_wl.values()), 'r.')
    plt.annotate("U: " + u + "\nB: " + b + "\nV: " +
                 v + "\nR: " + r, xy=(0.80, 0.80),
                 xycoords='axes fraction')
    plt.xticks(np.arange(wl_min, wl_max + 1, 100))
    plt.xlabel('Wavelength ($\mathregular{10^{-9}m}$)', fontsize=10)
    plt.ylabel('Flux (W$\mathregular{sr^{-1}}$$\mathregular{m^{-3}}$)',
               fontsize=10)
    plt.title('Flux vs Wavelength graph for a blackbody for constant temperature T\
= ' + str(T) + 'K', fontsize=7)
    plt.show()
    plt.savefig('Flux vs. Wavelength: ' + name + 'K.png')
    plt.close()
