
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
import astropy.units as u

G = 4.3e-6  # (kpc/M_solar) * (km/s)^2

# Pinfield et al. (1998), 'The mass of the Pleiades'
# https://ui.adsabs.harvard.edu/abs/1998MNRAS.299..955P/abstract

# Wikipedia
A, B = 15.3, -11.9  # (km/s)/kpc

# NGC2516
ra_N2516, de_N2516, d_NGC2516 = 119.51667, -60.75333, 400
c1 = coord.SkyCoord(ra=ra_N2516*u.degree, dec=de_N2516*u.degree,
                    distance=d_NGC2516*u.pc, frame='icrs')
gc1 = c1.transform_to(coord.Galactocentric)
RG = np.sqrt(gc1.x**2 + gc1.y**2)
# Piskunov et al. 2007
A0, B0, RG0 = 14.5, -13, 8.5*u.kpc
delta_Rg = (RG - RG0) / RG0
A = A0 - A0 * delta_Rg
B = A - A0 + B0 + 2*A0*delta_Rg

#
rt = np.arange(1, 25, .1)  # pc
rt_kpc = rt / 1000.

# M = (rt_kpc**3 * 2 * (A - B)**2) / G

# Piskunov et al. 2007
M = (4 * A * (A - B) * rt_kpc**3) / G

# Pinfield et al. (1998)
M2 = (1000 * rt_kpc / 1.46) ** 3

plt.title(r"Cluster mass from $r_{t}$ and Oort's constants")
plt.plot(rt, M)
plt.plot(rt, M2)
plt.xlabel("Tidal radius [pc]")
plt.ylabel(r"Mass $[M_{\odot}]$")
plt.grid(ls=":")
plt.show()
