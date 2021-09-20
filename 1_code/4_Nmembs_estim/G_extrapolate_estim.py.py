
import numpy as np
from astropy.io import ascii
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

"""
Use the King's profile estimation of total members up to G=15 mag, and
extrapolate this number to the G=18 mag limit.

This should give an estimate of the total number of members (up to G=18 mag)
bypassing the non-uniform field contamination below G=15 mag.

The result is almost 2100 members, which is a lot more than other methods
estimate
"""

# Members estimated up to the above tidal radius
N_membs_KP = 330

data = ascii.read("../../2_pipeline/1_data_filter/out/NGC2516.dat")

cent = np.array([[0.01], [-0.031]]).T
xy = np.array([data[_] for _ in ('_x', '_y')]).T
dist = cdist(cent, xy)

# Tidal radius
rt_max = 2

rads = np.linspace(.1, rt_max, 10)
N_membs_rad, prop = [], []
for rad in rads:

    # Stars within this radius
    msk_rad = dist[0] <= rad
    data_r = data[msk_rad]

    # Proportion of all stars to G<15 stars within this radius
    N_G15 = (data_r['Gmag'] <= 15.).sum()
    proportion = len(data_r) / N_G15

    # Total number of members within this radius
    N_membs_rad.append(proportion * N_membs_KP)
    prop.append(proportion)

plt.subplot(121)
plt.plot(rads, prop)
plt.subplot(122)
plt.title("N_tot = {:.0f}".format(N_membs_rad[-1]))
plt.plot(rads, N_membs_rad)
plt.show()
