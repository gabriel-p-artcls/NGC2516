
import numpy as np
import matplotlib.pyplot as plt


"""
Calculate approximate number of cluster members for a given King Profile,
from 0 up to rt. General form of Eq (3) in Froebrich et al. (2007);
374, 399-408.
"""

# Maximum density value, taken from the A2 plot
max_dens = 3.2 * 3600

# Core and tidal radius for the KP
rt = 90 / 60
for rc in (10, 12, 14, 16, 18):
    rc = rc / 60

    xy = []
    for fd in np.arange(6500, 8000, 100):

        # normalizing constant
        cd = (max_dens - fd) / (rc**2 * ((rc**2)**-.5 - (rc**2 + rt**2)**-.5)**2)

        x = 1 + (rt / rc) ** 2
        N_memb = int(round(
            (np.pi * cd * rc ** 2) * (
                np.log(x) - 4 + (4 * np.sqrt(x) + (x - 1)) / x)))

        xy.append([fd, N_memb])

    plt.title("rc={:.1f}".format(rc * 60))
    plt.plot(*np.array(xy).T)
    plt.axvline(7150, ls=":")
    plt.show()
