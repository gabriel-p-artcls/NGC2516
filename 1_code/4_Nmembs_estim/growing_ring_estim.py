
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

"""
A simple way to estimate the number of members. Starting from fixed center
values in xy+PMs+Plx, we count how many stars are simultaneously included
in growing 4D spheres. This is extended up to the estimated radius of the
cluster, i.e. 1.5 deg.
"""


fname = "NGC2516_7deg.dat"
in_folder = '../2_pipeline/1_data_filter/out/'
out_fig_folder = '../2_pipeline/5_members_estim/tmp/'

print("Reading data")
data = ascii.read(in_folder + fname)
print("Total number of stars:", len(data))

# Initial centers
RA_m, DE_m = 0.008489, -0.0408002
pmRA_m, pmDE_m = -4.65389, 11.225667
Plx_m = 2.431259

# CMD filter (points on either side of a straight line)
# Source: https://stackoverflow.com/a/3838398/1391441
x1, x2, y1, y2 = 0.8, 2.5, 13.5, 18.

# The radius in xy+PMs+Plx will increase by these step in each corresponding
# unit. The Plx radius is divided by 5 and has a maximum at 0.2 mas
rad = 0.1

N_membs_rad = []
for i in range(1, 16):
    rad_i = i * rad
    print(rad_i, rad_i / 5)

    # Coords, PMs, Plx filters
    msk1 = np.sqrt((data['_x'] - RA_m)**2 + (data['_y'] - DE_m)**2) < rad_i
    msk2 = np.sqrt(
        (data['pmRA'] - pmRA_m)**2 + (data['pmDE'] - pmDE_m)**2) < rad_i
    msk3 = abs(data['Plx'] - Plx_m) < min(rad_i / 5, 0.2)

    # CMD filter
    v1 = (x2 - x1, y2 - y1)
    v2 = (x2 - data['BP-RP'], y2 - data['Gmag'])
    xp = v1[0] * v2[1] - v1[1] * v2[0]  # Cross product
    msk4 = xp > 0.

    # Apply filters
    msk = msk1 & msk2 & msk3 & msk4
    N_membs_rad.append([rad_i, msk.sum()])

N_membs_rad = np.array(N_membs_rad).T

fig = plt.figure(figsize=(10, 5))
plt.title("Nmemb={}".format(N_membs_rad[1][-1]))
plt.plot(*N_membs_rad)
plt.xlabel("Rad [xy/PMs; 5*Plx]")
plt.ylabel("N_membs")
fig.tight_layout()
fout = out_fig_folder + fname.replace('.dat', '_2.png')
plt.savefig(fout, dpi=150, bbox_inches='tight')
