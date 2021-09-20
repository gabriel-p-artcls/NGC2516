
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval


fname = "NGC2516_VOR.dat"
in_file_folder = '../2_pipeline/6_pyUPMASK/out/'
out_fig = '../2_pipeline/paper_figs/NGC2516_4_plots.png'

# plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12


print("Reading data")
data = Table.read(in_file_folder + fname, format='ascii')
print("Total number of stars: {}".format(len(data)))

# Plot
gs_unit = 5
gs_x, gs_y = 4, 1
fig, axs = plt.subplots(gs_y, gs_x, figsize=(
    gs_unit * gs_x, gs_unit * gs_y))

N = len(data)
interval = ZScaleInterval()
zmin, zmax = interval.get_limits(data['Gmag'])
mag = np.clip(data['Gmag'], a_min=zmin, a_max=zmax)
factor = 5 * 500. * (1 - 1 / (1 + 150 / N ** 0.85))
sizes = .1 + factor * (10 ** ((mag - zmin) / -2.5))

ax = axs[0]
ax.set_title("A")
ax.minorticks_on()
ax.scatter(data['_x'], data['_y'], alpha=.5, marker='.',
           s=sizes, color='grey')
ax.set_xlabel(r"$\alpha^{*}$")
ax.set_ylabel(r"$\delta^{*}$")
ax.set_ylim(-2, 2)
ax.set_xlim(2, -2)
ax.set_xticks((-2, -1, 0, 1, 2))
ax.set_yticks((-2, -1, 0, 1, 2))

#
ax = axs[1]
ax.set_title("B")
ax.minorticks_on()
ax.scatter(
    data['BP-RP'], data['Gmag'], marker='x', s=10, lw=.25, alpha=.5)
# edgecolor='grey',            facecolor='none', s=10, alpha=.5)
ax.set_ylim(18.05, 4)
ax.set_xlim(-.2, 3.2)
ax.set_xlabel("BP-RP")
ax.set_ylabel("G")

#
ax = axs[2]
ax.set_title("C")
ax.minorticks_on()
cx, cy, length = -4.8, 11, 10
xmax, xmin = cx + length, cx - length
ymax, ymin = cy + length, cy - length
msk = (data['pmRA'] < xmax) & (data['pmRA'] > xmin) &\
    (data['pmDE'] < ymax) & (data['pmDE'] > ymin)
ax.scatter(
    data['pmRA'][msk], data['pmDE'][msk], alpha=.2, marker='x',
    lw=.1)
ax.set_xlabel(r"$\mu_{\alpha} cos\delta$ [mas/yr]")
ax.set_ylabel(r"$\mu_{\delta}$ [mas/yr]")
ax.set_xlim(xmax, xmin)
ax.set_ylim(ymin, ymax)

#
ax = axs[3]
ax.set_title("D")
pmin, pmax = -.3, 3.3
msk = (data['Plx'] > pmin) & (data['Plx'] < pmax)
ax.hist(data['Plx'][msk], 50, alpha=.75)
ax.set_xlim(pmin, pmax)
ax.set_xlabel("Plx [mas]")
ax.set_ylabel("N")

fig.tight_layout()
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
