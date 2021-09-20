
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

"""
Sanity check on the nine ASteCA output for the 5 clusters with different
analyzed CMDs and binary fractions
"""

in_folder = '../2_pipeline/XX_ASteCA/out/'
out_folder = '../2_pipeline/XX_ASteCA/tmp/'

data = ascii.read(in_folder + 'asteca_output.dat')

cc = {5: 'orange', 6: 'blue', 7: 'green'}
masses = ['1000', '1100', '1200', '1300', '1400', '1450', '1500', '1600',
          '1700', '1800']
plx_dists = np.array([
    407.4729, 407.3476, 407.26515, 407.22866, 407.1864, 407.05489, 407.16046,
    408.3266, 412.8531, 420.855])


def parPlot(data, par):
    for cl in data:
        run = cl['NAME'].split('/')[0].replace("b_fr_0", "")
        mass = cl['NAME'].split('/')[1].replace("NGC2516_", "")
        mass_idx = masses.index(mass)
        x = mass_idx + 1
        # x = xx[name]
        xoffset = .35 * (float(run) - 6)
        x = x + xoffset
        y = cl['{}_median'.format(par)]
        y_16 = cl['{}_16th'.format(par)]
        y_84 = cl['{}_84th'.format(par)]
        y_mean = cl['{}_mean'.format(par)]

        if par == 'd':
            y = 10**(.2 * (y + 5))
            y_16 = 10**(.2 * (y_16 + 5))
            y_84 = 10**(.2 * (y_84 + 5))
            y_mean = 10**(.2 * (y_mean + 5))

        yerr = np.array([[y - y_16, y_84 - y]]).T
        ax.errorbar(x, y, yerr=yerr, fmt='o', lw=2, color=cc[int(run)])

        # y, yerr = cl['{}_mean'.format(par)], cl['{}_std'.format(par)]
        # ax.errorbar(x - 0.05, y, yerr=yerr, fmt='x', c='k', lw=.7, ls=':')
        plt.scatter(x, y_mean, marker='x', c='k', zorder=4)
        ax.annotate(int(run), (x + .05, y))

        if par == 'd':
            # d_pc = 10**(.2 * (cl['d_median'] + 5))
            # d_pc_16 = 10**(.2 * (cl['d_16th'] + 5))
            # d_pc_84 = 10**(.2 * (cl['d_84th'] + 5))
            # print(mass, d_pc, d_pc_16, d_pc_84)
            ax.plot((mass_idx + .5, mass_idx + 1.5), (
                plx_dists[mass_idx], plx_dists[mass_idx]), ls='--', lw=1.5,
                c='r')


# 'd': r'$\mu$',
ylabels = {'z': 'z', 'a': 'log(age)', 'E': r'$E_{B-V}$)', 'd': 'd [pc]',
           'M': r'M (M$_{\odot}$)', 'b': r'b$_{fr}$'}
for par in ('z', 'a', 'E', 'd', 'M', 'b'):
    print(par)
    fig, ax = plt.subplots(figsize=(20, 10))
    parPlot(data, par)
    # ax.grid(ls=':')
    ax.set_xticks(np.arange(len(masses)) + 1)
    ax.set_xticklabels(masses)
    # ax.set_yscale('log')
    ax.axvline(1.5, ls=':', lw=1.5, c='k')
    ax.axvline(2.5, ls=':', lw=1.5, c='k')
    ax.axvline(3.5, ls=':', lw=1.5, c='k')
    ax.axvline(4.5, ls=':', lw=1.5, c='k')
    ax.axvline(5.5, ls=':', lw=1.5, c='k')
    ax.axvline(6.5, ls=':', lw=1.5, c='k')
    ax.axvline(7.5, ls=':', lw=1.5, c='k')
    ax.axvline(8.5, ls=':', lw=1.5, c='k')
    ax.axvline(9.5, ls=':', lw=1.5, c='k')
    plt.xlim(.5, 10.5)
    plt.ylabel(ylabels[par], fontsize=12)
    plt.xlabel("Number of members", fontsize=12)
    plt.savefig(
        out_folder + "{}.png".format(par), dpi=150, bbox_inches='tight')
