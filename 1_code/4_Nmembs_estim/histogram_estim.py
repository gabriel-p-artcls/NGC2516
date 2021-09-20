
from astropy.io import ascii
from scipy.spatial.distance import cdist
from scipy.stats import median_abs_deviation as MAD
import numpy as np
from scipy.stats import binned_statistic_2d as bs2d
from itertools import compress
import matplotlib.pyplot as plt


"""
1. Generate a 2D histogram for the cluster region in PMs
2. Generate a 2D histogram for a field region, using the same edges as above
3. Subtract both histograms to see which cells contain an excess of cluster
   stars
4. For those cells, find all the cluster stars within them
5. Select a random set of N cluster stars in each of those cells, where N is
   the number in the cell
6. Convert the number of times a given star survived the elimination process
   into a probability

Gives reasonable results but it depends on a manually selected probability
value to isolate the final cluster members, which invalidates the
goal of being user-independent.
"""

fname = "NGC2516_7deg.dat"
in_folder = '../../2_pipeline/1_data_filter/out/'
out_fig_folder = '../../2_pipeline/4_members_estim/tmp/'
# Fixed center and radius defining the cluster region
cent, rad = (0., 0.), 1.5

# Probability value to filter the cluster regions that survived the most
prob_max = .5

print("Reading data")
data = ascii.read(in_folder + fname)
print("Total number of stars:", len(data))

# Remove nans from PMs
msk = (data['pmRA'].mask) | (data['pmDE'].mask)
data = data[~msk]
print("Stars filtered by nans in PMs:", len(data))

MAD_pm = np.sqrt(MAD(data['pmRA'])**2 + MAD(data['pmDE'])**2)
pm = np.array([data['pmRA'], data['pmDE']]).T
pm_cent = (np.median(data['pmRA']), np.median(data['pmDE']))
pm_dists = cdist([pm_cent], pm)[0]
msk = pm_dists < 5 * MAD_pm
data = data[msk]
print("Stars filtered by large PMs:", len(data))

length = min(np.ptp(data['_x']), np.ptp(data['_y']))
# Maximum number of rings to use
i_ring = int((length / (2 * rad))**2)


def main(bins_min=10, bins_max=50, bins_step=10, Nruns=2):
    """
    """
    x, y = data['_x'], data['_y']
    xy = np.array([x, y]).T
    xy_dists = cdist([cent], xy)[0]
    msk_cl = xy_dists < rad
    clust_reg = np.array([data[_][msk_cl] for _ in ('pmRA', 'pmDE')])
    e_cl_reg = np.array([data[_][msk_cl] for _ in ('e_pmRA', 'e_pmDE')])
    Ncl = msk_cl.sum()
    print("Stars in cluster region: {}".format(Ncl))

    idx_survive, Ntot, rad_old = [], 0, rad
    # membs_lst
    for ring in range(2, i_ring):
        rad_ring = np.sqrt(ring) * rad
        print(rad_ring)
        msk_fl = (xy_dists > rad_old) & (xy_dists <= rad_ring)
        rad_old = rad_ring
        field_reg = np.array([data[_][msk_fl] for _ in ('pmRA', 'pmDE')])
        e_fl_reg = np.array([data[_][msk_fl] for _ in ('e_pmRA', 'e_pmDE')])
        Nfl = msk_fl.sum()

        for bins in range(bins_min, bins_max, bins_step):
            print(bins)
            for run in range(Nruns):
                # Randomly perturb the data according to their uncertainties
                cl_r = clust_reg + e_cl_reg * np.random.normal(
                    0., 1., (2, Ncl))
                fl_r = field_reg + e_fl_reg * np.random.normal(
                    0., 1., (2, Nfl))

                cl_res = bs2d(*cl_r, None, 'count', bins=bins,
                              expand_binnumbers=True)
                fr_res = np.histogram2d(
                    *fl_r, bins=(cl_res.x_edge, cl_res.y_edge))[0]

                cl_ids = cl_res.binnumber
                diff_hist = cl_res.statistic - fr_res
                # membs_lst.append(np.clip(diff_hist, a_min=0, a_max=None).sum())
                # membs_lst.append(diff_hist.sum())

                # # Find (x, y) coordinates of the bin with the maximum value
                # h_max = np.unravel_index(diff_hist.argmax(), diff_hist.shape)
                # xmax = cl_res.x_edge[h_max[0]:h_max[0] + 2].mean()
                # ymax = cl_res.y_edge[h_max[1]:h_max[1] + 2].mean()

                idx_survive_h = []
                for i, xi in enumerate(diff_hist):
                    for j, Nc in enumerate(xi):
                        if Nc > 0:
                            # All the stars within this cell
                            mxy = (cl_ids[0] == i + 1) & (cl_ids[1] == j + 1)
                            # Store indexes of cluster region stars that
                            # survived the process of removing the excess
                            # field stars
                            idxs_in_cell = list(compress(range(len(mxy)), mxy))
                            # Select random 'Nc' indexes
                            if len(idxs_in_cell) > Nc:
                                idx_survive_h += list(np.random.choice(
                                    idxs_in_cell, int(Nc), replace=False))
                            else:
                                idx_survive_h += idxs_in_cell

                idx_survive += idx_survive_h
                Ntot += 1

    # Extract the unique indexes of cluster region stars that survived at
    # least once, and the total number of times that they survived
    idxs, Nc = np.unique(idx_survive, return_counts=True)
    # Convert to probability
    probs_final = np.zeros(clust_reg.shape[1])
    probs_final[idxs] = Nc
    probs_final /= Ntot

    # Filter stars outside the general region of the PMs overdensity
    xc, yc, radPM = -4.6, 11.3, 2
    msk = (abs(clust_reg[0, :] - xc) <= radPM) &\
        (abs(clust_reg[1, :] - yc) <= radPM)

    plt.subplot(121)
    plt.hist(probs_final[msk], 50)
    # plt.subplot(122)
    # plt.hist(membs_lst, 50)
    # plt.show()
    plt.subplot(122)
    # msk = probs_final > prob_max
    # print(msk.sum())
    plt.scatter(*clust_reg[:, msk], c=probs_final[msk])
    plt.show()


if __name__ == '__main__':
    main()
