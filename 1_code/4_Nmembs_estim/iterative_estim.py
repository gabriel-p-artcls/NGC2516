
from astropy.io import ascii
import numpy as np
# from scipy.spatial.distance import cdist
from astropy.table import vstack
from scipy.stats import median_abs_deviation as MAD
import matplotlib.pyplot as plt

"""
Estimate the number of members by filtering stars beyond the approximated
cluster region in all dimensions.

Starting from reasonable approximations to the center coordinates in
(x, y, Plx, pmRA, pmDE), a small radius in xy, and large radii in PMs and Plx,
the algorithm proceeds as follows:

1. Filter stars outside the region defined by the centers+radii. Also apply
   a manually estimated CMD filter
2. Obtain new center estimations as the means
3. Check if these new centers deviate lees than 'tol' from the old ones
4. If they don't keep these new centers, and estimate new radii using the
   'Nstd' parameter
5. Repeat 1-4 until the condition is met

I've found that even using a small xy radius and large Plx+PMs radii, the
convergence is very reasonable for Nstd=3. For Nstd=2 the visibly few stars
are selected as members. For Nstd=4 there is a clear dispersion in the selected
members in the xy space.

Thus, Nstd=3 seems to be the more reasonable choice.

"""

fname = "NGC2516_15deg.dat"
in_folder = '../../2_pipeline/vizier_data_format/'
out_fig_folder = '../../2_pipeline/4_members_estim/tmp/'
out_dat_folder = '../../2_pipeline/4_members_estim/out/'

print("Reading data")
data = ascii.read(in_folder + fname)
print("Total number of stars:", len(data))

# Initial centers
RA_m, DE_m = 0., 0.
pmRA_m, pmDE_m = -4.5, 11.25
Plx_m = 2.4

# CMD filter (points on either side of a straight line)
# Source: https://stackoverflow.com/a/3838398/1391441
x1, x2, y1, y2 = 0.8, 2.8, 13.5, 19.

# Initial small xy radius
xyRad = 5
# Initial PMs & Plx radii
pmRad, Plx_rad = 5, 1

# Tolerance for the center estimation
tol = 0.001

# Number of standard deviations to use in the filter
Nstd = 3.5

full_IDs = [int(_['EDR3Name'].split(' ')[-1]) for _ in data]
# Iterative block to refine the centers and radii
while True:
    # # Plx filter
    # # msk = abs(data['Plx'] - Plx_m) < Nstd * data['e_Plx']
    # msk = abs(data['Plx'] - Plx_m) < Plx_rad
    # members, field = data[msk], data[~msk]
    # # PMs filter
    # msk = np.sqrt(
    #     (members['pmRA'] - pmRA_m)**2 + (members['pmDE'] - pmDE_m)**2)\
    #     < pmRad
    # field = vstack([field, members[~msk]])
    # members = members[msk]
    # # CMD filter
    # v1 = (x2 - x1, y2 - y1)
    # v2 = (x2 - members['BP-RP'], y2 - members['Gmag'])
    # xp = v1[0] * v2[1] - v1[1] * v2[0]  # Cross product
    # msk = xp > 0.
    # field = vstack([field, members[~msk]])
    # members = members[msk]

    # msk1 = abs(data['Plx'] - Plx_m) < (Nstd + 1) * data['e_Plx']
    msk1 = abs(data['Plx'] - Plx_m) < Plx_rad
    msk2 = np.sqrt(
        (data['pmRA'] - pmRA_m)**2 + (data['pmDE'] - pmDE_m)**2) < pmRad
    # CMD filter
    v1 = (x2 - x1, y2 - y1)
    v2 = (x2 - data['BP-RP'], y2 - data['Gmag'])
    xp = v1[0] * v2[1] - v1[1] * v2[0]  # Cross product
    msk3 = xp > 0.

    msk4 = np.sqrt((data['_x'] - RA_m)**2 + (data['_y'] - DE_m)**2) < xyRad

    msk = msk1 & msk2 & msk3 & msk4
    members, field = data[msk], data[~msk]

    RA_n, DE_n = np.mean(members['_x']), np.mean(members['_y'])
    pmRA_n, pmDE_n = np.mean(members['pmRA']), np.mean(members['pmDE'])
    Plx_n = np.mean(members['Plx'])
    print(("x_c={:.2f}, y_c={:.2f}, pmRA_c={:.2f}, pmDE_c={:.2f}, "
           + "Plx_c={:.2f}").format(RA_n, DE_n, pmRA_n, pmDE_n, Plx_n))

    # Tolerance check
    if abs(1 - RA_n / RA_m < tol) and abs(1 - DE_n / DE_m < tol) and\
            abs(1 - pmRA_n / pmRA_m < tol) and\
            abs(1 - pmDE_n / pmDE_m < tol) and abs(1 - Plx_n / Plx_m < tol):
        break
    else:
        print("New iteration")

    # Update center values
    RA_m, DE_m = RA_n, DE_n
    pmRA_m, pmDE_m = pmRA_n, pmDE_n
    Plx_m = Plx_n

    # Update radii
    # xyRad = Nstd * np.sqrt(
    #     np.std(data[msk]['_x'])**2 + np.std(data[msk]['_y'])**2)
    pmRad = Nstd * np.sqrt(MAD(members['pmRA'])**2 + MAD(members['pmDE'])**2)
    # pmRad = Nstd * np.sqrt(
    #     np.std(members['pmRA'])**2 + np.std(members['pmDE'])**2)
    Plx_rad = Nstd * np.std(members['Plx'])
    # Plx_rad = Nstd * MAD(members['Plx'])
    # Plx_rad = abs(data[msk]['Plx'] - Plx_m).max()
    print("Nm={}, pmRad={:.2f}, Plx_rad={:.2f}".format(
        len(members), pmRad, Plx_rad))

# Save members
fout = out_dat_folder + fname
members.write(fout, format='csv', overwrite=True)

# Plot results
fig = plt.figure(figsize=(15, 15))
plt.suptitle("Estimated number of members: {}".format(len(members)))

plt.subplot(221)
xyRad = np.sqrt((members['_x'] - RA_m)**2 + (members['_y'] - DE_m)**2).max()
plt.title("rad = {:.2f}".format(xyRad))
plt.scatter(field['_x'], field['_y'], c='grey', alpha=.3, marker='.', s=2)
plt.scatter(members['_x'], members['_y'], edgecolor='C0',
            facecolor='none', zorder=4)
plt.scatter(RA_m, DE_m, marker='x', color='r', s=50, zorder=5)

plt.subplot(222)
pmRad = np.sqrt(
    (members['pmRA'] - pmRA_m)**2 + (members['pmDE'] - pmDE_m)**2).max()
plt.title("rad = {:.2f}".format(pmRad))
plt.scatter(field['pmRA'], field['pmDE'], c='grey', alpha=.3, marker='.', s=5)
plt.scatter(members['pmRA'], members['pmDE'], edgecolor='C0',
            facecolor='none', zorder=4)
plt.scatter(pmRA_m, pmDE_m, marker='x', color='r', s=50, zorder=5)
plt.xlim(-8, -2.)
plt.ylim(8, 14)

plt.subplot(223)
plt.title("rad = {:.2f}".format(Plx_rad))
msk = (field['Plx'] > 1.) & (field['Plx'] < 3.)
plt.hist(field['Plx'][msk], 50, alpha=.5)
plt.hist(members['Plx'], 25, alpha=.5, zorder=5)
plt.axvline(Plx_m, c='r', ls='--')
plt.xlim(1.95, 3.09)

plt.subplot(224)
plt.scatter(field['BP-RP'], field['Gmag'], c='grey', alpha=.3, marker='.', s=5)
plt.scatter(members['BP-RP'], members['Gmag'], edgecolor='C0',
            facecolor='none', zorder=5)
plt.plot((x1, x2), (y1, y2), 'r-', zorder=5)
plt.xlim(-.2, 3.)
plt.ylim(19.1, 4)

# plt.show()
fig.tight_layout()
fout = out_fig_folder + fname.replace('.dat', '.png')
plt.savefig(fout, dpi=150, bbox_inches='tight')


"""
Estimate the number of members per quadrant, as the radius increases.
"""

# x, y = data['_x'], data['_y']
# cent_d = np.sqrt((x - RA_m)**2 + (y - DE_m)**2)
# xy = np.array([x, y]).T
# xy_dists = cdist([(RA_m, DE_m)], xy)[0]
# fr_area = np.ptp(x) * np.ptp(y)
# fr_quad_area = fr_area / 4.
# # print("Total area:  {:.2f} [deg^2]".format(fr_area))

# # Top left quadrant
# msk1c = (x < RA_m) & (y >= DE_m)
# # Top right quadrant
# msk2c = (x >= RA_m) & (y >= DE_m)
# # Bottom left quadrant
# msk3c = (x < RA_m) & (y < DE_m)
# # Bottom right quadrant
# msk4c = (x >= RA_m) & (y < DE_m)

# print("Q1 density: {:.2f} [arcmin^2]".format(
#     msk1c.sum() / (fr_quad_area * 3600)))
# print("Q2 density: {:.2f} [arcmin^2]".format(
#     msk2c.sum() / (fr_quad_area * 3600)))
# print("Q3 density: {:.2f} [arcmin^2]".format(
#     msk3c.sum() / (fr_quad_area * 3600)))
# print("Q4 density: {:.2f} [arcmin^2]".format(
#     msk4c.sum() / (fr_quad_area * 3600)))

# Q1234_members = [[], [], [], []]
# rads = np.arange(.2, 2, .05)
# for rad in rads:
#     cl_area = np.pi * rad**2
#     cl_quad_area = cl_area / 4.
#     quad_area = fr_quad_area - cl_quad_area

#     for i, msk in enumerate((msk1c, msk2c, msk3c, msk4c)):
#         # Stars in quadrant outside 'rad'
#         Q_out_rad = (msk & (xy_dists > rad)).sum()
#         # Quadrant density
#         Q_dens = Q_out_rad / quad_area
#         # Total stars in quadrant within 'rad'
#         Q_in_rad = (msk & (xy_dists < rad)).sum()
#         # Estimated members in quadrant within 'rad'
#         Q_membs = Q_in_rad - Q_dens * cl_quad_area

#         Q1234_members[i].append(Q_membs)
#         # Q1234_members[i].append(max(0, Q_membs))

# Q1234_members = np.array(Q1234_members)

# fig = plt.figure(figsize=(10, 5))

# plt.plot(rads, Q1234_members[0], label='Q1', ls=":")
# plt.plot(rads, Q1234_members[1], label='Q2', ls=":")
# plt.plot(rads, Q1234_members[2], label='Q3', ls=":")
# plt.plot(rads, Q1234_members[3], label='Q4', ls=":")
# plt.plot(rads, Q1234_members.sum(0), label='Sum', lw=2.)
# # plt.plot(rads, Q1234_members.mean(0), label='Mean')
# plt.xlabel("Radius")
# plt.ylabel("N members")
# # plt.ylim(0, 1700)
# plt.legend()

# plt.grid()

# fig.tight_layout()
# fout = out_fig_folder + fname.replace('.dat', '_2.png')
# plt.savefig(fout, dpi=150, bbox_inches='tight')
