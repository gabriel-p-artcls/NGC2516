
import matplotlib.pyplot as plt
from astropy.io import ascii

"""
Make plots comparing the members from Jeffries (column
'F4 Candidate membership'), to those selected from the pyUPMASK analysis.
"""
pyU_method = 'GMM_3'  # 'VOR'

file_pyUPMASK = '../2_pipeline/7_members_filter/out/NGC2516_{}.dat'.format(
    pyU_method)
file_jeffries = '../2_pipeline/4_extract_members/out/jeffries_EDR3_membs.dat'
out_folder = '../2_pipeline/8_jeffries_vs_pyUPMASK/out/'

jeff_membs = ascii.read(file_jeffries)
pyUPMASK_membs = ascii.read(file_pyUPMASK)

fig = plt.figure(figsize=(20, 10))
plt.subplot(121)
plt.scatter(jeff_membs['RA_ICRS'], jeff_membs['DE_ICRS'], c='r', marker='x',
            s=15, lw=.8, label="Jeffries")
plt.scatter(
    pyUPMASK_membs['RA_ICRS'], pyUPMASK_membs['DE_ICRS'], facecolor='none',
    edgecolor='g', marker='o', lw=.5, label="pyUPMASK")
plt.xlabel("RA")
plt.ylabel("DE")
plt.legend()
plt.gca().invert_xaxis()

plt.subplot(122)
plt.title("Zoom")
plt.scatter(jeff_membs['RA_ICRS'], jeff_membs['DE_ICRS'], c='r', marker='x',
            s=15, lw=.8)
plt.scatter(
    pyUPMASK_membs['RA_ICRS'], pyUPMASK_membs['DE_ICRS'], facecolor='none',
    edgecolor='g', marker='o', lw=.5)
plt.xlabel("RA")
plt.ylabel("DE")
plt.xlim(121, 118)
plt.ylim(-61.4, -60.1)
fig.tight_layout()
plt.savefig(out_folder + "coords.png", dpi=150, bbox_inches='tight')

fig = plt.figure(figsize=(10, 10))
plt.scatter(
    jeff_membs['BP-RP'], jeff_membs['Gmag'], c='r', marker='x', s=15, lw=.8,
    label="Jeffries, N={}".format(len(jeff_membs)))
plt.scatter(pyUPMASK_membs['BP-RP'], pyUPMASK_membs['Gmag'], facecolor='none',
            edgecolor='g', marker='o', lw=.5, label="pyUPMASK, N={}".format(
                len(pyUPMASK_membs)))
plt.xlabel("BP-RP")
plt.ylabel("G")
plt.legend()
plt.gca().invert_yaxis()
plt.savefig(out_folder + "cmd.png", dpi=150, bbox_inches='tight')

fig = plt.figure(figsize=(20, 10))
plt.subplot(121)
plt.scatter(jeff_membs['pmRA'], jeff_membs['pmDE'], c='r', marker='x', s=15,
            label="Jeffries")
plt.scatter(
    pyUPMASK_membs['pmRA'], pyUPMASK_membs['pmDE'], facecolor='none',
    edgecolor='g', marker='o', lw=.5, label="pyUPMASK")
plt.xlabel("pmRA")
plt.ylabel("pmDE")
plt.gca().invert_xaxis()
plt.legend()

plt.subplot(122)
plt.title("Zoom")
plt.scatter(jeff_membs['pmRA'], jeff_membs['pmDE'], c='r', marker='x', s=15,
            label="Jeffries")
plt.scatter(
    pyUPMASK_membs['pmRA'], pyUPMASK_membs['pmDE'], facecolor='none',
    edgecolor='g', marker='o', lw=.5, label="pyUPMASK")
plt.xlabel("pmRA")
plt.ylabel("pmDE")
plt.xlim(-2.5, -6.5)
plt.ylim(9.5, 13)
plt.savefig(out_folder + "PMs.png", dpi=150, bbox_inches='tight')

fig = plt.figure(figsize=(20, 10))
plt.subplot(121)
plt.hist(jeff_membs['Plx'], 100, color='r', histtype='step', label="Jeffries")
plt.hist(
    pyUPMASK_membs['Plx'], 20, color='g', histtype='step', label="pyUPMASK")
plt.xlabel('Plx')
plt.legend()

plt.subplot(122)
plt.title("Zoom")
plt.hist(jeff_membs['Plx'], 150, color='r', histtype='step', label="Jeffries")
plt.hist(
    pyUPMASK_membs['Plx'], 20, color='g', histtype='step', label="pyUPMASK")
plt.xlabel('Plx')
plt.xlim(2., 3.)

plt.savefig(out_folder + "Plx.png", dpi=150, bbox_inches='tight')
