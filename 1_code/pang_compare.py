
from astropy.io import ascii
import pathlib
import numpy as np
import matplotlib.pyplot as plt


"""
Compare our members selection with that from Pang
"""


data_pang = ascii.read('../2_pipeline/pang_data_format/ngc2516_pang.dat')
data = ascii.read('../2_pipeline/4_members_estim/out/NGC2516_15deg.dat')
print("Pang members: {}, our members: {}".
      format(len(data_pang), len(data)))

cent = (np.median(data_pang['Xcor']), np.median(data_pang['Ycor']), np.median(data_pang['Zcor']))
from scipy.spatial.distance import cdist
dist0=cdist(np.array([cent]), np.array([data_pang['Xcor'], data_pang['Ycor'], data_pang['Zcor']]).T)
breakpoint()
# plt.subplot(221)
# plt.scatter(data_pang['_x'], data_pang['_y'], c='r', alpha=.5,
#             marker='.')
# plt.scatter(data['_x'], data['_y'], c='g', alpha=.5, marker='.')
# plt.subplot(222)
# plt.scatter(data_pang['BP-RP'], data_pang['Gmag'], marker='.',
#             c='r', alpha=.5)
# plt.scatter(data['BP-RP'], data['Gmag'], c='g', alpha=.5, marker='.')
# plt.gca().invert_yaxis()
# plt.subplot(223)
# plt.scatter(data_pang['pmRA'], data_pang['pmDE'],
#             marker='.', c='r', alpha=.5)
# plt.scatter(data['pmRA'], data['pmDE'], c='g', alpha=.5, marker='.')
# plt.subplot(224)
# plt.hist(data_pang['Plx'], 25, color='r', alpha=.5)
# plt.hist(data['Plx'], 25, alpha=.5, color='green')
# plt.show()

# IDs in our dataset
full_IDs = [int(_['EDR3Name'].split(' ')[-1]) for _ in data]
pang_us_match, pang_no_match = [], []
for i, st in enumerate(data_pang['GaiaEDR3']):
    try:
        # Index of Pang's member in our data
        j = full_IDs.index(st)
        pang_us_match.append(j)
    except ValueError:
        # print(st, [data_pang[i][_] for _ in ('Plx', 'pmRA', 'pmDE', 'Gmag',
        # 'BP-RP')])
        # Index of Pang's member in Pang's data
        pang_no_match.append(i)
pang_us_match = data[np.array(pang_us_match)]
pang_us_no_match = data_pang[np.array(pang_no_match)]
print("Pang stars matched: {}".format(len(pang_us_match)))
print("Pang stars not matched: {}".format(len(pang_us_no_match)))

# IDs in Pang's dataset
full_IDs = [_ for _ in data_pang['GaiaEDR3']]
pang_us_match, us_pang_no_match = [], []
for i, st in enumerate(data['EDR3Name']):
    try:
        # Index of our member in Pang's data
        j = full_IDs.index(int(st.split(' ')[-1]))
        pang_us_match.append(j)
    except ValueError:
        # Index of our member in our data
        us_pang_no_match.append(i)
pang_us_match = data_pang[np.array(pang_us_match)]
us_pang_no_match = data[np.array(us_pang_no_match)]
print("Our stars matched: {}".format(len(pang_us_match)))
print("Our stars not matched: {}".format(len(us_pang_no_match)))

plt.subplot(221)
plt.scatter(pang_us_no_match['_x'], pang_us_no_match['_y'], c='r', alpha=.5,
            marker='.')
plt.scatter(us_pang_no_match['_x'], us_pang_no_match['_y'], c='g', alpha=.5,
            marker='.')
plt.subplot(222)
plt.scatter(pang_us_no_match['BP-RP'], pang_us_no_match['Gmag'], marker='.',
            c='r', alpha=.5)
plt.scatter(us_pang_no_match['BP-RP'], us_pang_no_match['Gmag'], c='g',
            alpha=.5, marker='.')
plt.gca().invert_yaxis()
plt.subplot(223)
plt.scatter(pang_us_no_match['pmRA'], pang_us_no_match['pmDE'],
            marker='.', c='r', alpha=.5)
plt.scatter(us_pang_no_match['pmRA'], us_pang_no_match['pmDE'], c='g',
            alpha=.5, marker='.')
plt.subplot(224)
plt.hist(pang_us_no_match['Plx'], 25, color='r', alpha=.5)
plt.hist(us_pang_no_match['Plx'], 25, alpha=.5, color='green')
plt.show()
