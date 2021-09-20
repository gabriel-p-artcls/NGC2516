
from astropy.io import ascii
import numpy as np


"""
Prepare Pang data to be processed by ASteCA. Add 'BP-RP' color, and missing
uncertainties for G and BP-RP
"""

data_pang = ascii.read('../0_data/pang_data.dat', delimiter=' ')

# Search for the error values for each star in the Pang dataset, in our large
# dataset. Not all of the Pang stars will be present; for these stars
# we use the value of the star in the large dataset with the closets Gmag.
data_large = ascii.read('../2_pipeline/vizier_data_format/NGC2516_10deg.dat')
full_IDs = [int(_['EDR3Name'].split(' ')[-1]) for _ in data_large]
e_G, e_BPRP = [], []
for i, st in enumerate(data_pang['GaiaEDR3']):
    print(i)
    try:
        j = full_IDs.index(st)
    except ValueError:
        # Store the value attached to the star with the closest Gmag
        j = np.argmin(abs(data_large['Gmag'] - data_pang['Gmag'][i]))
    e_G.append(data_large['e_Gmag'][j])
    e_BPRP.append(data_large['e_BP-RP'][j])

# Add columns to table
data_pang.add_column(e_G, name='e_Gmag')
data_pang.add_column(e_BPRP, name='e_BP-RP')
BPRP = data_pang['BPmag'] - data_pang['RPmag']
data_pang.add_column(BPRP, name='BP-RP')

# Apply magnitude cut
msk = data_pang['Gmag'] <= 19
data_pang[msk].write('../2_pipeline/pang_data_format/ngc2516_pang.dat',
                     format='csv', overwrite=True)
