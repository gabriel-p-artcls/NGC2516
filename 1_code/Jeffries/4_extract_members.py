from astropy.io import ascii
import numpy as np

"""
From the Jeffries file cross-matched with Gaia EDR3, generate a file with
members only. I.e.: those stars with a '1' in the F4 column
"""

in_folder = '../2_pipeline/3_jeffries_EDR3_cross_match/out/'
fname = 'jeffries_EDR3.dat'
out_folder = '../2_pipeline/4_extract_members/out/'

data = ascii.read(in_folder + fname)
# Add 'BP-RP' uncertainty
data['e_BP-RP'] = np.sqrt(data['e_BPmag']**2 + data['e_RPmag']**2)

cols_keep = (
    'ID', 'RAJ2000', 'DEJ2000', 'V', 'B-V', 'V-I', 'F1', 'F2', 'F3', 'F4',
    'F5', '_x', '_y', 'RA_ICRS', 'DE_ICRS', 'Source', 'Plx', 'e_Plx', 'PM',
    'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'BP-RP', 'e_BP-RP')
col_rmv = []
for col in data.columns:
    if col not in cols_keep:
        col_rmv.append(col)
print(col_rmv)
data.remove_columns(col_rmv)

msk = data['F4'] == 1
ascii.write(data[msk], out_folder + fname.replace('.dat', '_membs.dat'),
            delimiter=',')
