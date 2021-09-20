
from astropy.io import ascii
from astropy.table import Table

"""
Remove not used columns from a data file. Makes the file smaller and more
manageable.
"""


# file = "NGC2516"
# cols = (
#     'EDR3Name', '_x', '_y', 'RA_ICRS', 'DE_ICRS', 'Plx', 'e_Plx', 'pmRA',
#     'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'BP-RP', 'e_BP-RP')

cols = (
    'EDR3Name', '_x', '_y', 'RA_ICRS', 'DE_ICRS', 'Plx', 'e_Plx', 'pmRA',
    'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'BP-RP', 'e_BPmag-RPmag')
file = "NGC2516_7deg"

in_folder = '../0_data/'
out_folder = '../2_pipeline/1_data_filter/out/'

data = Table.read(in_folder + file + '.dat', format='ascii')

data = data[cols]
print("Total number of stars in file", len(data))

# # 2x2 deg square
# msk = (data['_x'] < 2.) & (data['_x'] > -2.) & (data['_y'] < 2.) &\
#     (data['_y'] > -2.)
# data = data[msk]
# print(len(data))

# # Mag filter
# msk = data['Gmag'] < 16
# data = data[msk]
# print(len(data))
# out_file = "NGC2516_16"

ascii.write(data, out_folder + file + '.dat', format='csv',
            overwrite=True)
