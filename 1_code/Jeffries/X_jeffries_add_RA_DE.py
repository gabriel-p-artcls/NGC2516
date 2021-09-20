
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import Angle

"""
Convert the (RA, DE) coordinates in the original data file (in hms and dms)
to degrees.
"""

in_folder = '../0_data/'
fname = 'jeffries_WEBDA.dat'
out_folder = '../2_pipeline/2_jeffries_add_RA_DE/out/'

data = ascii.read(in_folder + fname)

ra = Angle((data['col2'], data['col3'], data['col4']), unit=u.hourangle)
de = Angle((data['col5'], data['col6'], data['col7']), unit=u.degree)

# Remove columns: RAh RAm RAs  DE- DEd DEm DEs
data.remove_columns(['col2', 'col3', 'col4', 'col5', 'col6', 'col7'])
data.add_column(ra.deg, name='RA', index=1)
data.add_column(de.deg, name='DE', index=2)

ascii.write(
    data, out_folder + fname, format="fixed_width", delimiter=' ',
    names=('ID', 'RAJ2000', 'DEJ2000', 'Vmag', 'B-V', 'V-I', 'F1', 'F2',
           'F3', 'F4', 'F5'), overwrite=True)
