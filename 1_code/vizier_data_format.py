
from astropy.io import ascii
import numpy as np
from astropy.coordinates import SkyCoord

"""
Download a large region (15x15 deg) from Vizier using the TAPVizieR service
using this command:

-- output format : csv
SELECT "I/350/gaiaedr3".RA_ICRS,  "I/350/gaiaedr3".DE_ICRS,  "I/350/gaiaedr3".Plx,  "I/350/gaiaedr3".e_Plx,  "I/350/gaiaedr3".PM,  "I/350/gaiaedr3".pmRA,  "I/350/gaiaedr3".e_pmRA,
"I/350/gaiaedr3".pmDE,  "I/350/gaiaedr3".e_pmDE,  "I/350/gaiaedr3".Gmag,  "I/350/gaiaedr3".e_Gmag,  "I/350/gaiaedr3".e_BPmag,  "I/350/gaiaedr3".e_RPmag,  "I/350/gaiaedr3"."BP-RP",  "I/350/gaiaedr3".EDR3Name
FROM "I/350/gaiaedr3"
WHERE 1=CONTAINS(POINT('ICRS',"I/350/gaiaedr3".RA_ICRS,"I/350/gaiaedr3".DE_ICRS), BOX('ICRS', 119.527, -60.800, 15., 15.))
AND "I/350/gaiaedr3".Gmag<19
AND "I/350/gaiaedr3".Plx>2

(notice Gmag<19 and Plx>2)

I have to add the missing columns '_x, _y, e_BP-RP'. The BP-RP error is
estimated simply using the quadrature.
"""


def main():
    """
    """
    data = ascii.read('../0_data/NGC2516_TAPVizier.dat')
    print("Full data read")

    # Add _x, _y columns
    ra, de = data['RA_ICRS'], data['DE_ICRS']
    # Center of the frame
    ra_0, de_0 = .5 * (np.max(ra) + np.min(ra)), .5 * (np.max(de) + np.min(de))
    coord1 = SkyCoord(ra_0, de_0, frame='icrs', unit='deg')
    coord2 = SkyCoord(ra, de, frame='icrs', unit='deg')
    r = coord1.separation(coord2)
    theta = coord1.position_angle(coord2)
    x = r * np.cos(theta)
    y = r * np.sin(theta)

    data['_x'], data['_y'] = x, y

    data['e_BP-RP'] = np.sqrt(data['e_BPmag']**2 + data['e_RPmag']**2)

    data.write("../2_pipeline/vizier_data_format/NGC2516_15deg.dat",
               format='csv')
    print("Finished")


if __name__ == '__main__':
    main()
