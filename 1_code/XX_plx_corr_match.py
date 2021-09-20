
import numpy as np
from astropy.io import ascii
from scipy.optimize import differential_evolution as DE

"""
Match the file with the selected members, with the large file with the
corrected parallax. This way we can extract the selected members with their
corrected parallax values, without having to use the cross-match script.

Then, apply a simple bootstrap on the likelihood used by AsteCA to estimate
the parallax distance.
"""

data_plx_c = ascii.read("../2_pipeline/XX_plx_corr/NGC2516_plx_corr.dat")
IDs_plx_c = list(data_plx_c['EDR3Name'])
masses = ['1000', '1100', '1200', '1300', '1400', '1450', '1500', '1600',
          '1700', '1800']


def main(Nboot=100):
    """
    """
    for mass in masses:
        data_membs = ascii.read(
            "../2_pipeline/XX_ASteCA/in/NGC2516_{}.dat".format(mass))
        plx_c = []
        for st in data_membs:
            idx = IDs_plx_c.index(st['EDR3Name'])
            plx_c.append(data_plx_c['Plx'][idx])

        plx_c = np.array(plx_c)

        dmin, dmax = 1. / plx_c.max(), 1. / plx_c.min()
        bounds = [[dmin, dmax]]

        dvals = []
        for _ in range(Nboot):
            # if _ % 100 == 0:
            #     print(_)
            plx_c_r = plx_c + data_membs['e_Plx'] * np.random.normal(
                0., 1., len(plx_c))
            result = DE(lnlike, bounds, args=(plx_c_r, data_membs['e_Plx']))
            dvals.append(result.x[0])

        dvals = 1000. * np.array(dvals)
        print(mass, np.mean(dvals), np.std(dvals))


def lnlike(d, plx, e_plx):
    """
    Simple likelihood used in Cantat-Gaudin et al. (2018), 'A Gaia DR2 view of
    the Open Cluster population in the Milky Way'

    The final estimated value is almost always equivalent to the weighted
    average.
    """
    return np.sum(((plx - 1. / d) / e_plx)**2)


if __name__ == '__main__':
    main()
