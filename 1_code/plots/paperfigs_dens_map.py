
import numpy as np
from scipy import stats
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

# Select the magnitude cut
mag_lim = 15

in_file = '../2_pipeline/1_data_filter/out/NGC2516.dat'
out_fig = '../2_pipeline/paper_figs/dens_map.png'
out_file = '../2_pipeline/paper_figs/NGC2516_{}.dat'.format(mag_lim)

data = Table.read(in_file, format='ascii')

def main():

    # xx = np.random.uniform(-2, 10, 2)
    # xmin, xmax = min(xx), max(xx)
    # yy = np.random.uniform(-5, 7, 2)
    # ymin, ymax = min(yy), max(yy)
    # xv = np.random.uniform(xmin, xmax, 500)
    # yv = np.random.uniform(ymin, ymax, 500)
    # values = np.array([xv, yv])
    # dataMirror(values)

    x_data, y_data = data['_x'], data['_y']

    # Use half of Scotts factor (scipy's default).
    values = np.vstack([x_data, y_data])
    kernel = stats.gaussian_kde(values)
    bdw = kernel.covariance_factor() * np.max(values.std(axis=1)) * .5
    print("Bandwidth: ", bdw)

    mag_msk = data['Gmag'] < mag_lim

    ascii.write(data[mag_msk], out_file)

    gs_unit = 5
    gs_x, gs_y = 3, 1
    fig = plt.figure(figsize=(gs_unit * gs_x, gs_unit * gs_y))
    gs = gridspec.GridSpec(gs_y, gs_x)

    def mplot(ax, x, y, symb, mag_lim, bdw=None):
        im = mapPlotKDE(x, y, symb, mag_lim, bdw)
        # im = mapPlotHist(x, y, symb, mag_lim)

        ax.set_xticks((-2, -1, 0, 1, 2))
        ax.set_yticks((-2, -1, 0, 1, 2))
        plt.gca().invert_xaxis()
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="2%", pad=0.05)
        # plt.colorbar(im, cax=cax)
        ax.set_xlabel(r"$\alpha^{*}$")
        ax.set_ylabel(r"$\delta^{*}$")

    ax = plt.subplot(gs[0:1, 0:1])
    mplot(ax, x_data, y_data, None,
          (data['Gmag'].min(), data['Gmag'].max()), bdw)

    ax = plt.subplot(gs[0:1, 1:2])
    mplot(ax, x_data[mag_msk], y_data[mag_msk], "<", mag_lim, bdw)

    ax = plt.subplot(gs[0:1, 2:3])
    mplot(ax, x_data[~mag_msk], y_data[~mag_msk], r"\geq", mag_lim, bdw)

    # ax = plt.subplot(gs[0:1, 2:3])
    # plt.hist(data['Gmag'], 50)
    # plt.xlabel("G")

    fig.tight_layout()
    plt.savefig(out_fig, dpi=150)


def mapPlotHist(x, y, symbol, mag_lim, gd=25):
    """
    """
    plt.title(r"G${}${:.1f}  |  N={}".format(symbol, mag_lim, len(x)))

    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    hist2d = np.histogram2d(x, y, gd)[0]
    im = plt.imshow(np.rot90(hist2d), cmap=plt.get_cmap('RdYlBu_r'),
                    extent=(xmin, xmax, ymin, ymax), interpolation='bilinear')

    # # Same as above but using plt.hist2d() and showing contour lines
    # hist2d, ybins, xbins, im = plt.hist2d(
    #     x, y, gd, cmap=plt.get_cmap('RdYlBu_r'))
    # plt.contour(
    #     hist2d, colors='#551a8b', linewidths=.5, levels=3,
    #     extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()])

    xd = (xmax - xmin) / gd
    yd = (ymax - ymin) / gd
    H_dens = hist2d / (xd * yd) * (1 / 3600)
    Hmin, Hmax = np.min(H_dens), np.max(H_dens)
    print(Hmin, Hmax)

    return im


def mapPlotKDE(xd, yd, symbol, mag_lim, bdw):
    kde, x, y, k_pos, ext_range = kde_center(xd, yd, bdw)

    if symbol is not None:
        plt.title(r"G${}${:.1f}  |  N={}".format(symbol, mag_lim, len(xd)))
    else:
        plt.title(r"{:.1f}$\leq$G$\leq${:.1f}  |  N={}".format(
            mag_lim[0], mag_lim[1], len(xd)))
    kde = np.reshape(k_pos.T, x.shape)
    im = plt.imshow(
        np.rot90(kde), cmap=plt.get_cmap('RdYlBu_r'), extent=ext_range)
    plt.contour(x, y, kde, colors='#551a8b', linewidths=0.5)

    return im


def kde_center(x_data, y_data, bdw, gd=50):
    """
    Find the KDE maximum value pointing to the center coordinates.
    """
    xmin, xmax = np.min(x_data), np.max(x_data)
    ymin, ymax = np.min(y_data), np.max(y_data)
    values = np.vstack([x_data, y_data])

    # Fix KDE artifact around the borders
    values = dataMirror(values)

    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(
        values, bw_method=bdw / np.max(values.std(axis=1)))

    # Grid density (number of points).
    gd_c = complex(0, gd)
    # Define x,y grid.
    x_grid, y_grid = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x_grid.ravel(), y_grid.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)

    return kernel, x_grid, y_grid, k_pos, (xmin, xmax, ymin, ymax)


def dataMirror(data, perc=.05):
    """
    Mirror a small percentage of data in all borders to remove the KDE artifact
    that lowers the density close to the edges.

    Source: https://stackoverflow.com/a/33602171/1391441
    """

    # Reshape to: (N_points, 2)
    data = data.T
    # Box
    xmin, ymin = np.min(data, 0)
    xmax, ymax = np.max(data, 0)
    midy = (ymax - ymin) * .5

    points_left_up = np.copy(data)
    points_left_up[:, 0] = xmin - (points_left_up[:, 0] - xmin)
    points_left_up[:, 1] += midy

    points_left_down = np.copy(data)
    points_left_down[:, 0] = xmin - (points_left_down[:, 0] - xmin)
    points_left_down[:, 1] -= midy

    points_right_up = np.copy(data)
    points_right_up[:, 0] = xmax + (xmax - points_right_up[:, 0])
    points_right_up[:, 1] += midy

    points_right_down = np.copy(data)
    points_right_down[:, 0] = xmax + (xmax - points_right_down[:, 0])
    points_right_down[:, 1] -= midy

    points_down = np.copy(data)
    points_down[:, 1] = ymin - (points_down[:, 1] - ymin)
    points_up = np.copy(data)
    points_up[:, 1] = ymax + (ymax - points_up[:, 1])

    # plt.scatter(*data.T, c='g', alpha=.3)
    # plt.scatter(*points_left_up.T, c='b', alpha=.2)
    # plt.scatter(*points_left_down.T, c='r', alpha=.2)
    # plt.scatter(*points_right_up.T, c='b', alpha=.2)
    # plt.scatter(*points_right_down.T, c='r', alpha=.2)
    # plt.scatter(*points_up.T, c='cyan', alpha=.2)
    # plt.scatter(*points_down.T, c='orange', alpha=.2)
    # plt.show()

    # Mirrored data
    mirror_data = np.append(
        np.append(points_left_down, points_left_up, axis=0),
        np.append(points_right_down, points_right_up, axis=0), axis=0)
    mirror_data = np.append(
        mirror_data, np.append(points_down, points_up, axis=0), axis=0)

    # Combine all data
    points = np.append(data, mirror_data, axis=0)

    # Trim mirrored frame to within a 'perc' pad of the full (x, y) range
    xr, yr = np.ptp(data.T[0]) * perc, np.ptp(data.T[1]) * perc
    xmin, xmax = xmin - xr, xmax + xr
    ymin, ymax = ymin - yr, ymax + yr
    msk = (points[:, 0] > xmin) & (points[:, 0] < xmax) &\
        (points[:, 1] > ymin) & (points[:, 1] < ymax)

    return points[msk].T


if __name__ == '__main__':
    main()
