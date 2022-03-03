import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u, astropy.constants as c
from astropy.wcs import WCS
from astropy.io import fits
import json

#matplotlib.rcParams['figure.dpi'] = 300

# micron = u.def_unit('micron', u.um)
# u.add_enabled_units([micron])

def main(name):

    fn = name + '.fits'
    fn_mod = name + '_mod.fits'
    if not os.path.isfile(fn_mod):
        with fits.open(fn) as hdul:
             hdul[0].header['CTYPE3'] = ''
             hdul[0].header['CUNIT3'] = 'um'
             hdul.writeto(fn_mod)

    hdu = fits.open(fn_mod)[0]
    hdu1 = fits.open(fn_mod)[1]
    wcs = WCS(hdu.header)
    print("hdu header:")
    print(repr(hdu.header))

    # find 1500 A
    wavel = 0.15
    wavelengths = np.array([wavel[0] for wavel in hdu1.data])
    idx = np.argmax(wavelengths >= wavel)

    # unit_to_mag = 2.518e-8
    unit_to_mag = 3631e-26
    #w = 120
    #plot_data = hdu.data[idx, (125-w//2):(125+w//2), (125-w//2):(125+w//2)].transpose()
    # plot_data = hdu.data[idx]
    plot_data = np.mean(hdu.data, axis=0)
    plot_data_mag = -2.5 * np.log10(plot_data / unit_to_mag)

    print('plot_data shape:', plot_data.shape)

    # vmax = plot_data_mag[plot_data_mag < 1e20].max()
    # vmin = plot_data_mag.min()
    # print('vmin, vmax =', vmin, vmax)

    #with plt.style.context(['science']):
    fig = plt.figure()
    ax = plt.subplot(projection=wcs, slices=('x', 'y', 200))
    cmap = plt.cm.gray_r
    cmap = plt.cm.viridis
    # cmap = plt.cm.ScalarMappable(cmap=plt.cm.gray)
    # cmap.set_clim(vmin, vmax)
    cmap.set_bad('k', 1.0)
    im = ax.imshow(plot_data_mag, cmap=cmap) #, vmax=vmax, vmin=vmin)
    cb = plt.colorbar(im, extend='max')
    cb.set_label('Magnitude/arcsec$^2$/micron')
    ax.set_title("{:.0f} Angstrom".format(10000 * hdu1.data[idx][0]))
    ax.autoscale(tight=True)
    fig.savefig(name + '.png', dpi=300)

if __name__ == "__main__":
    main(sys.argv[1])
