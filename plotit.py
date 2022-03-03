import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import astropy.units as u, astropy.constants as c
from astropy.wcs import WCS
from astropy.io import fits
import json
import warnings

from CIE_color_system import CIE_color_fast as CIE

#matplotlib.rcParams['figure.dpi'] = 300

# micron = u.def_unit('micron', u.um)
# u.add_enabled_units([micron])


def rect(y, x):
    """Integrate y over x.

    y can be multidimension arrays but only the first dimension is integrated.

    """

    dx = x[1:] - x[:-1]
    yleft = y[:-1, ...]
    ret = np.dot(np.transpose(yleft), dx)  # np.dot sums over the last
                                           # dimension
    return np.transpose(ret)


def abmag(fl, l):
    """Return the AB magnitude given the spectrum intensity and wavelengths.

    fl: array, f_lambda (W/m2/micron)
    l: array, wavelengths (micron)

    Ref: https://en.wikipedia.org/wiki/AB_magnitude

    """

    light_speed = c.c.to('um/s').value  # in micron
    jy = 1e-26                        # in W/m2/Hz
    v = light_speed / l               # in Hz

    #fv = l**2 / light_speed * fl      # in W/m2/Hz
    fv = np.zeros(fl.shape)
    fv_over_v = np.zeros(fl.shape)
    for i in range(len(l)):
        fv[i, ...] = fl[i, ...] * l[i]**2
        fv_over_v[i, ...] = fv[i, ...] / v[i]
    fv /= light_speed

    #nume = trapz(fv / v, v)           # in W/m2/Hz
    #denom = trapz(1. / v, v)          # in 1
    nume = rect(fv_over_v, v)
    denom = trapz(1. / v, v)

    #l_piv = np.sqrt(trapz(np.ones(len(l)), l) / trapz(l**-2, l))  # pivot wavelengths
    #print('l_piv =', l_piv)
    denom *= 3631 * jy
    mAB = -2.5 * np.log10(nume / denom)
    return mAB


def color_filter(data, wavelengths, scaledown=0.1, enhancement=3.0):
    """A color filter that converts spectrums into RGB.

    Args:
        data (array of shape (n, m1, m2)): The spectrum data
        wavelengths (array of shape (n, )): wavelengths in micron
        scaledown (float, default=0.1): set Y=1 at scaledown * Y.max()
        enhancement (float, default=3.0): enhancement factor for the dim pixels

    Returns:
        RGB (array of shape (m1, m2, 3))
    """

    n_lam, m1, m2 = data.shape
    assert n_lam == len(wavelengths)
    lams_nm = 1000. * wavelengths
    spec = data.T.reshape([m2 * m1, n_lam])
    XYZ = CIE.XYZ_from_spec(lams_nm, spec) # (m2 * m1, 3)
    Y = XYZ[:, 1]
    Ymax = scaledown * Y.max()
    scale = 10**(np.log10(Y/Ymax) / enhancement) * 1 / Y
    # scale[scale > 1] = 1.0
    for i in range(3):
        XYZ[:, i] *= scale
    RGB = CIE.XYZ_to_sRGB(XYZ.T, is_norm=0)  # (3, m2 * m1)
    pick = np.any(RGB > 1, axis=0)
    RGB[:, pick] /= np.max(RGB[:, pick], axis=0)
    RGB[RGB < 0] = 0
    RGB = RGB.reshape([3, m2, m1]).T
    return RGB


def combine_spec(specs, wavelengths, spec2, wavelength2, width2=None):
    """Combine pixels of spectrum with pixels of single-line surface brightness.

    Note that the units of spec2 / wavelength2 should be the same as that of
    specs.

    Parameters
    ----------
    specs: (m, n1, n2, ...) array
        The pixeled spectrum 
    wavelengths: (m, ) array
        List of wavelengths
    spec2: (n1, n2, ...) array 
        The single-line surface brightness (e.g. Ly-alpha, H-alpha)
    wavelength2: double
        The wavelength of spec2
    width2: double (Default: None)
        The width of spec2. If set to None (default), use the resolution at the
        corresponding wavelength in wavelengths.

    """

    assert specs.shape[1:] == spec2.shape
    pick = np.argmax(wavelengths >= wavelength2)
    dlam_instru = wavelengths[pick + 1] - wavelengths[pick]  # micron
    if width2 is None:
        width2 = dlam_instru
    # elif width2 > dlam_instru:
    #     warnings.warn("Your width2 is larger than the width of the instrument. "
    #                   "This is not supported right now. Using the instrumental width.")
    #     width2 = dlam_instru
    else:
        raise SystemExit("width2 is not supported yet.")
    specs_cmb = specs.copy()
    specs_cmb[pick, ...] += spec2 / width2
    return specs_cmb


def main(name, lambd=-1):

    fn = name + '.fits'
    fn_mod = name + '_mod.fits'
    if not os.path.isfile(fn_mod):
        with fits.open(fn) as hdul:
             hdul[0].header['CTYPE3'] = ''
             hdul[0].header['CUNIT3'] = 'um'
             hdul.writeto(fn_mod)

    hdu = fits.open(fn_mod)[0]
    data = hdu.data
    hdu1 = fits.open(fn_mod)[1]
    wcs = WCS(hdu.header)
    print("hdu header:")
    print(repr(hdu.header))

    wavelengths = np.array([wavel[0] for wavel in hdu1.data])
    if lambd > 0:
        # at a given wavelength
        wavel = lambd
        idx = np.argmax(wavelengths >= wavel)
        plot_data = data[idx]
        unit_to_mag = 3631e-26
        mAB = -2.5 * np.log10(plot_data / unit_to_mag)
    else:
        # bandpass mAB
        mAB = abmag(data, wavelengths)
        #mAB = np.zeros((data.shape[1], data.shape[2]))
        #for i in range(mAB.shape[0]):
        #    for j in range(mAB.shape[1]):
        #        mAB[i, j] = abmag(data[:, i, j], wavelengths)

    print('mAB shape:', mAB.shape)

    #w = 120
    #mAB = mAB[(125-w//2):(125+w//2), (125-w//2):(125+w//2)].transpose()

    # vmax = mAB[mAB < 1e20].max()
    # vmin = mAB.min()
    # print('vmin, vmax =', vmin, vmax)

    #with plt.style.context(['science']):
    fig = plt.figure()
    ax = plt.subplot(projection=wcs, slices=('x', 'y', 200))
    cmap = plt.cm.gray_r
    # cmap = plt.cm.ScalarMappable(cmap=plt.cm.gray)
    # cmap.set_clim(vmin, vmax)
    cmap.set_bad('k', 1.0)
    im = ax.imshow(mAB, cmap=cmap) #, vmax=vmax, vmin=vmin)
    cb = plt.colorbar(im, extend='max')
    cb.set_label('magnitude/arcsec$^2$')
    # ax.set_title("{:.0f} Angstrom".format(10000 * hdu1.data[idx][0]))
    if lambd > 0:
        ax.set_title("{} micron".format(lambd))
    else:
        ax.set_title("Bandpass AB magnitude")
    ax.autoscale(tight=True)
    fig.savefig(name + '.png', dpi=300)


if __name__ == "__main__":
    main(sys.argv[1])
