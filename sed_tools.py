import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

# matplotlib.use("Agg")

logger = logging.getLogger('dev')
logger.setLevel(logging.DEBUG)

def get_sed(prefix, kind="sed", is_int=0):
    """Read SED from fits or .dat file
    Given a fits data with grids of spectra, return a spatially integrated spectrum

    Args
        prefix (str): prefix to get the data filename prefix_total.fits or prefix_sed.dat
        kind (str): default "sed". Either "sed" or "fits"
        is_int (bool): default 0. Effective when kind = "fits". Integrate over solid
            angle and output in the unit W/m2/um. Other wise output in unit W/m2/um/arcsec2.

    """

    assert kind in ["sed", "fits"]
    if kind == "fits":
        fn = f'{prefix}_total.fits'
        with fits.open(fn) as _hdus:
            _wavelengths = np.array([wavel[0] for wavel in _hdus[1].data])  # micron
            _data = _hdus[0].data
            _spectrum = np.mean(np.mean(_data, axis=-1), axis=-1)
            if is_int:
                hd = _hdus[0].header
                solid_angle_per_pixel = hd['CDELT1'] * hd['CDELT2']  # arcsec2
                logger.info("Units:" + str(hd['BUNIT']+'*'+hd['CUNIT1']+'*'+hd['CUNIT2']))
                _spectrum *= solid_angle_per_pixel * (_data.shape[1] * _data.shape[2])
            else:
                logger.info("Units:" + str(_hdus[0].header['BUNIT']))
    elif kind == "sed":
        fn = f'{prefix}_sed.dat'
        _wavelengths, _spectrum = np.loadtxt(fn, unpack = True)
        with open(fn, 'r') as _f:
            _f.readline()
            logger.info("Units:" + _f.readline().split(' ')[-1])
    return _wavelengths, _spectrum


def combine_spec(specs, wavelengths, spec2, wavelength2, width2='resol'):
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
        The single-line surface brightness (e.g. Ly-alpha, H-alpha). Should have the same units
        as specs EXPECT lacking a wavelength component (W/m2/arcsec2)
    wavelength2 (double)
        The wavelength of spec2 (same units as wavelengths)
    width2 ('resol' or double):
        (Default: None) The width of spec2. If set to None (default), use the
        resolution at the corresponding wavelength in wavelengths.

    """

    assert specs.shape[1:] == spec2.shape, "specs: {}, spec2: {}".format(specs.shape, spec2.shape)
    pick = np.argmax(wavelengths >= wavelength2)
    dlam_instru = wavelengths[pick + 1] - wavelengths[pick]  # micron
    ret = specs.copy()
    if width2 is "resol":
        logger.warn('Warning: Using resol of the wavelengths. This is not physical unless '\
                       'wavelengths represents real instrument.')
        ret[pick, ...] += spec2 / dlam_instru
    else:
        assert type(width2) is float
        # if width2 > dlam_instru:
        #     raise SystemExit("Your width2 is larger than the width of the instrument."\
        #                      "This is not supported right now.")
        lam_l, lam_r = wavelength2 - width2/2, wavelength2 + width2/2
        pick_l = np.argmax(wavelengths >= lam_l)
        pick_r = np.argmax(wavelengths >= lam_r)
        assert pick_l and pick_r  # either of them should not be 0
        if pick_r == pick_l:
            pick_r += 1
        ret[pick_l:pick_r, ...] += spec2 / width2
    return ret
