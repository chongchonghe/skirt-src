import os
import numpy as np
from astropy.io import fits
from sed_tools import combine_spec, get_sed
from ramtools import ramses
from h_alpha import halpha_sb

DATA = "../data/new_2022"

def combine_skirt_with_ha(fn_ha, fn_dust):
    """ combine SKIRT dust continuum with H-alpha """

    lam_halpha = 0.65628   # um
    with fits.open(fn_dust) as hdul:
        hdr = hdul[0].header
        data = hdul[1].data
        wavel2, spec2 = get_sed(f'{root_path}/{feature}', kind="fits", is_int=1)
        data_with_ha = combine_spec(
            spec_per_arcsec2_mat, wavel, halpha_sb_mean_arr,
            lam_halpha, width2=ha_width)

        # Todo: figure out the units of the return of get_sed(). Is it the same as the
        # original data?

def combine_fits_with_ha(hdul, ha):
    """ hdul: fits data of the stellar continuum. The dimensions of the data are (wavelength, x,
    y)

    Args:
        hdul (fits data): continuum data, should have the unit W/m2/micron/arcsec2
        ha (numpy array): h-alpha, should have the unit W/m2/arcsec2

    Return:
        (numpy array) the new combined data
    """

    assert hdul[0].header['CUNIT3'] == 'micron'
    wavelengths = np.array([wavel[0] for wavel in hdul[1].data])  # micron
    data = hdul[0].data     # W/m2/micron/arcsec2
    lam_halpha = 0.65628    # micron
    spec2 = ha               # W/m2/arcsec2
    return combine_spec(data, wavelengths, spec2, lam_halpha, width2='resol')

def main():

    # 'crab'
    center = (0.32061501, 0.49432314, 0.45865159)
    box_fraction = 0.1008558850    # box width in boxlen unit, crab
    jobid = "2.2.2.v2"
    axis = 1

    # 'crab-small'
    center = (0.32061501, 0.49432314, 0.45865159)
    box_fraction = 0.063035        # box width in boxlen unit, crab_small
    jobid = "2.2.2.v2"
    axis = 1

    # skirt = "sam19_max_1e9ph_sedsb_total.fits"  # crab
    skirt = "sam19_max_small_1e9ph_total_total.fits"  # crab-small

    # calculate H-alpha
    r = ramses.Ramses(jobid=jobid, ram_dir="/startrek/chongchong/Sam")  #TODO: change it to read_param_from_nml
    width = box_fraction * r.unit_l   # cm


    halpha_fns = "/startrek2nb/chongchong/Sam/test_amr2cube/max_{f}_l13.dat"    # crab_small
    den = halpha_fns.format(f='den')
    xHII = halpha_fns.format(f='xHII')
    surfb = halpha_sb(den, xHII, width, axis=axis)  # erg s-1 cm-2 arcsec-2
    # data_ha = np.loadtxt(fn_ha)     # erg s-1 cm-2 arcsec-2
    sb_SI_to_cgs = 1000.            # from W/m2/arcsec2 to erg/s/cm2/arcsec2
    # data_ha_SI = data_ha / sb_SI_to_cgs      # W/m2/arcsec2
    data_ha_SI = surfb / sb_SI_to_cgs      # W/m2/arcsec2

    # fn_ha = f"{DATA}/{thedir}/ha/{ha}"
    diro = "skirt+ha"
    os.makedirs(diro, exist_ok=True)
    fn_dust_and_ha = f"{diro}/{skirt.replace('.fits', '+ha.fits')}"
    fn_dust = f"skirt/{skirt}"
    with fits.open(fn_dust) as hdul:
        assert hdul[0].header['BUNIT'] == 'W/m2/micron/arcsec2', \
                f"BUNIT = {hdul[0].header['BUNIT']}"
        hdul[0].data = combine_fits_with_ha(hdul, data_ha_SI)
        hdul.writeto(fn_dust_and_ha, overwrite=True)
    return

if __name__ == "__main__":

    main()
