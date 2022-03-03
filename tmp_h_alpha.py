import logging
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
# from astropy.wcs import WCS
from astropy.io import fits
import json
from glob import glob
from IPython.display import display
try:
    from academicpython.json_encoder import MyEncoder, NoIndent
    # from academicpython import json_encoder as jen
    print("MyEncoder imported")
except ModuleNotFoundError:
    MyEncoder = json.JSONEncoder
    NoIndent = lambda x: x

from academicpython import plottools as pt
import imp
imp.reload(pt)
from sed_tools import combine_spec, get_sed
sys.path.append("/startrek2nb/chongchong/Sam/coding")
from pkgs.radiation_module import DustExtinction
from pkgs import ramses

logger = logging.getLogger('dev')
logger.setLevel(logging.DEBUG)

class MyLog:

    def __init__(self):
        self.count = 0
        self.dp = None

    def print(self, *args):
        if self.count == 0:
            print(*args)
        elif self.count == 1:
            print(".")
            self.dp = display(*args, display_id=True)
        else:
            self.dp.update(*args)
            # sys.stdout.write("\033[F") #back to previous line
            # sys.stdout.write("\033[K") #clear line
            # sys.stdout.write("\033[F") #back to previous line
            # sys.stdout.write("\033[K") #clear line
            # print(*args)
        self.count += 1
        return

    def print(self, *args):
        print(*args)


def read_dat(fn):
    f = FortranFile(fn, 'r')
    datasize = f.read_ints(dtype='int32')
    data = f.read_reals(dtype='float32')
    f.close()
    assert data.shape[0] == np.product(datasize), \
        f"shape: {data.shape}, size: {datasize}"
    return data.reshape(datasize)


def lum_halpha(n, V):
    """get the H-alpha luminosity of each grid
    n (array): ionized-hydrogen (electron) number density, cm^-3
    V (float, or array with same shape as n): volume, cm^3

    Reference: MBW Eq. (10.96)
    """

    hu_halpha = 5.447e-12  # in erg, = 3.4eV
    alphaB = 2.56e-13     # case B recombination coeff, in cm3 s-1
    ndot = alphaB * n * n * V
    return 0.450 * hu_halpha * ndot   # cgs


def halpha_sb(fn_den, fn_xHII, width, axis=0, with_dust=1):
    """
    fn_den (str): filename of the density data (output of amr2cube.f90)
    fn_xHII (str): filename of the xHII data (output of amr2cube.f90)
    width (float): half-width of the box (in cm)??? --> side of the box

    return: surfb_ext, units: erg s-1 cm-2 arcsec-2
    """

    den = read_dat(fn_den)
    xHII = read_dat(fn_xHII)

    #<<< H-alpha radition density
    nside = den.shape[0]
    dx = width / nside  # cm
    vol = dx**3   # cgs
    surfb_grids = np.float32(lum_halpha(den * xHII, vol)) # erg s^-1
    surfb_grids /= (4 * np.pi * dx**2)  # erg s-1 cm-2 sr-1, grids. !!! large
    del xHII
    # surfb_noext = np.sum(surfb_grids, axis=axis)

    if with_dust:
        #<<< Dust extinction
        de = DustExtinction(Z=1.0, cloud='SMC')   # solar Z, SMC
        lam_halpha = 0.65628   # um
        de_halpha = de.get_tau(lam_halpha, 1e21)   # tau per 10^21 cm^-2
        NH21 = np.cumsum(den, axis=axis) * dx * 1e-21 # !!! large
        del den
        ext_grids = np.exp(-1. * de_halpha * NH21) # !!! large
        del NH21
        surfb_ext = np.sum(surfb_grids * ext_grids, axis=axis) # erg s-1 cm-2 sr-1
        del ext_grids
    else:
        surfb_ext = np.sum(surfb_grids, axis=axis) # erg s-1 cm-2 sr-1
    surfb_ext *= 2.35044e-11                               # erg s-1 cm-2 arcsec-2

    # # save data
    # fn = f"../data/run_v5/spectrum_Job{jobid}/out{out}.txt"
    # os.makedirs(os.path.dirname(fn), exist_ok=1)
    # np.savetxt(fn, surfb_ext)

    return surfb_ext


def do_full_spec(jobid, data_base_dir, plot_base_dir, halpha_data_base_dir,
                 outs=None, feature='main_CK_1e6ph_sedsb'):
    """Combine H-alpha with continuum from dust extinction

    Example:
    >>> do_full_spec(
            jobid='2.2.2',
            data_base_dir="../data/yorp07/run_v6/Job2.2.2",
            plot_base_dir="../results/fig/run_v6/Halpha",
            halpha_data_base_dir="/startrek2nb/chongchong/Sam/test_amr2cube",
            feature='main_CK_corrected_1e6ph_sedsb',
        )

    """

    # jobid = "2.2.2"
    # pt.set_plotdir(f"../results/fig/run_v6/Halpha/Job{jobid}")
    pt.set_plotdir(f"{plot_base_dir}/Job{jobid}")
    VERSION = "yorp07/run_v6"

    # halpha_data_dir = f"/startrek2nb/chongchong/Sam/test_amr2cube/Job{jobid}"
    halpha_data_dir = f"{halpha_data_base_dir}/Job{jobid}"
    fn_json = f"{plot_base_dir}/Job{jobid}/data.json"
    lmax = 9
    axis = 1
    radius = 0.49
    r = ramses.Ramses(jobid, jobdir="/startrek2nb/chongchong/Sam")
    width = 2 * radius * r.boxlen * r.unit_l_code   # cm
    ha_width = 0.6 / 8  # micron, JWST, source at z=8

    if outs is None:
        _outs = r.get_all_outputs()
        outs = []
        for i in _outs:
            if os.path.isdir(f"{data_base_dir}/out{i}") and\
               os.path.isfile(f"{halpha_data_dir}/out{i}_den_l{lmax}.dat"):
                if len(glob(f"{data_base_dir}/out{i}/{feature}*")) > 0:
                    outs.append(i)
    if outs == []:
        raise SystemExit("outs is empty")
    # if os.path.isfile(fn_json):
    #     with open(fn_json, 'r') as fin:
    #         dic = json.load(fin)
    #     halpha_strength = dic['halpha_strength']
    #     halpha_strength_dict = dic['halpha_strength_dict']
    # else:
    dic = {
        'outs': [],
        'dust_wavelength': {},
        'dust_spec': {},
        'dust_spec_unit': 'erg / (s cm2 arcsec2 micron)',
        'halpha_strength': [],
        'halpha_strength_dict': {},
        'combined_spec': {},
        'combined_spec_unit': 'erg / (s cm2 arcsec2 micron)',
    }
    halpha_strength = []
    halpha_strength_dict = {}
    halpha_no_dust_strength = []
    halpha_no_dust_strength_dict = {}
    dic['outs'] = NoIndent(list(outs))
    # spectrum width of H-alpha (for combining with dust continuum)
    dic['halpha_width'] = ha_width
    dic['halpha_width_unit'] = "micron"
    dic['halpha_width_about'] = "Highest resolution of source at z = 8 from JWST. "\
        "0.6 micron is the highest spectrum resolution of JWST."
    mylog = MyLog()
    is_plot = 0

    #outs = [44]
    outs = [26]
    outs = [38]

    for out in outs:
        # if os.path.isfile(f"{pt.PLOT_DIR}/sed_with_ha_out{out}.pdf"):
        #     mylog.print(f"{pt.PLOT_DIR}/sed_with_ha_out{out}.pdf exists. Skipped")
        #     continue
        #-------- spatially integrated spectrum of continuum --------#
        # root_path = f"../data/{VERSION}/sam_{jobid}_out{out}"
        root_path = f"{data_base_dir}/out{out}"
        #feature = 'main_CK_1e6ph_sedsb'
        fn_fits = f"{root_path}/{feature}_total.fits"
        #wavel, spec = get_sed(f'{root_path}/{feature}')
        #wavel2, spec2 = get_sed(f'{root_path}/{feature}', kind="fits", is_int=1)

        ## plot comparison between spec and spec2
        #compare_spec12 = 0
        #if compare_spec12:
        #    with plt.style.context(['science', 'no-latex']):
        #        plt.figure()
        #        plt.plot(wavel, spec, 'C0', label="SED")
        #        plt.plot(wavel2, spec2, 'C1', label="From pixel", zorder=-1)
        #        ax = plt.gca()
        #        ax.legend()
        #        ax.set(ylabel=r'$F_{\lambda}$ (erg/s/cm2/arcsec2/micron)',
        #            xlabel=r'wavelength (micron)',
        #            yscale='log', xscale='log')
        #        set_y_decades(5, ax)
        #        plt.show()
        #    return

        ## check the convergence of spec and spec2. They have to be close.
        #assert abs(spec.max() / spec2.max() - 1.0) < 0.01
        #flag = 0
        #for i in range(len(spec)):
        #    if spec[i] is np.nan:
        #        continue
        #    if np.abs(spec[i]) < 1e-20 and np.abs(spec2[i]) < 1e-20:
        #        continue
        #    # if np.abs(spec[i] / spec2[i]) - 1.0 > 1e-1:
        #    #     flag = 1
        #    assert np.abs(spec[i] / spec2[i]) < 3 and np.abs(spec[i] / spec2[i]) > 1/3, \
        #            "spec:\n{}, \nspec2:\n{}".format(spec, spec2)
        ## assert not flag

        ## integrate over space
        #with fits.open(fn_fits) as _hdus:
        #    hd = _hdus[0].header
        #    solid_angle_per_pixel = hd['CDELT1'] * hd['CDELT2']  # arcsec2
        #    total_arcsec2 = hd['NAXIS1'] * hd['NAXIS2'] * solid_angle_per_pixel
        #    assert hd['CUNIT1'] == "arcsec" and hd['CUNIT2'] == "arcsec"
        #spec_per_arcsec2 = spec / total_arcsec2   # W / (m2 arcsec2 micron)
        #spec_per_arcsec2 *= 1000.0                # erg / (s cm2 arcsec2 micron)
        #dic['dust_wavelength'][out] = NoIndent(list(wavel))
        #dic['dust_spec'][out] = NoIndent(list(spec_per_arcsec2))

        # H-alpha
        suffix = "_w"
        den = f"{halpha_data_dir}/out{out}_den_l{lmax}{suffix}.dat"
        xHII = f"{halpha_data_dir}/out{out}_xHII_l{lmax}{suffix}.dat"
        surfb_ext = halpha_sb(den, xHII, width, axis=axis) # erg s-1 cm-2 arcsec-2
        surfb_no_dust = halpha_sb(den, xHII, width, axis=axis, with_dust=0) # erg s-1 cm-2 arcsec-2
        halpha_sb_no_dust_mean = np.mean(np.mean(surfb_no_dust, axis=-1), axis=-1)  # erg s-1 cm-2 arcsec-2
        print(halpha_sb_no_dust_mean)
        return

        # plot RGB image
        is_plot_im = 0
        if is_plot_im:
            f, ax = plt.subplots()
            extent = [-radius, radius, -radius, radius]
            im = ax.imshow(np.log10(surfb_ext)[::-1, :], extent=extent,
                           vmin=-18, vmax=-12, cmap='Reds')
            ax.set(xlabel='y (code unit)', ylabel='z (code unit)', )
            ax.set_title('H-alpha flux per sr')
            cb = f.colorbar(im, ax=ax)
            cb.set_label(r'log F (erg/s/cm2/sr)')
            pt.save(f"out{out}.pang")
            plt.show()

            # save data
            fn = f"../data/{VERSION}/spectrum_Job{jobid}/out{out}.txt"
            os.makedirs(os.path.dirname(fn), exist_ok=1)
            np.savetxt(fn, surfb_ext)

        # Combine H-alpha with dust extinction
        lam_halpha = 0.65628   # um
        # dlambda_lambda = 10 / 3e5  # * km/s / c
        # lam_halpha = 0.65628   # um
        # dlam_halpha = dlambda_lambda * lam_halpha  # um
        # erg s-1 cm-2 arcsec-2
        halpha_sb_mean = np.mean(np.mean(surfb_ext, axis=-1), axis=-1)  # erg s-1 cm-2 arcsec-2
        halpha_strength.append(halpha_sb_mean.astype(float))
        halpha_strength_dict[out] = halpha_sb_mean.astype(float)
        halpha_sb_no_dust_mean = np.mean(np.mean(surfb_no_dust, axis=-1), axis=-1)  # erg s-1 cm-2 arcsec-2
        halpha_no_dust_strength.append(halpha_sb_no_dust_mean.astype(float))
        halpha_no_dust_strength_dict[out] = halpha_sb_no_dust_mean.astype(float)
        halpha_sb_mean_arr = np.array([halpha_sb_mean]).reshape([1,1])

        spec_per_arcsec2_mat = spec_per_arcsec2[:, np.newaxis, np.newaxis]
        data_with_ha = combine_spec(
            spec_per_arcsec2_mat, wavel, halpha_sb_mean_arr,
            lam_halpha, width2=ha_width)
        dic['combined_spec'][out] = NoIndent(list(data_with_ha[:, 0, 0]))
        # fmt = "{:13s}: {:.2e}, {:.2e}, {:.2e}"
        # print(fmt.format("data_with_ha", data_with_ha.max(),
        #                  data_with_ha.min(), data_with_ha.mean()))
        is_plot = 1
        if is_plot:
            plt.close()         # show only the last figure in Jupyter Notebook
            with plt.style.context(['science', 'no-latex']):
                f, ax = plt.subplots()
                plt.plot(wavel, data_with_ha[:, 0, 0])
                ax.set(xscale='log', xlabel=r"$\lambda$ ($\mu$m)", yscale='log',
                       ylabel=r'$F_{\lambda}$ [erg/(s cm$^2$ arcsec$^2$ $\mu$m)]',
                       ylim=[6e-16, 6e-11],
                )
                pt.save(f"sed_with_ha_out{out}.pdf", isprint=0)
                mylog.print(f"sed_with_ha_out{out}.pdf saved")
    if is_plot:
        print("Showing the combined spectrum from the last output:")
        plt.show()

    dic["halpha_strength"] = NoIndent(halpha_strength)
    dic["halpha_strength_dict"] = NoIndent(halpha_strength_dict)
    dic["halpha_strength_unit"] = 'erg / (s cm2 arcsec2)'
    dic["halpha_no_dust_strength"] = NoIndent(halpha_no_dust_strength)
    dic["halpha_no_dust_strength_dict"] = NoIndent(halpha_no_dust_strength_dict)
    dic["halpha_no_dust_strength_unit"] = 'erg / (s cm2 arcsec2)'
    jdump = json.dumps(dic, indent=2, cls=MyEncoder, sort_keys=True)
    with open(fn_json, 'w') as fout:
        fout.write(jdump)
    # with open(fn_json, 'w') as fout:
    #     json.dump(dic, fout, indent=2, cls=MyEncoder, sort_keys=True)
    print(fn_json, 'saved!')
    return

if __name__ == "__main__":
    jobid = '3.2.2'
    data_base_dir = f"../data/yorp07/run_v6/Job{jobid}"
    plot_base_dir = "../results/fig/run_v6/Halpha"
    do_full_spec(
        jobid=jobid,
        data_base_dir=data_base_dir,
        plot_base_dir=plot_base_dir,
        halpha_data_base_dir="../data/data_amr2cube",
        outs=None,
        feature='main_w_1e6ph_sedsb',
    )

