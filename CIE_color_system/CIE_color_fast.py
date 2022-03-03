# -*- coding: utf-8 -*-
"""
CIE 1931 color space utilities. All written in Python Vectorization. Very fast!
Author: Chong-Chong He
Date: 2020-08-09

"""
from math import exp
import numpy as np


def gaussian(x, alpha, mu, sigma1, sigma2):
    """ x is an array. Output: an array with same length. """

    def gaussian_single(_x):
        sr = (_x - mu) / (sigma1 if _x < mu else sigma2)
        return alpha * exp(-(sr * sr) / 2)
    return np.array([gaussian_single(xi) for xi in x])


def xyzfunc_from_wavelength(wavelength):
    """Input: wavelength, (n,) array of wavelengths (in nm).
    Return: (3, n) array representing x(l), y(l), z(l) """

    l_Angstrom = wavelength * 10.
    x = gaussian(l_Angstrom,  1.056, 5998, 379, 310) + \
        gaussian(l_Angstrom,  0.362, 4420, 160, 267) + \
        gaussian(l_Angstrom, -0.065, 5011, 204, 262)
    y = gaussian(l_Angstrom,  0.821, 5688, 469, 405) + \
        gaussian(l_Angstrom,  0.286, 5309, 163, 311)
    z = gaussian(l_Angstrom,  1.217, 4370, 118, 360) + \
        gaussian(l_Angstrom,  0.681, 4590, 260, 138)
    return np.stack((x, y, z))


def rect(y, x):
    """ Integrate y over x.

    x: (n,) array
    y: (n,) or (m, n) array 
        y is integrated over the last dimension. 
    
    Return: (m, 1) array
    """

    dx = x[1:] - x[:-1]
    yleft = y[..., :-1]
    return np.matmul(yleft, dx)  # np.dot sums over the last dimension


def XYZ_from_spec(lam_nm, spec):
    """ 
    lam_nm: (l,) array
        wavelength (in nm) 
    spec: (n, l) array
    
    Return: (n, 3) array

    """

    # cut wavelength in 380-780
    bigger = np.where(lam_nm >= 380)[0]
    pick = bigger[np.where(lam_nm[bigger] <= 780)]  # (n_pick,)
    # m1 = np.argmax(lam_nm >= 380)
    # m2 = np.argmax(lam_nm > 780)
    # print("m1 m2:", m1, m2)
    # assert m2 - m1 > 10
    wavelength = lam_nm[pick]                    # (m2 - m1, )
    spec_pick = spec[:, pick]                    # (n, m2 - m1)
    
    xyzfunc = xyzfunc_from_wavelength(wavelength)  # (3, m2 - m1)
    XYZ = np.zeros([spec.shape[0], 3])
    for i in range(3):
        XYZ[:, i] = rect(spec_pick * xyzfunc[i], wavelength)
    # return np.sum(xyzfunc * spec, axis=1)
    return XYZ


def XYZ_to_RGB(XYZ, is_shift=True, is_norm=True):
    """ Convert XYZ to RGB, with rescaling of XYZ and possible
    shifting (to reduce negatives) and rgb normalization.

    Ref: https://scipython.com/blog/converting-a-spectrum-to-a-colour/
    xyz: (3, n) array
    
    Return: (3, n) array
    """
    
    xyz = XYZ / np.sum(XYZ, axis=0)  # chromaticity coordinates x, y, z
    XYZ2RGB = np.array([[2.3706743, -0.9000405, -0.4706338],
                        [-0.5138850,  1.4253036,  0.0885814],
                        [0.0052982, -0.0146949,  1.0093968]])
    rgb = np.matmul(XYZ2RGB, xyz)
    if not is_shift and not is_norm:
        return rgb
    if is_shift:
        nega = np.any(rgb < 0, axis=0)  # (n,)
        rgb[:, nega] = rgb[:, nega] - np.min(rgb[:, nega], axis=0)
    if is_norm:
        rgb /= np.max(rgb, axis=0)
    return rgb


def spec_to_RGB(wavelength, spec, is_shift=True, is_norm=True):
    """ 
    wavelength: (l,) array 
    spec: (n, l) array
    
    Return: (3, n) array
    """

    XYZ = XYZ_from_spec(wavelength, spec)  # (n, 3)
    return XYZ_to_RGB(XYZ.T, is_shift=is_shift, is_norm=is_norm)


def XYZ_to_sRGB(XYZ, is_norm=True):
    """ 
    xyz: (3, n) array. Y is the luminance in (0, 1).
    
    Output: (3, n) array

    """
    
    if is_norm:
        xyz = XYZ / np.sum(XYZ, axis=0)  # chromaticity coordinates x, y, z
    else:
        xyz = XYZ
    XYZ2sRGB = np.array([[3.2404542, -1.5371385, -0.4985314],
                         [-0.9692660,  1.8760108,  0.0415560],
                         [0.0556434, -0.2040259,  1.0572252]])
    sRGB = np.matmul(XYZ2sRGB, xyz)

    # adjust function
    adj = np.zeros(sRGB.shape)
    adj = (211 * sRGB**(5/12) - 11) / 200
    pick = sRGB < 0.0031308
    adj[pick] = 323 / 25 * sRGB[pick]
    return adj


def spec_to_sRGB(wavelength, spec):
    """ Do not use! Need investigation """

    print("Do not use this function. Need investigation.")
    return

    XYZ = XYZ_from_spec(wavelength, spec)
    rgb = XYZ_to_sRGB(XYZ)

    # # We're not in the RGB gamut: approximate by desaturating
    # nega = np.any(rgb < 0, axis=0)  # [n, ]
    # rgb[:, nega] = rgb[:, nega] - np.min(rgb[:, nega], axis=0)
    return rgb / np.max(rgb, axis=0)


def plot_color_image(spec_cut, wavelengths, ):
    """ Given spectrum, plot a CIE color image. """

    #     wavelengths = np.array([wavel[0] for wavel in hdu1.data])
    assert spec_cut.ndim == 2
    assert spec_cut.shape[-1] == len(wavelengths)
    lams_nm = 1000. * wavelengths
    colors = spec_to_RGB(lams_nm, spec_cut, is_shift=True, is_norm=True)

    Y = XYZ_from_spec(lams_nm, spec_cut)[1, :]
    print("Y min:", Y.min())
    print("Y max:", Y.max())

    f, ax = plt.subplots(figsize=[4, 4])
    colors_swap = np.swapaxes(np.swapaxes(colors, 0, 1), 1, 2)
    colors_swap = colors_swap[::-1, :, :]
    ax.imshow(colors_swap)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{feature}_color.pdf")
    plt.savefig(f"{plot_dir}/{feature}_color.png", dpi=300)
    plt.show()
    
    return
