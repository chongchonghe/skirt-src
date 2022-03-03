# colour_system.py
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps, trapz

def gaussian_sing(x, alpha, mu, sigma1, sigma2):
    sr = (x - mu) / (sigma1 if x < mu else sigma2)
    return alpha * np.exp( -(sr * sr)/2 )

gaussian = np.vectorize(gaussian_sing)
gaussian.excluded.add((1,2,3,4))

def xyz_from_wavelength(wavelength):
    """Input: wavelength, (n,) array of wavelengths (in nm). 
    Return: (3, n) array of xyz"""

    x = gaussian(wavelength,  1.056, 5998, 379, 310) + \
        gaussian(wavelength,  0.362, 4420, 160, 267) + \
        gaussian(wavelength, -0.065, 5011, 204, 262)
    y = gaussian(wavelength,  0.821, 5688, 469, 405) + \
        gaussian(wavelength,  0.286, 5309, 163, 311)
    z = gaussian(wavelength,  1.217, 4370, 118, 360) + \
        gaussian(wavelength,  0.681, 4590, 260, 138)
    return np.stack((x, y, z))

xyz_from_wavelength


def rect(y, x):
    """ integrate y over x, which are arrays with same dimension. y can
    be multidimension arrays, and only the first dimension is integrated. """

    dx = x[1:] - x[:-1]
    yleft = y[:-1, ...]
    ret = np.dot(np.transpose(yleft), dx)  # np.dot sums over the last dimension
    return np.transpose(ret)


def xyz_from_xy(x, y):
    """Return the vector (x, y, 1-x-y)."""
    return np.array((x, y, 1-x-y))

class ColourSystem:
    """A class representing a colour system.

    A colour system defined by the CIE x, y and z=1-x-y coordinates of
    its three primary illuminants and its "white point".

    TODO: Implement gamma correction

    """

    # The CIE colour matching function for 380 - 780 nm in 5 nm intervals
    cmf = np.loadtxt('cie-cmf.txt', usecols=(1,2,3))

    def __init__(self, red, green, blue, white):
        """Initialise the ColourSystem object.

        Pass vectors (ie NumPy arrays of shape (3,)) for each of the
        red, green, blue  chromaticities and the white illuminant
        defining the colour system.

        """

        # Chromaticities
        self.red, self.green, self.blue = red, green, blue
        self.white = white
        # The chromaticity matrix (rgb -> xyz) and its inverse
        self.M = np.vstack((self.red, self.green, self.blue)).T 
        self.MI = np.linalg.inv(self.M)
        # White scaling array
        self.wscale = self.MI.dot(self.white)
        # xyz -> rgb transformation matrix
        self.T = self.MI / self.wscale[:, np.newaxis]

    def xyz_to_rgb(self, xyz, out_fmt=None, globalnorm=True):
        """Transform from xyz to rgb representation of colour.

        The output rgb components are normalized on their maximum
        value. If xyz is out the rgb gamut, it is desaturated until it
        comes into gamut.

        By default, fractional rgb components are returned; if
        out_fmt='html', the HTML hex string '#rrggbb' is returned.

        """

        # CC: make it work for ndim=2 array xyz

        # rgb = self.T.dot(xyz)
        # self.T: [3, 3]; xyz: [3, n]
        # equiv to np.dot(self.T, xyz)
        rgb = np.matmul(self.T, xyz)  # [3, n]

        # if np.any(rgb < 0):
        #     # We're not in the RGB gamut: approximate by desaturating
        #     w = - np.min(rgb)
        #     rgb += w
        # We're not in the RGB gamut: approximate by desaturating
        nega = np.any(rgb < 0, axis=0)  # [n, ]
        rgb[:, nega] = rgb[:, nega] - np.min(rgb[:, nega], axis=0)

        # if not np.all(rgb==0):
        #     # Normalize the rgb vector
        #     rgb /= np.max(rgb)
        if globalnorm:
            rgb /= rgb.max()
        else:
            notallzero = np.invert(np.all(rgb==0, axis=0))
            rgb[:, notallzero] /= np.max(rgb[:, notallzero], axis=0)

        if out_fmt == 'html':
            return self.rgb_to_hex(rgb)
        return rgb

    def rgb_to_hex(self, rgb):
        """Convert from fractional rgb values to HTML-style hex string."""

        hex_rgb = (255 * rgb).astype(int)
        return '#{:02x}{:02x}{:02x}'.format(*hex_rgb)

    def spec_to_xyz_multiD_version(self, spec, lams=None):
        """Convert a spectrum to an xyz point.

        The spectrum must be on the same grid of points as the colour-matching
        function, self.cmf: 380-780 nm in 5 nm steps.

        """

        # CC: apply to any wavelength samples instead of 380-780-5
        assert spec.ndim == 2, "try to reshape spec as n_lambd by n_pixel array"
        if lams is None:
            if spec.ndim == 1:
                XYZ = np.sum(spec[:, np.newaxis] * self.cmf, axis=0)
            else:
                XYZ = np.transpose(np.dot(np.transpose(spec), self.cmf))
        else:
            # arbitrary lambda samples
            # pick the lambdas in the range 380-780
            bigger = np.where(lams >= 380)[0]
            pick = bigger[np.where(lams[bigger] <= 780)]
            lams_pick = lams[pick]
            spec_pick = spec[pick]
            if spec.ndim == 1:
                XYZ = np.zeros(3)
                for i in range(3):
                    lams0 = np.arange(380, 781, 5)
                    f = interp1d(lams0, self.cmf[:, i], kind='cubic')
                    cmf = f(lams_pick)
                    # integrate cmf over lams_pick
                    # XYZ[i] = np.sum(spec_pick * cmf)   # only when lams_pick are evenly spaced
                    XYZ[i] = trapz(spec_pick * cmf, lams_pick)
            else:
                XYZ = np.zeros([3] + list(spec.shape[1:]))
                for i in range(3):
                    lams0 = np.arange(380, 781, 5)
                    f = interp1d(lams0, self.cmf[:, i], kind='cubic')
                    cmf = f(lams_pick)

                    # integrate cmf over lams_pick
                    cmf_tile = np.tile(cmf, spec_pick.shape[1:][::-1] + (1,))
                    cmf_tile = np.transpose(cmf_tile)
                    XYZ[i] = rect(spec_pick * cmf_tile, lams_pick)

                    # # the following may replace the three lines above
                    # product = np.transpose(spec_pick.transpose() * cmf)
                    # XYZ[i] = rect(product, lams_pick)

        # print('XYZ =', XYZ)
        den = np.sum(XYZ, axis=0)
        if spec.ndim == 1:
            if den == 0.:
                return XYZ
            return XYZ / den
        else:
            # nonzero = den != 0.
            den[den == 0] = 1.
            dim = (3,) + (1,) * den.ndim
            XYZ = XYZ / np.tile(den, dim)
            return XYZ

    def spec_to_xyz(self, spec, lams=None):
        """Convert a spectrum to an xyz point. spec has ndim = 2

        The spectrum must be on the same grid of points as the colour-matching
        function, self.cmf: 380-780 nm in 5 nm steps.

        """

        # CC: apply to any wavelength samples instead of 380-780-5
        if lams is None:
            XYZ = np.transpose(np.dot(np.transpose(spec), self.cmf))
        else:
            # arbitrary lambda samples
            # pick the lambdas in the range 380-780
            n_pixel = spec.shape[1]
            bigger = np.where(lams >= 380)[0]
            pick = bigger[np.where(lams[bigger] <= 780)]  # (n_pick,)
            lams_pick = lams[pick]                        # (n_pick,)
            spec_pick = spec[pick]                        # (n_pick, n_pixel)
            XYZ = np.zeros([3, n_pixel])                  # (3, n_pixel)
            for i in range(3):
                lams0 = np.arange(380, 781, 5)
                f = interp1d(lams0, self.cmf[:, i], kind='cubic')
                cmf = f(lams_pick)
                # integrate cmf over lams_pick
                product = np.transpose(spec_pick.transpose() * cmf)  # (n_pick, n_pixel)
                XYZ[i] = rect(product, lams_pick)

        # print('XYZ =', XYZ)
        den = np.sum(XYZ, axis=0)  # (n_pixel,)
        den[den == 0] = 1.
        XYZ *= 1. / den
        return XYZ

    def spec_to_rgb(self, spec, lams=None, out_fmt=None, globalnorm=True):
        """Convert a spectrum to an rgb value."""

        assert spec.ndim == 2, "try to reshape spec as n_lambd by n_pixel array"
        xyz = self.spec_to_xyz(spec, lams=lams)
        return self.xyz_to_rgb(xyz, out_fmt, globalnorm=globalnorm)
            

illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
cs_hdtv = ColourSystem(red=xyz_from_xy(0.67, 0.33),
                       green=xyz_from_xy(0.21, 0.71),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)

cs_smpte = ColourSystem(red=xyz_from_xy(0.63, 0.34),
                        green=xyz_from_xy(0.31, 0.595),
                        blue=xyz_from_xy(0.155, 0.070),
                        white=illuminant_D65)

cs_srgb = ColourSystem(red=xyz_from_xy(0.64, 0.33),
                       green=xyz_from_xy(0.30, 0.60),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)

