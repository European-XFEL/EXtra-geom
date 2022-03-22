"""
The following is the attempt to write a module for the dssc detector, that is
similar to the extra-geometry package, but includes the hex-cart conversion
that Andreas Scherz wrote in 2020. It should eventually be merged into the
aforementioned package.
"""
import joblib
import numpy as np
import math as math
from skimage.util import pad as skpad
from scipy.ndimage import gaussian_filter1d


# -----------------------------------------------------------------------------
# helper methods
# -----------------------------------------------------------------------------
def is_even(num):
    """
    is num an even number?
    :param num: integer
    :return: boolean True = num is even number
    """
    return not (num & 1)


def is_odd(num):
    """
    is num an even number?
    :param num: integer
    :return: boolean True = num is even number
    """
    return num & 1


def next_power_of_2(x=1):
    """
    Returns the next higher power of 2 compared to x
    :param x: current dimension
    :return: next higher power of 2
    """
    return 1 << (x - 1).bit_length()


def get_im_center(im: np.array) -> (int, int):
    """
    the center of image is the DC components with respect to discrete Fourier
    transformation
    :returns center pos (int):
    """
    return int(im.shape[0] / 2), int(im.shape[1] / 2)


def pad4fft(im: np.array = None,
            dim: [None, bool, int] = None,
            fill: [None, str, float] = None) -> np.array:
    """
    pad image using dim parameter and pad image with given fill
    :param im:
    :param dim: None == up to next power of 2, True/False == oversampling, or
    int with given dimension
    :param fill: None == zeros, str == special format see skimage.pad, float ==
    fill value
    :return: object itself
    """
    if im is not None:
        im = np.array(im).squeeze()
    else:
        im = skdata.lena()
    # print('image size is ' + str(im.shape))
    # print('dimension is supposed to be ' + str(dim))

    if dim is None:
        dim = max((next_power_of_2(im.shape[0]), next_power_of_2(im.shape[1])))
        dim = (dim, dim)
    elif type(dim) is bool:
        oversampling = dim
        print('is boolean option')
        dim = max((next_power_of_2(im.shape[0]), next_power_of_2(im.shape[1])))
        dim *= 2 if oversampling else 1  # another factor 2 for oversampling
        dim = (dim, dim)
    elif (dim[0] < im.shape[0]) and (dim[1] < im.shape[1]):
        return im

    # print('dimension is ' + str(dim))

    w0 = ((dim[0] - im.shape[0]) // 2)
    w1 = ((dim[1] - im.shape[1]) // 2)

    width0 = (w0, w0) if is_even(im.shape[0]) else (w0, w0 + 1)
    width1 = (w1, w1) if is_even(im.shape[1]) else (w1, w1 + 1)

    if fill is None:
        im = skpad(im, (width0, width1), mode='constant', constant_values=0)
    elif type(fill) is str:
        im = skpad(im, (width0, width1), mode=fill)
    else:
        im = skpad(im, (width0, width1), mode='constant', constant_values=fill)

    return im


def n_choose_k(n, k):
    """
    :n choose k problem: returns the number of ways to choose a subset k from a
    set of n elements this is eqivalent to:
        factorial(n)/factorial(k)/factorial(n-k)

    written by A. Scherz
    """
    return math.factorial(n) / math.factorial(k) / math.factorial(n - k)


# -----------------------------------------------------------------------------
# Thiran Filter class
# -----------------------------------------------------------------------------
def thiran_type_I(order=3):
    return Thiran(1, order)


def thiran_type_II(order=2):
    return Thiran(2, order)


class Thiran:
    """
    there are two types of 1D fractional filters based on Thiran filters. both
    types can have an order N=1 to inf. some types, provide simple algorithms
    for applying the delay to a signal. for type II and orders larger than N.
    he transfer and impulse response fct have to be calculated and are more
    computational expensive
    The generator defaults the Thiran filter of type I and order N=3
    For more details see of the filters and algorithms,
        see pg. 685, IEEE Trans. on Im. Proc. 17, 679 (2008)

    written by A. Scherz, XFEL, 2012-11-26
    updated by A. Scherz, XFEL, 2019-02-11 for python
    added to SpeckleToolbox by A. Scherz, XFEL, 2020-08-06
    """

    type_I = staticmethod(thiran_type_I)

    type_II = staticmethod(thiran_type_II)

    def __init__(self, filter_type=1, order=3, delay=0):
        self._type = filter_type
        self._n = order
        self._delay = delay
        self._b = np.zeros(order)
        self._a_minus_delay = 0
        self._a_plus_delay = 0

    def __str__(self):
        print('Thiran 1D Fractional filter')
        print('  Type : ' + self._type.__str__())
        print(' Order : ' + self._n.__str__())
        print(' Delay : ' + self._delay.__str__())

        return ''

    def apply_delay(self, signal):
        # periodic boundaries are important. add padding of at least the order
        # N left and right
        # this should be applied before calling this method for type I
        shifted_signal = np.array(signal)

        if self._type == 1:     # type I Thiran
            for j in range(len(shifted_signal) - self._n - 1, self._n - 1, -1):
                for k in range(1, self._n + 1):
                    shifted_signal[j] += self._b[k - 1] \
                        * (shifted_signal[j - k] - shifted_signal[j + k])
        else:   # type II of order N = 2
            for j in range(self._n - 1, len(shifted_signal) - self._n - 1, 1):
                shifted_signal[j] += self._a_plus_delay \
                    * (shifted_signal[j - 1] - shifted_signal[j + 1])
            for j in range(len(shifted_signal) - self._n - 1, self._n - 1, -1):
                shifted_signal[j] += self._a_minus_delay \
                    * (shifted_signal[j + 1] - shifted_signal[j - 1])

        return shifted_signal  # shifted signal

    def delay(self, tau):
        if (tau < -0.5) or (tau > 0.5):
            raise ValueError('WARNING out of range: delay: 0.0 - 0.5')
        self._delay = tau
        # update coefficients for delay
        self._coefficients()
        return self

    def _coefficients(self):
        for k in range(1, self._n + 1):
            p_k = 1
            for l in range(self._n + 1):
                p_k *= (self._delay - l) / (self._delay - l - k)
            # n_choose_k function is factorial(n)/factorial(k)/factorial(n-k)
            self._b[k - 1] = (-1) ** k * p_k * n_choose_k(self._n, k)
        self._a_plus_delay = \
            (-4 + self._delay**2 + np.sqrt(12 - 3*self._delay**2))
        self._a_minus_delay = self._a_plus_delay
        self._a_plus_delay /= (3*self._delay + 2 + self._delay**2)
        self._a_minus_delay /= (-3*self._delay + 2 + self._delay**2)


# -----------------------------------------------------------------------------
# shear methods using the Thiran filter
# -----------------------------------------------------------------------------
def shear_im(im: np.array, delay: float,
             im_axis: int = 1,
             thiran: Thiran = thiran_type_I()) -> np.array:
    """
    : direction (str): 1=='hor' or 0=='ver', required
    : thiran: Thiran type 1 of order 3 is default (default)
    """
    im = np.float64(im)
    cx, cy = get_im_center(im)
    # print(im.shape)
    # print(cx, cy)
    if im_axis == 0:
        for k in range(2 * cy):
            shift = delay * (k - cy)
            # full pixel shifts
            im[:, k] = np.roll(im[:, k], round(shift))
            # now shift a pixel fraction
            im[:, k] = thiran.delay(shift-round(shift)).apply_delay(im[:, k])
            im[:, k] = gaussian_filter1d(im[:, k], 0.6)
            #  thiran_filter(im[:, k], thiran(order, shift - round(shift)))
        # print('apply vertical shear')
    elif im_axis == 1:
        for k in range(2 * cx):
            shift = delay * (k - cx)
            # full pixel shifts
            im[k, :] = np.roll(im[k, :], round(shift))
            # now shift a pixel fraction: shift between +/- [0-0.5]
            im[k, :] = thiran.delay(shift - round(shift)).apply_delay(im[k, :])
            im[k, :] = gaussian_filter1d(im[k, :], 0.6)
            # thiran_filter(im[k, :], thiran(order, shift - round(shift)))
        # print('apply horizontal shear')
    else:
        # TODO: throw exception
        print('direction undefined')

    return im


def convert_hex2cart(im_hex, thiran=thiran_type_II()):
    """
    [ imH ] = ConversionCart2Hex( imC, orderN, showImage )
    imC: image in cartesian coord.
    thiran: Thiran filter type 1 of order 3 is default
    show: optional, if true an image representation in hex is shown.
    note: show hex image may take considerable amount of time for larger
    images!)
    Nov.28, 2012, A. Scherz
    """
    # first step: zero pad the image for extending boundaries

    im = pad4fft(im_hex,
                 dim=[2*im_hex.shape[0], 2*im_hex.shape[1]],
                 fill='reflect')
    cx, cy = get_im_center(im)

    # need to skew image to hexagonal shape before apply shear operations
    height = im_hex.shape[0]
    y0 = height % 2
    for y in range(- height // 2, height // 2 + y0 + 1):
        im[y - y0 + cx, :] = np.roll(im[y - y0 + cx, :],
                                     - math.floor((y + y0) / 2))

    sqrt3 = math.sqrt(3)
    shear1 = sqrt3 - math.sqrt(6 / sqrt3)
    shear2 = math.sqrt(sqrt3 / 6)
    shear3 = 2 - math.sqrt(6 / sqrt3)
    # this is conversion from cartesian to hexagonal lattice
    im = shear_im(im, shear3, im_axis=0, thiran=thiran)
    im = shear_im(im, shear2, im_axis=1, thiran=thiran)
    im = shear_im(im, shear1, im_axis=0, thiran=thiran)

    # the image should now form a rectangular box of following size:
    height = int(np.round(hex2cart([im_hex.shape[0], 0]))[0])
    width = int(np.round(hex2cart([0, im_hex.shape[1]]))[1])
    y0 = height % 2

    im = im[-height // 2 + cx: height // 2 + cx + y0,
            -width // 2 + cy: width // 2 + cy + 1]

    return im


def hex2cart(r_hex):
    sqrt3 = math.sqrt(3)
    rhex = math.sqrt(2 / sqrt3) * np.array([[1, 1 / 2], [0, sqrt3 / 2]])
    rcart = np.linalg.inv(rhex)

    return np.dot(rcart, r_hex)


# -----------------------------------------------------------------------------
# The DSSC assembler class
# -----------------------------------------------------------------------------
class DSSCAssembler:
    """
    This class splits the dssc detector modules into tiles, performs a
    hex-cart coordinate transformation over all 32 tiles in parallel, and puts
    the tiles in the specified location.

    Notes
    -----
    - The module is not quite finished yet. It is missing high-level methods to
    make it more compatible with the extra-geometry package (conversion of flat
    coordinate structure into the nested one in the aforementioned package).
    - There is a high frequency ringing originating from the conversion. This
    should eventually be dealt with by applying appropriate hp-filters in fs.
    At the moment I just apply a gaussian filter in fs, which is not quite
    correct.
    """
    generic_parameters = {
        # module | (x,y) | orientation 
        '0': [[571.5, 294.0], 'lr'],
        '1': [[571.5, 11.0], 'lr'],
        '2': [[710.5, 295.0], 'lr'],
        '3': [[710.5, 11.0], 'lr'],
        '4': [[849.5, 295.0], 'lr'],
        '5': [[849.5, 11.0], 'lr'],
        '6': [[987.0, 295.0], 'lr'],
        '7': [[987.5, 11.0], 'lr'],
        '8': [[9.0, 282.5], 'lr'],
        '9': [[9.0, 0.0], 'lr'],
        '10': [[147.0, 283.5], 'lr'],
        '11': [[147.0, 0.0], 'lr'],
        '12': [[287.0, 284.0], 'lr'],
        '13': [[287.0, 1.0], 'lr'],
        '14': [[425.5, 284.0], 'lr'],
        '15': [[426.5, 1.0], 'lr'],
        '16': [[417.5, 579.0], 'ud'],
        '17': [[417.5, 862.0], 'ud'],
        '18': [[278.5, 579.0], 'ud'],
        '19': [[278.5, 862.0], 'ud'],
        '20': [[139.5, 579.0], 'ud'],
        '21': [[139.5, 862.0], 'ud'],
        '22': [[0.0, 578.0], 'ud'],
        '23': [[0.5, 861.0], 'ud'],
        '24': [[979.0, 589.0], 'ud'],
        '25': [[979.0, 872.0], 'ud'],
        '26': [[840.5, 588.0], 'ud'],
        '27': [[840.0, 872.0], 'ud'],
        '28': [[701.0, 588.0], 'ud'],
        '29': [[701.0, 872.0], 'ud'],
        '30': [[562.5, 588.0], 'ud'],
        '31': [[562.5, 871.0], 'ud']
    }

    @staticmethod
    def find_limits(pos_list):
        """
        Finds maximum parameters in parameter list. Used to find image
        dimensions of image.

        Parameters
        ----------
        pos_list: list
            A dictionary containing the positions of all tiles.

        Returns
        -------
        lim_x: list
        lim_y: list
        """
        vals_x = [x for x, _ in pos_list]
        vals_y = [y for _, y in pos_list]
        lim_x = [np.min(vals_x), np.max(vals_x)]
        lim_y = [np.min(vals_y), np.max(vals_y)]
        return lim_x, lim_y

    @staticmethod
    def read_positions(geometry_parameters):
        """
        Convert the parameter list in the dictionary to a list

        Parameters
        ----------
        geometry_parameters: dict

        Returns
        -------
        positions: list
        """
        positions = []
        for key in geometry_parameters.keys():
            positions.append(
                geometry_parameters[key][0])
        return positions

    @staticmethod
    def read_geometry_file(filename):
        """
        Read and set internal parameter list from hdf5 file.

        Parameters
        ----------
        filename

        Notes
        -----
        This method could be used to make the class fully compatible with the
        extra-geometry package.
        """
        raise NotImplementedError

    @staticmethod
    def write_positions(parameter_dic, position_list):
        """
        Fill internal parameter dictionary.

        Parameters
        ----------
        parameter_dic: dict
            old dictionary
        position_list: list
            list to be inserted into dictionary

        Returns
        -------
        dict
            new dictionary to be used in DSSCAssembler class.
        """
        dic_new = dict(parameter_dic)
        for key, ind in zip(dic_new.keys(), range(len(position_list))):
            dic_new[key][0] = position_list[ind]
        return dic_new

    @staticmethod
    def write_geometry_file(filename, params):
        """
        Write parameters into hdf5 file.

        Parameters
        ----------
        filename: str
        params: list

        Notes
        -----
        This method could be used to make the class fully compatible with the
        extra-geometry package by saving the parameters in the same format as
        implemented within the latter.
        """
        raise NotImplementedError

    @staticmethod
    def shift_tile(tile, delay, axis=0):
        """
        Apply sub pixel shift

        Parameters
        ----------
        tile: numpy.ndarray
        delay: float
        axis: {'0', '1'}, optional

        Returns
        -------
        arr_out: numpy.ndarray
        """
        arr_out = np.copy(tile)
        if delay == 0:
            return arr_out

        sh = np.shape(tile)
        thiran = thiran_type_II()
        thiran.delay(delay)

        if axis == 0:
            for i in range(sh[0]):
                arr_out[i, :] = thiran.apply_delay(arr_out[i, :])
        elif axis == 1:
            for i in range(sh[1]):
                arr_out[:, i] = thiran.apply_delay(arr_out[:, i])

        return arr_out

    @staticmethod
    def split_module(module_hex):
        """
        Return individual tiles of module

        Parameters
        ----------
        module_hex: numpy.ndarray

        Returns
        -------
        tile_hex: list, [numpy.ndarray, numpy.ndarray]
        """
        tile_hex = np.zeros((2, 128, 256), dtype='float64')
        tile_hex[0, :, :] = module_hex[:, 0:256]
        tile_hex[1, :, :] = module_hex[:, 256:512]
        return tile_hex

    @classmethod
    def clean_geometry_dict(cls, geometry_parameters):
        """
        Shape parameter list, such that there are no empty spaces around image.

        Parameters
        ----------
        geometry_parameters: dict
        """
        positions = cls.read_positions(geometry_parameters)
        [limx, limy] = cls.find_limits(positions)
        for i in range(len(positions)):
            positions[i][0] -= limx[0]
            positions[i][1] -= limy[0]
        return cls.write_positions(geometry_parameters, positions)

    @classmethod
    def shift_quadrant(cls, param_dict, q, shift):
        """
        Shift all tiles of a quadrant individually.

        Parameters
        ----------
        param_dict: dict
        q: int, {0, 1, 2, 3}
            quadrant to be shifted
        shift: list, [float, float]
            amount by which the quadrant should be shifted

        Returns
        -------
        positions: dict
            updated parameter dictionary
        """
        quadrants = {
            '0': range(0, 8),
            '1': range(8, 16),
            '2': range(16, 24),
            '3': range(24, 32)
        }
        positions = cls.read_positions(param_dict)
        for i in quadrants[str(q)]:
            positions[i][0] += shift[0]
            positions[i][1] += shift[1]
        return cls.write_positions(param_dict, positions)

    def __init__(self, geom_params):
        """
        Constructor of DSSCAssembler class. Similarly to the extra-geometry
        package one usually constructs an object of this class by the
        respective classmethods.

        Parameters
        ----------
        geom_params: dict
            dicionary indicating the position and orientation of all tiles.
        """
        self.geom_params = geom_params

    @classmethod
    def convert_assemble_single(cls, data, custom_geometry=None):
        """
        construct converted image from raw DSSC-detector data and given
        geometry file.

        Parameters
        ----------
        data: numpy.ndarray
            raw detector data
        custom_geometry: dict
            parameter dictionary

        Returns
        -------
        numpy.ndarray:
            image in cartesian coordinates
        """
        if custom_geometry:
            dssc_obj = cls(geom_params=custom_geometry)
        else:
            dssc_obj = cls(geom_params=cls.generic_parameters)

        data_conv = dssc_obj.convert_all(data)
        return dssc_obj.assemble_image(data_conv)

    @classmethod
    def assemble_single(cls, data,
                        custom_geometry=None):
        """
        construct image from raw DSSC-detector data and given
        geometry file.

        Parameters
        ----------
        data: numpy.ndarray
            raw detector data
        custom_geometry: dict
            parameter dictionary
        """
        raise NotImplementedError

    def assemble_image(self, tiles_conv):
        positions = self.read_positions(self.geom_params)
        [limx, limy] = self.find_limits(positions)
        s_x, s_y = int(np.floor(limx[1])), int(np.floor(limy[1]))
        assembled_image = np.zeros((s_x+120, s_y+276))
        for i in range(32):
            [s_x, s_y] = self.geom_params[str(i)][0]
            d_x, s_x = float(s_x - np.floor(s_x)), int(np.floor(s_x))
            d_y, s_y = float(s_y - np.floor(s_y)), int(np.floor(s_y))
            tile = tiles_conv[i]
            tile = self.shift_tile(tile, d_x, 
                                   axis=1)
            tile = self.shift_tile(tile, d_y, 
                                   axis=0)

            assembled_image[s_x:s_x+120, s_y:s_y+276] = tile
        return assembled_image

    @classmethod
    def assemble_from_file(cls, data, geometry_file):
        """
        Assemble image from parameters given in hdf5 file.

        Parameters
        ----------
        data: numpy.ndarray
            raw detector data
        geometry_file: str
            path to geometry file
        """
        raise NotImplementedError

    def convert_all(self, hex_data):
        """
        Split and convert all modules.

        Parameters
        ----------
        hex_data: numpy.ndarray
            raw detector data

        Returns
        -------
        list:
            list containing data of 32 tiles.
        """
        tiles = []
        for i in range(16):
            split = self.split_module(hex_data[i])
            tiles.append(split[0])
            tiles.append(split[1])

        cart_data = joblib.Parallel(n_jobs=32) \
            (joblib.delayed(self.convert_single_tile)
                (n, tiles[n]) for n in range(32))

        return cart_data

    def convert_single_tile(self, m, tile_hex):
        """
        Convert a single tile.

        Parameters
        ----------
        m: int
            tile number
        tile_hex: numpy.ndarray
            raw data of single tile

        Returns
        -------
        tile_cart: numpy.ndarray
            tile in cartesian coordinates
        """
        tile_cart = convert_hex2cart(tile_hex)
        key = str(m)

        if self.geom_params[key][1] == 'lr':
            tile_cart = np.flip(tile_cart, axis=1)
        elif self.geom_params[key][1] == 'ud':
            tile_cart = np.flip(tile_cart, axis=0)

        return tile_cart
