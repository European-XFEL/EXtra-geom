"""Detector geometry handling."""
import warnings
from itertools import product
from typing import List, Tuple

import h5py
import numpy as np

from .base import DetectorGeometryBase, GeometryFragment


class GenericGeometry(DetectorGeometryBase):
    """A generic detector layout based either on the CrystFEL geom file or on a set of parameters.

    The coordinates used in this class are 3D (x, y, z), and represent metres.

    The :attr:`expected_data_shape` is a following triple :

        1. the number of modules which is the length of the :attr:`corner_coordinates` array
        2. the number of tiles in a module along the slow-scan direction multiplied by
           the number of slow-scan pixels per tile :attr:`frag_ss_pixels`
        3. the number of tiles in a module along the fast-scan direction multiplied by
           the number of fast-scan pixels per tile :attr:`frag_fs_pixels`
    """
    detector_type_name = 'Generic Detector'

    @classmethod
    def from_simple_description(cls, pixel_size: float, slow_pixels: int, fast_pixels: int,
                                corner_coordinates: List[np.ndarray] = [np.zeros(3)],
                                ss_vec: np.ndarray = np.array([1, 0, 0]),
                                fs_vec: np.ndarray = np.array([0, 1, 0]),
                                n_tiles_per_module: int = 1,
                                tile_gap: float = None,
                                tile_vec: np.ndarray = None,
                                ):
        """ Creates a generic detector from a dictionary.

        Parameters
        __________

        pixel_size: float
            the size of a pixel in meters (reversed CrystFEL's `res`)
        slow_pixels, fast_pixels: int
            the size of a tile along the slow- and the fast-scan axes
        corner_coordinates: ndarray
            3D coordinates of the first pixel of each module
        ss_vec, fs_vec: ndarray
            3D unit vectors of the slow- and the fast-scan directions in the lab coordinates (the X-axis
            points to the left looking along the beam, the Y-axis points up, and the Z-axis goes with the beam).
            Example: np.array([0, 1, 0])
        n_tiles_per_module: int
            the number of tiles in each module, default=1
        tile_gap: float
            the gap between two tiles in metres, default=pixel_size
        tile_vec: ndarray
            the direction of tile replication, default=[1, 0, 0]

        """

        modules = []

        tile_vec = np.array(tile_vec) if tile_vec else ss_vec

        # Get the tile shift: it is a multiple of either `fast_pixels` or `slow_pixels`
        tile_offset_value = np.abs(np.inner(fs_vec * fast_pixels + ss_vec * slow_pixels, tile_vec))
        tile_gap = tile_gap if tile_gap else pixel_size

        for corner_coord in corner_coordinates:
            module = []
            for t in range(n_tiles_per_module):
                tile = GeometryFragment(corner_coord +
                                        tile_vec * t * (tile_offset_value * pixel_size + tile_gap),
                                        ss_pixels=slow_pixels, fs_pixels=fast_pixels,
                                        ss_vec=ss_vec * pixel_size,
                                        fs_vec=fs_vec * pixel_size)
                module += [tile]
            modules += [module]

        geom = cls(modules)

        geom.pixel_size = pixel_size
        geom.frag_fs_pixels = fast_pixels
        geom.frag_ss_pixels = slow_pixels
        geom.n_modules = len(corner_coordinates)
        geom.n_tiles_per_module = n_tiles_per_module

        assert np.inner(fs_vec, ss_vec) == 0    # scan vectors are perpendicular
        assert np.linalg.norm(fs_vec) == np.linalg.norm(ss_vec) == 1    # unit vectors

        # the numbers of tiles per module in the fast- and slow-scan directions respectively:
        geom.fs_tiles = abs(np.inner(fs_vec, tile_vec)) * n_tiles_per_module or 1
        geom.ss_tiles = abs(np.inner(ss_vec, tile_vec)) * n_tiles_per_module or 1

        geom.expected_data_shape = (geom.n_modules,
                                    geom.ss_tiles * geom.frag_ss_pixels,
                                    geom.fs_tiles * geom.frag_fs_pixels)

        return geom

    def inspect(self, axis_units='px', frontview=True, aspect='auto'):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------

        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        aspect : str, float
          Set the aspect ratio of the plot, possible values:
            * 'auto' (default): automatic; fill the position rectangle with data
            * 'equal': same scaling from data to plot units for x and y
            * a number: a figure will be stretched such that the height is num times the width.
              aspect=1 is the same as 'equal'.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)

        ax.set_title(f"{self.detector_type_name} geometry ({self.filename})")

        scale = self._get_plot_scale_factor(axis_units)

        # Label modules and tiles
        for ch, module in enumerate(self.modules):
            label = f"M{ch}" if len(self.modules) > 1 else ""
            if self.n_tiles_per_module > 1:
                for t in [0, self.n_tiles_per_module - 1]:
                    cx, cy, _ = module[t].centre() * scale
                    ax.text(cx, cy, label + f"T{t}",
                            verticalalignment='center',
                            horizontalalignment='center')
            elif label:
                # one tile per each of multiple modules: only modules are labelled
                cx, cy, _ = module[0].centre() * scale
                ax.text(cx, cy, label, verticalalignment='center', horizontalalignment='center')
        ax.set_aspect(aspect)
        return ax

    def _tile_slice(self, tileno: int) -> Tuple[slice]:
        # Since python 3.9 it is legal to annotate the output simply as  `-> tuple[slice]`
        """ Which part of the data array is this tile?"""
        if self.fs_tiles > 1:
            fs_slice = slice(tileno * self.frag_fs_pixels, (tileno + 1) * self.frag_fs_pixels)
        else:
            fs_slice = slice(0, self.frag_fs_pixels)
        if self.ss_tiles > 1:
            ss_slice = slice(tileno * self.frag_ss_pixels, (tileno + 1) * self.frag_ss_pixels)
        else:
            ss_slice = slice(0, self.frag_ss_pixels)
        return ss_slice, fs_slice

    def split_tiles(self, module_data: np.ndarray) -> List[np.ndarray]:
        # Since python 3.9 it is legal to annotate the output simply as  `-> list[np.ndarray]`

        if self.n_tiles_per_module == 1:
            return [module_data]
        elif self.ss_tiles > 1:
            return [module_data[..., s * self.frag_ss_pixels: (s + 1) * self.frag_ss_pixels, :]
                    for s in range(self.n_tiles_per_module)]
        else:
            return [module_data[..., s * self.frag_fs_pixels: (s + 1) * self.frag_fs_pixels]
                    for s in range(self.n_tiles_per_module)]

    def _module_coords_to_tile(self, slow_scan: np.ndarray, fast_scan: np.ndarray):
        """Positions in module to tile numbers & pos in tile.

        `slow_scan` and `fast_scan` are arrays of equal size, they contain
        the coordinates of points of interest along the slow and the fast scan axes
        respectively.

        Returned values are three arrays, and they contain:
        1. the number of the tile in a module the point belongs to,
        2. the slow-scan axis coordinate within the tile, i.e. `mod(ss, frag_ss_pixel)`
        3. the fast-scan axis coordinate within the tile, i.e. `mod(fs, frag_fs_pixel)`
        """
        if self.n_tiles_per_module == 1:
            return 0, slow_scan, fast_scan
        elif self.ss_tiles > 1:
            tileno, tile_ss = np.divmod(slow_scan, self.frag_ss_pixels)
            return tileno.astype(np.int16), tile_ss, fast_scan
        else:
            tileno, tile_fs = np.divmod(fast_scan, self.frag_fs_pixels)
            return tileno.astype(np.int16), slow_scan, tile_fs

    @classmethod
    def from_crystfel_geom(cls, filename):
        raise NotImplementedError


class AGIPD_1MGeometry(DetectorGeometryBase):
    """Detector layout for AGIPD-1M

    The coordinates used in this class are 3D (x, y, z), and represent metres.

    You won't normally instantiate this class directly:
    use one of the constructor class methods to create or load a geometry.
    """
    detector_type_name = 'AGIPD-1M'
    pixel_size = 2e-4  # 2e-4 metres == 0.2 mm
    frag_ss_pixels = 64
    frag_fs_pixels = 128
    expected_data_shape = (16, 512, 128)
    n_quads = 4
    n_modules = 16
    n_tiles_per_module = 8

    @classmethod
    def from_quad_positions(cls, quad_pos, asic_gap=2, panel_gap=29,
                            unit=pixel_size):
        """Generate an AGIPD-1M geometry from quadrant positions.

        This produces an idealised geometry, assuming all modules are perfectly
        flat, aligned and equally spaced within their quadrant.

        The quadrant positions are given in pixel units, referring to the first
        pixel of the first module in each quadrant, corresponding to data
        channels 0, 4, 8 and 12.

        The origin of the coordinates is in the centre of the detector.
        Coordinates increase upwards and to the left (looking along the beam).

        To give positions in units other than pixels, pass the *unit* parameter
        as the length of the unit in metres.
        E.g. ``unit=1e-3`` means the coordinates are in millimetres.
        """
        asic_gap_px = asic_gap * unit / cls.pixel_size
        panel_gap_px = panel_gap * unit / cls.pixel_size

        # How much space one tile takes up, including the gaps
        # separating it from its neighbour.
        # In the y dimension, 128 px + gap between modules
        module_height = (cls.frag_fs_pixels + panel_gap_px) * cls.pixel_size
        # In x, 64 px + gap between tiles (asics)
        tile_width = (cls.frag_ss_pixels + asic_gap_px) * cls.pixel_size

        quads_x_orientation = [1, 1, -1, -1]
        quads_y_orientation = [-1, -1, 1, 1]
        modules = []
        for p in range(16):
            quad = p // 4
            quad_corner = quad_pos[quad]
            x_orient = quads_x_orientation[quad]
            y_orient = quads_y_orientation[quad]
            p_in_quad = p % 4
            corner_y = (quad_corner[1] * unit)\
                       - (p_in_quad * module_height)

            tiles = []
            modules.append(tiles)

            for a in range(8):
                corner_x = (quad_corner[0] * unit)\
                           + x_orient * tile_width * a
                tiles.append(GeometryFragment(
                    corner_pos=np.array([corner_x, corner_y, 0.]),
                    ss_vec=np.array([x_orient, 0, 0]) * unit,
                    fs_vec=np.array([0, y_orient, 0]) * unit,
                    ss_pixels=cls.frag_ss_pixels,
                    fs_pixels=cls.frag_fs_pixels,
                ))
        return cls(modules)

    def quad_positions(self):
        """Retrieve the coordinates of the first pixel in each quadrant

        The coordinates returned are 2D and in pixel units, compatible with
        :meth:`from_quad_positions`.
        """
        return np.array([
            self.modules[q * 4][0].corner_pos[:2] for q in range(4)
        ]) / self.pixel_size

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------

        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        # Label modules and tiles
        for ch, module in enumerate(self.modules):
            s = 'Q{Q}M{M}'.format(Q=(ch // 4) + 1, M=(ch % 4) + 1)
            cx, cy, _ = module[4].centre() * scale
            ax.text(cx, cy, s, fontweight='bold',
                    verticalalignment='center',
                    horizontalalignment='center')

            for t in [0, 7]:
                cx, cy, _ = module[t].centre() * scale
                ax.text(cx, cy, 'T{}'.format(t + 1),
                        verticalalignment='center',
                        horizontalalignment='center')

        ax.set_title('AGIPD-1M detector geometry ({})'.format(self.filename))
        return ax

    def position_modules_interpolate(self, data):
        """Assemble data from this detector according to where the pixels are.

        This performs interpolation, which is very slow.
        Use :meth:`position_modules_fast` to get a pixel-aligned approximation
        of the geometry.

        Parameters
        ----------

        data : ndarray
          The three dimensions should be channelno, pixel_ss, pixel_fs
          (lengths 16, 512, 128). ss/fs are slow-scan and fast-scan.

        Returns
        -------
        out : ndarray
          Array with the one dimension fewer than the input.
          The last two dimensions represent pixel y and x in the detector space.
        centre : ndarray
          (y, x) pixel location of the detector centre in this geometry.
        """
        from scipy.ndimage import affine_transform
        assert data.shape == (16, 512, 128)
        size_yx, centre = self._get_dimensions()
        tmp = np.empty((16 * 8,) + size_yx, dtype=data.dtype)

        for i, (module, mod_data) in enumerate(zip(self.modules, data)):
            tiles_data = np.split(mod_data, 8)
            for j, (tile, tile_data) in enumerate(zip(module, tiles_data)):
                # We store (x, y, z), but numpy indexing, and hence affine_transform,
                # work like [y, x]. Rearrange the numbers:
                fs_vec_yx = tile.fs_vec[:2][::-1]
                ss_vec_yx = tile.ss_vec[:2][::-1]

                # Offset by centre to make all coordinates positive
                corner_pos_yx = tile.corner_pos[:2][::-1] + centre

                # Make the rotation matrix
                rotn = np.stack((ss_vec_yx, fs_vec_yx), axis=-1)

                # affine_transform takes a mapping from *output* to *input*.
                # So we reverse the forward transformation.
                transform = np.linalg.inv(rotn)
                offset = np.dot(rotn, corner_pos_yx)  # this seems to work, but is it right?

                affine_transform(
                    tile_data,
                    transform,
                    offset=offset,
                    cval=np.nan,
                    output_shape=size_yx,
                    output=tmp[i * 8 + j],
                )

        # Silence warnings about nans - we expect gaps in the result
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out = np.nanmax(tmp, axis=0)

        return out, centre

    def _get_dimensions(self):
        """Calculate appropriate array dimensions for assembling data.

        Returns (size_y, size_x), (centre_y, centre_x)
        """
        corners = []
        for module in self.modules:
            for tile in module:
                corners.append(tile.corners())
        corners = np.concatenate(corners)[:, :2] / self._pixel_shape

        # Find extremes, add 1 px margin to allow for rounding errors
        min_xy = corners.min(axis=0).astype(int) - 1
        max_xy = corners.max(axis=0).astype(int) + 1

        size = max_xy - min_xy
        centre = -min_xy
        # Switch xy -> yx
        return tuple(size[::-1]), centre[::-1]

    @staticmethod
    def split_tiles(module_data):
        # Split into 8 tiles along the slow-scan axis
        # This simple slicing is faster than np.split().
        return [module_data[..., s:s+64, :] for s in range(0, 512, 64)]

    @classmethod
    def _tile_slice(cls, tileno):
        # Which part of the array is this tile?
        # tileno = 0 to 7
        tile_offset = tileno * cls.frag_ss_pixels
        ss_slice = slice(tile_offset, tile_offset + cls.frag_ss_pixels)
        fs_slice = slice(0, cls.frag_fs_pixels)  # Every tile covers the full 128 pixels
        return ss_slice, fs_slice

    @classmethod
    def _module_coords_to_tile(cls, slow_scan, fast_scan):
        tileno, tile_ss = np.divmod(slow_scan, cls.frag_ss_pixels)
        return tileno.astype(np.int16), tile_ss, fast_scan

    def to_distortion_array(self, allow_negative_xy=False):
        """Return distortion matrix for AGIPD detector, suitable for pyFAI.

        Parameters
        ----------

        allow_negative_xy: bool
          If False (default), shift the origin so no x or y coordinates are
          negative. If True, the origin is the detector centre.

        Returns
        -------
        out: ndarray
            Array of float 32 with shape (8192, 128, 4, 3).
            The dimensions mean:

            - 8192 = 16 modules * 512 pixels (slow scan axis)
            - 128 pixels (fast scan axis)
            - 4 corners of each pixel
            - 3 numbers for z, y, x
        """
        # Overridden only for docstring
        return super().to_distortion_array(allow_negative_xy)


class AGIPD_500K2GGeometry(DetectorGeometryBase):
    """Detector layout for AGIPD-500k

    The coordinates used in this class are 3D (x, y, z), and represent metres.

    You won't normally instantiate this class directly:
    use one of the constructor class methods to create or load a geometry.
    """
    detector_type_name = 'AGIPD-500K2G'
    pixel_size = 2e-4  # 2e-4 metres == 0.2 mm
    frag_ss_pixels = 64
    frag_fs_pixels = 128
    expected_data_shape = (8, 512, 128)
    n_modules = 8
    n_tiles_per_module = 8

    @classmethod
    def from_origin(cls, origin=(0, 0), asic_gap=2, panel_gap=(16, 30),
                    unit=pixel_size):
        """Generate an AGIPD-500K2G geometry from origin position.

        This produces an idealised geometry, assuming all modules are perfectly
        flat, aligned and equally spaced within the detector.

        The default origin (0, 0) of the coordinates is the bottom-right corner
        of the detector. If another coordinate is given as the origin, it is
        relative to the bottom-right corner. Coordinates increase upwards and
        to the left (looking along the beam).

        To give positions in units other than pixels, pass the *unit* parameter
        as the length of the unit in metres. E.g. ``unit=1e-3`` means the
        coordinates are in millimetres.
        """
        asic_gap_px = asic_gap * unit / cls.pixel_size
        panel_gap_x = panel_gap[0] * cls.pixel_size
        panel_gap_y = panel_gap[1] * cls.pixel_size

        # How much space one tile takes up, including the gaps
        # separating it from its neighbour.
        # In the y dimension, 128 px
        module_height = cls.frag_fs_pixels * cls.pixel_size
        # In x, 64 px + gap between tiles (asics)
        tile_width = (cls.frag_ss_pixels + asic_gap_px) * cls.pixel_size
        module_width = 8 * tile_width - asic_gap_px * cls.pixel_size  # 8 tiles + 7 gaps

        # coordinates relative to the first pixel of the first module
        # detector's bottom-right corner
        ref = (
            - module_width,  # x
            - (3 * (cls.frag_fs_pixels + panel_gap[1]) * unit)  # y
        )
        # origin
        ref = (- (origin[0] * unit + ref[0]), - (origin[1] * unit + ref[1]))

        modules = []
        for p in range(cls.n_modules):
            panel_corner_y = ref[1] - ((p // 2) * (module_height + panel_gap_y))
            panel_corner_x = ref[0] + ((p % 2) * (module_width + panel_gap_x))

            tiles = []
            modules.append(tiles)

            for a in range(cls.n_tiles_per_module):
                corner_x = panel_corner_x - a * tile_width

                tiles.append(GeometryFragment(
                    corner_pos=np.array([corner_x, panel_corner_y, 0.]),
                    ss_vec=np.array([-1, 0, 0]) * unit,
                    fs_vec=np.array([0, 1, 0]) * unit,
                    ss_pixels=cls.frag_ss_pixels,
                    fs_pixels=cls.frag_fs_pixels,
                ))

        return cls(modules)

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------

        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        # Label modules and tiles
        for ch, module in enumerate(self.modules):
            quad, mod = (1 + v for v in divmod(ch, 4))
            cx, cy, _ = module[4].centre() * scale
            ax.text(cx, cy, f'Q{quad}M{mod}', fontweight='bold',
                    verticalalignment='center',
                    horizontalalignment='center')

            for t in [0, 7]:
                cx, cy, _ = module[t].centre() * scale
                ax.text(cx, cy, 'T{}'.format(t + 1),
                        verticalalignment='center',
                        horizontalalignment='center')

        ax.set_title(f'AGIPD-500K2G detector geometry ({self.filename})')
        ax.set_aspect(1)
        return ax

    @staticmethod
    def split_tiles(module_data):
        # Split into 8 tiles along the slow-scan axis
        # This simple slicing is faster than np.split().
        return [module_data[..., s:s+64, :] for s in range(0, 512, 64)]

    @classmethod
    def _tile_slice(cls, tileno):
        # Which part of the array is this tile?
        # tileno = 0 to 7
        tile_offset = tileno * cls.frag_ss_pixels
        ss_slice = slice(tile_offset, tile_offset + cls.frag_ss_pixels)
        fs_slice = slice(0, cls.frag_fs_pixels)  # Every tile covers the full 128 pixels
        return ss_slice, fs_slice

    @classmethod
    def _module_coords_to_tile(cls, slow_scan, fast_scan):
        tileno, tile_ss = np.divmod(slow_scan, cls.frag_ss_pixels)
        return tileno.astype(np.int16), tile_ss, fast_scan

    def to_distortion_array(self, allow_negative_xy=False):
        """Return distortion matrix for AGIPD500K detector, suitable for pyFAI.

        Parameters
        ----------

        allow_negative_xy: bool
          If False (default), shift the origin so no x or y coordinates are
          negative. If True, the origin is the detector centre.

        Returns
        -------
        out: ndarray
            Array of float 32 with shape (4096, 128, 4, 3).
            The dimensions mean:

            - 8192 = 8 modules * 512 pixels (slow scan axis)
            - 128 pixels (fast scan axis)
            - 4 corners of each pixel
            - 3 numbers for z, y, x
        """
        # Overridden only for docstring
        return super().to_distortion_array(allow_negative_xy)

    def write_crystfel_geom(self, *args, **kwargs):
        super().write_crystfel_geom(*args, nquads=1, **kwargs)


def agipd_asic_seams():
    """Make a boolean array marking the double-width pixels in an AGIPD module

    This returns a (512, 128) array with False for normal (square) pixels, and
    True for the 400 x 200 µm pixels at the horizontal joins between ASICs.
    """
    arr = np.zeros((512, 128), dtype=np.bool_)
    # The outer edges (indexes 0 & 511) appear to be normal pixels
    arr[64::64] = True
    arr[63:511:64] = True
    return arr


class LPD_1MGeometry(DetectorGeometryBase):
    """Detector layout for LPD-1M

    The coordinates used in this class are 3D (x, y, z), and represent metres.

    You won't normally instantiate this class directly:
    use one of the constructor class methods to create or load a geometry.
    """
    detector_type_name = 'LPD-1M'
    pixel_size = 5e-4  # 5e-4 metres == 0.5 mm
    frag_ss_pixels = 32
    frag_fs_pixels = 128
    n_quads = 4
    n_modules = 16
    n_tiles_per_module = 16
    expected_data_shape = (16, 256, 256)
    _draw_first_px_on_tile = 8  # The first pixel in stored data is on tile 8

    @classmethod
    def from_quad_positions(cls, quad_pos, *, unit=1e-3, asic_gap=None,
                            panel_gap=None):
        """Generate an LPD-1M geometry from quadrant positions.

        This produces an idealised geometry, assuming all modules are perfectly
        flat, aligned and equally spaced within their quadrant.

        The quadrant positions refer to the corner of each quadrant
        where module 4, tile 16 is positioned.
        This is the corner of the last pixel as the data is stored.
        In the initial detector layout, the corner positions are for the top
        left corner of the quadrant, looking along the beam.

        The origin of the coordinates is in the centre of the detector.
        Coordinates increase upwards and to the left (looking along the beam).

        Parameters
        ----------
        quad_pos: list of 2-tuples
          (x, y) coordinates of the last corner (the one by module 4) of each
          quadrant.
        unit: float, optional
          The conversion factor to put the coordinates into metres.
          The default 1e-3 means the numbers are in millimetres.
        asic_gap: float, optional
          The gap between adjacent tiles/ASICs. The default is 4 pixels.
        panel_gap: float, optional
          The gap between adjacent modules/panels. The default is 4 pixels.
        """
        px_conversion = unit / cls.pixel_size
        asic_gap_px = 4 if (asic_gap is None) else asic_gap * px_conversion
        panel_gap_px = 4 if (panel_gap is None) else panel_gap * px_conversion

        # How much space one panel/module takes up, including the 'panel gap'
        # separating it from its neighbour.
        # In the x dimension, we have only one asic gap (down the centre)
        panel_width = (256 + asic_gap_px + panel_gap_px) * cls.pixel_size
        # In y, we have 7 gaps between the 8 ASICs in each column.
        panel_height = (256 + (7 * asic_gap_px) + panel_gap_px) * cls.pixel_size

        # How much space does one tile take up, including gaps to its neighbours?
        tile_width = (cls.frag_fs_pixels + asic_gap_px) * cls.pixel_size
        tile_height = (cls.frag_ss_pixels + asic_gap_px) * cls.pixel_size

        # Size of a tile from corner to corner, excluding gaps
        tile_size = np.array([cls.frag_fs_pixels, cls.frag_ss_pixels, 0]) * cls.pixel_size

        panels_across = [-1, -1, 0, 0]
        panels_up = [0, -1, -1, 0]
        modules = []
        for p in range(cls.n_modules):
            quad = p // 4
            quad_corner_x = quad_pos[quad][0] * unit
            quad_corner_y = quad_pos[quad][1] * unit

            p_in_quad = p % 4
            # Top beam-left corner of panel
            panel_corner_x = (quad_corner_x +
                              (panels_across[p_in_quad] * panel_width))
            panel_corner_y = (quad_corner_y +
                              (panels_up[p_in_quad] * panel_height))

            tiles = []
            modules.append(tiles)

            for a in range(cls.n_tiles_per_module):
                if a < 8:
                    up = -a
                    across = -1
                else:
                    up = -(15 - a)
                    across = 0

                tile_last_corner = (
                    np.array([panel_corner_x, panel_corner_y, 0.0])
                    + np.array([across, 0, 0]) * tile_width
                    + np.array([0, up, 0]) * tile_height
                )
                tile_first_corner = tile_last_corner - tile_size

                tiles.append(GeometryFragment(
                    corner_pos=tile_first_corner,
                    ss_vec=np.array([0, 1, 0]) * cls.pixel_size,
                    fs_vec=np.array([1, 0, 0]) * cls.pixel_size,
                    ss_pixels=cls.frag_ss_pixels,
                    fs_pixels=cls.frag_fs_pixels,
                ))
        return cls(modules)

    @classmethod
    def from_h5_file_and_quad_positions(cls, path, positions, unit=1e-3):
        """Load an LPD-1M geometry from an XFEL HDF5 format geometry file

        By default, both the quadrant positions and the positions
        in the file are measured in millimetres; the unit parameter controls
        this. The passed positions override quadrant positions from the file, if
        it contains them: see :meth:`from_h5_file` to use them.

        The origin of the coordinates is in the centre of the detector.
        Coordinates increase upwards and to the left (looking along the beam).

        This version of the code only handles x and y translation,
        as this is all that is recorded in the initial LPD geometry file.

        Parameters
        ----------

        path : str
          Path of an EuXFEL format (HDF5) geometry file for LPD.
        positions : list of 2-tuples
          (x, y) coordinates of the last corner (the one by module 4) of each
          quadrant.
        unit : float, optional
          The conversion factor to put the coordinates into metres.
          The default 1e-3 means the numbers are in millimetres.
        """
        assert len(positions) == 4
        modules = []
        with h5py.File(path, 'r') as f:
            for Q, M in product(range(1, 5), range(1, 5)):
                quad_pos = np.array(positions[Q - 1])
                mod_grp = f['Q{}/M{}'.format(Q, M)]
                mod_offset = mod_grp['Position'][:2]

                tiles = []
                for T in range(1, cls.n_tiles_per_module+1):
                    corner_pos = np.zeros(3)
                    tile_offset = mod_grp['T{:02}/Position'.format(T)][:2]
                    corner_pos[:2] = quad_pos + mod_offset + tile_offset

                    # Convert units (mm) to metres
                    corner_pos *= unit

                    # LPD geometry is measured to the last pixel of each tile.
                    # Subtract tile dimensions for the position of 1st pixel.
                    ss_vec = np.array([0, 1, 0]) * cls.pixel_size
                    fs_vec = np.array([1, 0, 0]) * cls.pixel_size
                    first_px_pos = (corner_pos
                                    - (ss_vec * cls.frag_ss_pixels)
                                    - (fs_vec * cls.frag_fs_pixels))

                    tiles.append(GeometryFragment(
                        corner_pos=first_px_pos,
                        ss_vec=ss_vec,
                        fs_vec=fs_vec,
                        ss_pixels=cls.frag_ss_pixels,
                        fs_pixels=cls.frag_fs_pixels,
                    ))
                modules.append(tiles)

        return cls(modules, filename=path)

    @classmethod
    def from_h5_file(cls, path):
        """Load an LPD-1M geometry from an XFEL HDF5 format geometry file

        This requires a file containing quadrant positions, which not all
        XFEL geometry files do. Use :meth:`from_h5_file_and_quad_positions` to
        load a file which does not have them.

        Parameters
        ----------

        path : str
          Path of an EuXFEL format (HDF5) geometry file for LPD.
        """
        with h5py.File(path, 'r') as f:
            try:
                quadpos = [f[f'Q{Q}/Position'][:2] for Q in range(1, 5)]
            except KeyError:
                raise ValueError(
                    "This HDF5 geometry file does not include quadrant positions. "
                    "You can use it with separately specified positions by "
                    "calling from_h5_file_and_quad_positions()"
                )

        return cls.from_h5_file_and_quad_positions(path, quadpos)

    def to_h5_file_and_quad_positions(self, path):
        """Write this geometry to an XFEL HDF5 format geometry file

        The quadrant positions are stored in the file, but also returned.
        These and the numbers in the file are in millimetres.

        The file and quadrant positions produced by this method are compatible
        with :meth:`from_h5_file_and_quad_positions`.
        """

        quad_pos = []

        for q in range(4):
            quad_fragmts_corners = []
            for mod in self.modules[q * 4: (q + 1) * 4]:
                quad_fragmts_corners.extend(f.corners() for f in mod)

            quad_points_xy = np.concatenate(quad_fragmts_corners)[:, :2]
            quad_pos.append(quad_points_xy.max(axis=0))

        quad_pos = np.stack(quad_pos)

        module_offsets = []
        tile_offsets = []

        for m, module in enumerate(self.modules):
            tile_positions = np.stack([f.corners().max(axis=0)[:2] for f in module])
            module_position = tile_positions.max(axis=0)
            tile_offsets.append(tile_positions - module_position)
            module_offsets.append(module_position - quad_pos[m // 4])

        with h5py.File(path, 'w') as hf:
            for q in range(4):
                Q = q + 1
                hf[f'Q{Q}/Position'] = quad_pos[q] * 1000  # m -> mm

            for m in range(16):
                Q, M = (m // 4) + 1, (m % 4) + 1
                mod_grp = hf.create_group(f'Q{Q}/M{M}')
                mod_grp['Position'] = module_offsets[m] * 1000  # m -> mm

                for t in range(self.n_tiles_per_module):
                    T = t + 1
                    mod_grp[f'T{T:02}/Position'] = tile_offsets[m][t] * 1000  # m -> mm

        return quad_pos * 1000  # m -> mm

    def quad_positions(self, h5_file=None):
        """Get the positions of the 4 quadrants

        Quadrant positions are returned as (x, y) coordinates in millimetres.
        Their meaning is as in :meth:`from_h5_file_and_quad_positions`.

        To use the returned positions with an existing XFEL HDF5 geometry file,
        the path to that file should be passed in. In that case, the offsets of
        M4 T16 in each quadrant are read from the file to calculate a suitable
        quadrant position. The geometry in the file is not checked against this
        geometry object at all.
        """
        positions = np.zeros((4, 2), dtype=np.float64)

        if h5_file is None:
            for q in range(4):
                quad_fragmts_corners = []
                for mod in self.modules[q * 4: (q + 1) * 4]:
                    quad_fragmts_corners.extend(f.corners() for f in mod)

                quad_points_xy = np.concatenate(quad_fragmts_corners)[:, :2]
                positions[q] = quad_points_xy.max(axis=0) * 1000  # m -> mm
        else:
            with h5py.File(h5_file, 'r') as f:
                for q in range(4):
                    # XFEL HDF5 geometry files for LPD record the position of
                    # the high-x, high-y corner of each tile. Instead of
                    # identifying which corner this is, we'll just take a
                    # maximum over all 4 corners.
                    # This assumes the tile is axis-aligned - for now, the XFEL
                    # geometry format has no way to express rotation anyway.
                    m4t16 = self.modules[(q * 4) + 3][15]
                    m4t16_max_corner = m4t16.corners().max(axis=0)

                    mod_grp = f[f'Q{q + 1}/M4']
                    mod_offset = mod_grp['Position'][:2]
                    tile_offset = mod_grp['T16/Position'][:2]

                    tile_pos = m4t16_max_corner[:2] * 1000  # m (xyz) -> mm (xy)
                    positions[q] = tile_pos - tile_offset - mod_offset

        return positions

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------

        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        # Label modules and tiles
        for ch, module in enumerate(self.modules):
            s = 'Q{Q}M{M}'.format(Q=(ch // 4) + 1, M=(ch % 4) + 1)
            cx, cy, _ = module[0].centre() * scale
            ax.text(cx, cy, s, fontweight='bold',
                    verticalalignment='center',
                    horizontalalignment='center')

            for t in [7, 8, 15]:
                cx, cy, _ = module[t].centre() * scale
                ax.text(cx, cy, 'T{}'.format(t + 1),
                        verticalalignment='center',
                        horizontalalignment='center')

        ax.set_title('LPD-1M detector geometry ({})'.format(self.filename))
        return ax

    @staticmethod
    def split_tiles(module_data):
        # This slicing is faster than using np.split()
        return [
            # Tiles 1-8 numbered top to bottom. Data starts at bottom, so
            # count backwards through them.
            module_data[..., y-32:y, :128] for y in range(256, 0, -32)
        ] + [
            # Tiles 9-16 numbered bottom to top.
            module_data[..., y:y+32, 128:] for y in range(0, 256, 32)
        ]

    @classmethod
    def _tile_slice(cls, tileno):
        # Which part of the array is this tile?
        if tileno < 8:  # First half of module (0 <= t <= 7)
            fs_slice = slice(0, 128)
            tiles_up = 7 - tileno
        else:  # Second half of module (8 <= t <= 15)
            fs_slice = slice(128, 256)
            tiles_up = tileno - 8
        tile_offset = tiles_up * 32
        ss_slice = slice(tile_offset, tile_offset + cls.frag_ss_pixels)
        return ss_slice, fs_slice

    @classmethod
    def _module_coords_to_tile(cls, slow_scan, fast_scan):
        tiles_across, tile_fs = np.divmod(fast_scan, cls.frag_fs_pixels)
        tiles_up, tile_ss = np.divmod(slow_scan, cls.frag_ss_pixels)

        # Each tiles_across is 0 or 1. To avoid iterating over the array with a
        # conditional, multiply the number we want by 1 and the other by 0.
        tileno = (
            (1 - tiles_across) * (7 - tiles_up)  # tileno 0-7
            + tiles_across * (tiles_up + 8)      # tileno 8-15
        )
        return tileno.astype(np.int16), tile_ss, tile_fs

    def to_distortion_array(self, allow_negative_xy=False):
        """Return distortion matrix for LPD detector, suitable for pyFAI.

        Parameters
        ----------

        allow_negative_xy: bool
          If False (default), shift the origin so no x or y coordinates are
          negative. If True, the origin is the detector centre.

        Returns
        -------
        out: ndarray
            Array of float 32 with shape (4096, 256, 4, 3).
            The dimensions mean:

            - 4096 = 16 modules * 256 pixels (slow scan axis)
            - 256 pixels (fast scan axis)
            - 4 corners of each pixel
            - 3 numbers for z, y, x
        """
        # Overridden only for docstring
        return super().to_distortion_array(allow_negative_xy)


def invert_xfel_lpd_geom(path_in, path_out):
    """Invert the coordinates in an XFEL geometry file (HDF5)

    The initial geometry file for LPD was recorded with the coordinates
    increasing down and to the right (looking in the beam direction), but the
    standard XFEL coordinate scheme is the opposite, increasing upwards and to
    the left (looking in beam direction).

    This utility function reads one file, and writes a second with the
    coordinates inverted.
    """
    with h5py.File(path_in, 'r') as fin, h5py.File(path_out, 'x') as fout:
        src_ds = fin['DetectorDescribtion']
        dst_ds = fout.create_dataset('DetectorDescription', data=src_ds)
        for k, v in src_ds.attrs.items():
            dst_ds.attrs[k] = v

        for Q, M in product(range(1, 5), range(1, 5)):
            path = 'Q{}/M{}/Position'.format(Q, M)
            fout[path] = -fin[path][:]
            for T in range(1, 17):
                path = 'Q{}/M{}/T{:02}/Position'.format(Q, M, T)
                fout[path] = -fin[path][:]


class DSSC_1MGeometry(DetectorGeometryBase):
    """Detector layout for DSSC-1M

    The coordinates used in this class are 3D (x, y, z), and represent metres.

    You won't normally instantiate this class directly:
    use one of the constructor class methods to create or load a geometry.
    """
    # Hexagonal pixels, 236 μm step in fast-scan axis, 204 μm in slow-scan
    detector_type_name = 'DSSC-1M'
    pixel_size = 236e-6
    frag_ss_pixels = 128
    frag_fs_pixels = 256
    n_quads = 4
    n_modules = 16
    n_tiles_per_module = 2
    expected_data_shape = (16, 128, 512)
    # This stretches the dimensions for the 'snapped' geometry so that its pixel
    # grid matches the aspect ratio of the detector pixels.
    _pixel_shape = np.array([1., 1.5/np.sqrt(3)], dtype=np.float64) * pixel_size

    # Pixel corners described clockwise from the top, assuming the reference
    # point for a pixel is outside it, aligned with the top point & left edge.
    # The unit is the width of a pixel, 236 μm.
    # The 4/3 extends the hexagons into the next row to correctly tessellate.
    _pixel_corners = np.stack([
        (np.array([0, 0.25, 0.75, 1, 0.75, 0.25]) * 4 / 3),
        [0.5, 1, 1, 0.5, 0, 0]
    ])

    @classmethod
    def from_quad_positions(cls, quad_pos, *, unit=1e-3, asic_gap=None,
                            panel_gap=None):
        """Generate a DSSC-1M geometry from quadrant positions.

        This produces an idealised geometry, assuming all modules are perfectly
        flat, aligned and equally spaced within their quadrant.

        The position given should refer to the bottom right (looking
        along the beam) corner of the quadrant.

        The origin of the coordinates is in the centre of the detector.
        Coordinates increase upwards and to the left (looking along the beam).

        Parameters
        ----------
        quad_pos: list of 2-tuples
          (x, y) coordinates of the last corner (the one by module 4) of each
          quadrant.
        unit: float, optional
          The conversion factor to put the coordinates into metres.
          The default 1e-3 means the numbers are in millimetres.
        asic_gap: float, optional
          The gap between adjacent tiles/ASICs. The default is 2 mm.
        panel_gap: float, optional
          The gap between adjacent modules/panels. The default is 4 mm.
        """
        assert len(quad_pos) == 4
        asic_gap_m = 2e-3 if (asic_gap is None) else asic_gap * unit
        panel_gap_m = 4e-3 if (panel_gap is None) else panel_gap * unit

        quads_x_orientation = [-1, -1, 1, 1]
        quads_y_orientation = [1, 1, -1, -1]

        frag_width = cls._pixel_shape[0] * cls.frag_fs_pixels
        frag_height = cls._pixel_shape[1] * cls.frag_ss_pixels
        module_width = (2 * frag_width) + asic_gap_m
        quad_height = (4 * frag_height) + (3 * panel_gap_m)

        module_step_vec = np.array([0, frag_height + panel_gap_m, 0])
        tile_step_vec = np.array([frag_width + asic_gap_m, 0, 0])

        modules = []

        for p in range(cls.n_modules):
            Q = p // 4
            x_orient = quads_x_orientation[Q]
            y_orient = quads_y_orientation[Q]
            quad_corner_x = quad_pos[Q][0] * unit
            quad_corner_y = quad_pos[Q][1] * unit

            p_in_quad = p % 4

            # Measuring in terms of the step within a row, the
            # step to the next row of hexagons is 1.5/sqrt(3).
            ss_vec = np.array([0, y_orient, 0]) * cls.pixel_size * 1.5 / np.sqrt(3)
            fs_vec = np.array([x_orient, 0, 0]) * cls.pixel_size

            # Corner position is measured at low-x, low-y corner (bottom
            # right as plotted). We want the position of the corner
            # with the first pixel, which is either high-x low-y or
            # low-x high-y.
            if x_orient == -1:
                quad_start_x = quad_corner_x + module_width
                quad_start_y = quad_corner_y
            else:  # y_orient == -1
                quad_start_x = quad_corner_x
                quad_start_y = quad_corner_y + quad_height

            quad_start = np.array([quad_start_x, quad_start_y, 0.])
            module_start = quad_start + (y_orient * p_in_quad * module_step_vec)

            modules.append([
                GeometryFragment(
                    corner_pos=module_start + (x_orient * t * tile_step_vec),
                    ss_vec=ss_vec,
                    fs_vec=fs_vec,
                    ss_pixels=cls.frag_ss_pixels,
                    fs_pixels=cls.frag_fs_pixels,
                ) for t in range(cls.n_tiles_per_module)
            ])

        return cls(modules)

    @classmethod
    def from_h5_file_and_quad_positions(cls, path, positions, unit=1e-3):
        """Load a DSSC geometry from an XFEL HDF5 format geometry file

        The position given should refer to the bottom right (looking
        along the beam) corner of the quadrant. The passed positions override
        quadrant positions from the file, if it contains them:
        see :meth:`from_h5_file` to use them.

        By default, both the quadrant positions and the positions
        in the file are measured in millimetres; the unit parameter controls
        this.

        The origin of the coordinates is in the centre of the detector.
        Coordinates increase upwards and to the left (looking along the beam).

        This version of the code only handles x and y translation,
        as this is all that is recorded in the initial LPD geometry file.

        Parameters
        ----------

        path : str
          Path of an EuXFEL format (HDF5) geometry file for DSSC.
        positions : list of 2-tuples
          (x, y) coordinates of the corner of each quadrant (the one with lowest
          x and y coordinates).
        unit : float, optional
          The conversion factor to put the coordinates into metres.
          The default 1e-3 means the numbers are in millimetres.
        """
        assert len(positions) == 4
        modules = []

        quads_x_orientation = [-1, -1, 1, 1]
        quads_y_orientation = [1, 1, -1, -1]

        with h5py.File(path, 'r') as f:
            for Q, M in product(range(1, 5), range(1, 5)):
                quad_pos = np.array(positions[Q - 1])
                mod_grp = f['Q{}/M{}'.format(Q, M)]
                mod_offset = mod_grp['Position'][:2]

                # Which way round is this quadrant
                x_orient = quads_x_orientation[Q - 1]
                y_orient = quads_y_orientation[Q - 1]

                tiles = []
                for T in range(1, 3):
                    corner_pos = np.zeros(3)
                    tile_offset = mod_grp['T{:02}/Position'.format(T)][:2]
                    corner_pos[:2] = quad_pos + mod_offset + tile_offset

                    # Convert units (mm) to metres
                    corner_pos *= unit

                    # Measuring in terms of the step within a row, the
                    # step to the next row of hexagons is 1.5/sqrt(3).
                    ss_vec = np.array([0, y_orient, 0]) * cls.pixel_size * 1.5/np.sqrt(3)
                    fs_vec = np.array([x_orient, 0, 0]) * cls.pixel_size

                    # Corner position is measured at low-x, low-y corner (bottom
                    # right as plotted). We want the position of the corner
                    # with the first pixel, which is either high-x low-y or
                    # low-x high-y.
                    if x_orient == -1:
                        first_px_pos = corner_pos - (fs_vec * cls.frag_fs_pixels)
                    else:
                        first_px_pos = corner_pos - (ss_vec * cls.frag_ss_pixels)

                    tiles.append(GeometryFragment(
                        corner_pos=first_px_pos,
                        ss_vec=ss_vec,
                        fs_vec=fs_vec,
                        ss_pixels=cls.frag_ss_pixels,
                        fs_pixels=cls.frag_fs_pixels,
                    ))
                modules.append(tiles)

        return cls(modules, filename=path)

    @classmethod
    def from_h5_file(cls, path):
        """Load a DSSC geometry from an XFEL HDF5 format geometry file

        This requires a file containing quadrant positions, which not all
        XFEL geometry files do. Use :meth:`from_h5_file_and_quad_positions` to
        load a file which does not have them.

        Parameters
        ----------

        path : str
          Path of an EuXFEL format (HDF5) geometry file for DSSC.
        """
        with h5py.File(path, 'r') as f:
            try:
                quadpos = [f[f'Q{Q}/Position'][:2] for Q in range(1, 5)]
            except KeyError:
                raise ValueError(
                    "This HDF5 geometry file does not include quadrant positions. "
                    "You can use it with separately specified positions by "
                    "calling from_h5_file_and_quad_positions()"
                )

        return cls.from_h5_file_and_quad_positions(path, quadpos)

    def to_h5_file_and_quad_positions(self, path):
        """Write this geometry to an XFEL HDF5 format geometry file

        The quadrant positions are stored in the file, but also returned.
        These and the numbers in the file are in millimetres.

        The file and quadrant positions produced by this method are compatible
        with :meth:`from_h5_file_and_quad_positions`.
        """

        quad_pos = []

        for q in range(4):
            quad_fragmts_corners = []
            for mod in self.modules[q * 4: (q + 1) * 4]:
                quad_fragmts_corners.extend(f.corners() for f in mod)

            quad_points_xy = np.concatenate(quad_fragmts_corners)[:, :2]
            quad_pos.append(quad_points_xy.min(axis=0))

        quad_pos = np.stack(quad_pos)

        module_offsets = []
        tile_offsets = []

        for m, module in enumerate(self.modules):
            tile_positions = np.stack([f.corners().min(axis=0)[:2] for f in module])
            module_position = tile_positions.min(axis=0)
            tile_offsets.append(tile_positions - module_position)
            module_offsets.append(module_position - quad_pos[m // 4])

        with h5py.File(path, 'w') as hf:
            for q in range(4):
                Q = q + 1
                hf[f'Q{Q}/Position'] = quad_pos[q] * 1000  # m -> mm

            for m in range(16):
                Q, M = (m // 4) + 1, (m % 4) + 1
                mod_grp = hf.create_group(f'Q{Q}/M{M}')
                mod_grp['Position'] = module_offsets[m] * 1000  # m -> mm

                for t in range(self.n_tiles_per_module):
                    T = t + 1
                    mod_grp[f'T{T:02}/Position'] = tile_offsets[m][t] * 1000  # m -> mm

        return quad_pos * 1000  # m -> mm

    def quad_positions(self, h5_file=None):
        """Get the positions of the 4 quadrants

        Quadrant positions are returned as (x, y) coordinates in millimetres.
        Their meaning is as in :meth:`from_h5_file_and_quad_positions`.

        To use the returned positions with an existing XFEL HDF5 geometry file,
        the path to that file should be passed in. In that case, the offsets of
        M1 T1 in each quadrant are read from the file to calculate a suitable
        quadrant position. The geometry in the file is not checked against this
        geometry object at all.
        """
        positions = np.zeros((4, 2), dtype=np.float64)

        if h5_file is None:
            for q in range(4):
                quad_fragmts_corners = []
                for mod in self.modules[q * 4: (q + 1) * 4]:
                    quad_fragmts_corners.extend(f.corners() for f in mod)

                quad_points_xy = np.concatenate(quad_fragmts_corners)[:, :2]
                positions[q] = quad_points_xy.min(axis=0) * 1000  # m -> mm
        else:
            with h5py.File(h5_file, 'r') as f:
                for q in range(4):
                    # XFEL HDF5 geometry files record the position of the low-x,
                    # low-y corner of each tile. Instead of identifying which
                    # corner this is, we'll just take a minimum over all 4 corners.
                    # This assumes the tile is axis-aligned - for now, the XFEL
                    # geometry format has no way to express rotation anyway.
                    m1t1_min_corner = self.modules[q * 4][0].corners().min(axis=0)

                    mod_grp = f[f'Q{q + 1}/M1']
                    mod_offset = mod_grp['Position'][:2]
                    tile_offset = mod_grp['T01/Position'][:2]

                    tile_pos = m1t1_min_corner[:2] * 1000  # m (xyz) -> mm (xy)
                    positions[q] = tile_pos - tile_offset - mod_offset

        return positions

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------

        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        # Label modules and tiles
        for ch, module in enumerate(self.modules):
            s = 'Q{Q}M{M}'.format(Q=(ch // 4) + 1, M=(ch % 4) + 1)
            cx, cy, _ = module[0].centre() * scale
            ax.text(cx, cy, s, fontweight='bold',
                    verticalalignment='center',
                    horizontalalignment='center')

            for t in [1]:
                cx, cy, _ = module[t].centre() * scale
                ax.text(cx, cy, 'T{}'.format(t + 1),
                        verticalalignment='center',
                        horizontalalignment='center')

        ax.set_title('DSSC detector geometry ({})'.format(self.filename))
        return ax

    @staticmethod
    def split_tiles(module_data):
        # Split into 2 tiles along the fast-scan axis
        # This simple slicing is faster than np.split().
        return [module_data[..., :256], module_data[..., 256:]]

    def plot_data_fast(self,
                       data, *,
                       axis_units='px',
                       frontview=True,
                       ax=None,
                       figsize=None,
                       colorbar=False,
                       **kwargs):

        ax = super().plot_data_fast(data,
                                    axis_units=axis_units,
                                    frontview=frontview,
                                    ax=ax,
                                    figsize=figsize,
                                    colorbar=colorbar,
                                    **kwargs)
        if axis_units == 'px':
            # Squash image to physically equal aspect ratio, so a circle projected
            # on the detector looks like a circle on screen.
            ax.set_aspect(204/236.)
        return ax

    @classmethod
    def _tile_slice(cls, tileno):
        tile_offset = tileno * cls.frag_fs_pixels
        fs_slice = slice(tile_offset, tile_offset + cls.frag_fs_pixels)
        ss_slice = slice(0, cls.frag_ss_pixels)  # Every tile covers the full pixel range
        return ss_slice, fs_slice

    def to_distortion_array(self, allow_negative_xy=False):
        """Return distortion matrix for DSSC detector, suitable for pyFAI.

        Parameters
        ----------

        allow_negative_xy: bool
          If False (default), shift the origin so no x or y coordinates are
          negative. If True, the origin is the detector centre.

        Returns
        -------
        out: ndarray
            Array of float 32 with shape (2048, 512, 6, 3).
            The dimensions mean:

            - 2048 = 16 modules * 128 pixels (slow scan axis)
            - 512 pixels (fast scan axis)
            - 6 corners of each pixel
            - 3 numbers for z, y, x
        """
        # Overridden only for docstring
        return super().to_distortion_array(allow_negative_xy=allow_negative_xy)

    @classmethod
    def _adjust_pixel_coords(cls, ss_coords, fs_coords, centre):
        # Shift odd-numbered rows by half a pixel.
        fs_coords[1::2] -= 0.5
        if centre:
            # Vertical (slow scan) centre is 2/3 of the way to the start of the
            # next row of hexagons, because the tessellating pixels extend
            # beyond the start of the next row.
            ss_coords += 2/3
            fs_coords += 0.5


class DSSC_Geometry(DSSC_1MGeometry):
    """DEPRECATED: Use DSSC_1MGeometry instead"""
    def __init__(self, modules, filename='No file', metadata=None):
        super().__init__(modules, filename, metadata)
        warnings.warn(
            "DSSC_Geometry has been renamed to DSSC_1MGeometry.", stacklevel=2
        )


class JUNGFRAUGeometry(DetectorGeometryBase):
    """Detector layout for flexible Jungfrau arrangements

     The base JUNGFRAU unit (and rigid group) in combined arrangements is the
     JF-500K module, which is an independent detector unit of 2 x 4 ASIC tiles.

     In the default orientation, the slow-scan dimension is y and the fast-scan
     dimension is x, so the data shape for one module is (y, x).
    """
    detector_type_name = 'JUNGFRAU'
    pixel_size = 7.5e-5   # 7.5e-5 metres = 75 micrometer = 0.075 mm
    frag_ss_pixels = 256  # pixels along slow scan axis within tile
    frag_fs_pixels = 256  # pixels along fast scan axis within tile
    expected_data_shape = (0, 512, 1024)  # num modules filled at instantiation
    n_tiles_per_module = 8

    def __init__(self, modules, filename='No file', metadata=None):
        super().__init__(modules, filename, metadata)
        self.expected_data_shape = (len(modules), 512, 1024)
        self.n_modules = len(modules)

    @classmethod
    def from_module_positions(cls,offsets=((0,0),), orientations=None,
                              asic_gap=2, unit=pixel_size):
        """Generate a Jungfrau geometry object from module positions

        Parameters
        ----------

        offsets: iterable of tuples
          iterable of length n_modules containing a coordinate tuple (x,y)
          for each offset to the global origin. Coordinates are in pixel units
          by default.

          These offsets are positions for the bottom, beam-right corner of each
          module, regardless of its orientation.

        orientations: iterable of tuples
          list of length n_modules containing a unit-vector tuple (x,y) for
          each orientation wrt. the axes

          Orientations default to (1,1) for each module if this optional
          keyword argument is lacking; if not, the number of elements must
          match the number of modules as per offsets

        asic_gap: float
          The gap between the 8 ASICs within each module. This is in pixel units
          by default.

        unit: float
          The unit for *offsets* and *asic_gap*, in metres. Defaults to the
          pixel size (75 um).
        """
        px_conversion = unit / cls.pixel_size
        # fill orientations with defaults to match number of offsets
        if orientations is None:
            orientations = [(1,1) for _ in range(len(offsets))]
        else:
            if len(offsets) != len(orientations):
                print("Offsets and orientations have different number!")
        asic_gap *= px_conversion
        module_width = (4 * cls.frag_fs_pixels) + (3 * asic_gap)
        module_height = (2 * cls.frag_ss_pixels) + asic_gap
        modules = []
        for orientation, offset in zip(orientations, offsets):
            x_orient, y_orient = orientation
            # Correct corner-offsets in case of flipped modules
            if x_orient == 1:
                x_offset = offset[0]
            else:
                x_offset = offset[0] + module_width
            if y_orient == 1:
                y_offset = offset[1]
            else:
                y_offset = offset[1] + module_height
            tiles = []
            for a in range(8):
                row = a // 4     # 0, 1
                column = a % 4   # 0, 1, 2, 3
                corner_y = (y_offset * px_conversion)\
                           + y_orient * (cls.frag_fs_pixels + asic_gap) * row
                corner_x = (x_offset * px_conversion)\
                           + x_orient * (cls.frag_ss_pixels + asic_gap) * column
                tiles.append(GeometryFragment(
                    corner_pos=np.array([corner_x, corner_y, 0.]) * cls.pixel_size,
                    fs_vec=np.array([x_orient, 0, 0]) * cls.pixel_size,
                    ss_vec=np.array([0, y_orient, 0]) * cls.pixel_size,
                    ss_pixels=cls.frag_ss_pixels,
                    fs_pixels=cls.frag_fs_pixels,
                ))
            modules.append(tiles)
        return cls(modules)

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------

        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        for m in range(len(self.modules)):
            tiles = self.modules[m]

            # Label tiles in the module: A0 to A8
            for t, tile in enumerate(tiles):
                s = 'M{M}A{T}'.format(T=t, M=m)
                cx, cy, _ = tile.centre() * scale
                ax.text(cx, cy, s, fontweight='bold',
                        verticalalignment='center',
                        horizontalalignment='center')

        ax.set_title('Jungfrau detector geometry ({})'.format(self.filename))
        print(' Expected data shape:', self.expected_data_shape)
        return ax

    @staticmethod
    def split_tiles(module_data):
        # 2 rows of 4 ASICs each. This slicing is faster than np.split().
        return [
            module_data[..., :256, x:x+256] for x in range(0, 1024, 256)
        ] + [
            module_data[..., 256:, x:x+256] for x in range(0, 1024, 256)
        ]

    @classmethod
    def _tile_slice(cls, tileno):
        # Which part of the array is this tile?
        # tileno = 0 to 7
        tile_ss_offset = (tileno // 4) * cls.frag_ss_pixels
        tile_fs_offset = (tileno % 4) * cls.frag_fs_pixels
        ss_slice = slice(tile_ss_offset, tile_ss_offset + cls.frag_ss_pixels)
        fs_slice = slice(tile_fs_offset, tile_fs_offset + cls.frag_fs_pixels)
        return ss_slice, fs_slice


class PNCCDGeometry(DetectorGeometryBase):
    """Detector layout for pnCCD

    The large-area, pn-junction Charge Coupled Device detector consists
    of two movable modules with a single tile each.

    In its default configuration, the complete detector frame is read
    out and written to file as a single image, with the modules split
    along the slow-scan dimension y. The public methods of this type
    support both the combined image array as well as separated module
    with :attr:`expected_data_shape`.
    """

    detector_type_name = 'PNCCD1MP'
    pixel_size = 75e-6
    frag_ss_pixels = 512
    frag_fs_pixels = 1024
    n_modules = 2
    n_tiles_per_module = 1
    expected_data_shape = (2, 512, 1024)

    # Each module has a rectangular cutout around the intended beam
    # position, i.e. at the bottom center of the top module (towards
    # negative y) and the top center of the bottom module (towards
    # positive y).
    cutout_width = 60 * pixel_size
    cutout_height = 22 * pixel_size

    @staticmethod
    def split_tiles(module_data):
        return [module_data]

    @classmethod
    def _tile_slice(cls, tileno):
        return np.s_[0:cls.frag_ss_pixels], np.s_[0:cls.frag_fs_pixels]

    @classmethod
    def _module_coords_to_tile(cls, slow_scan, fast_scan):
        return 0, slow_scan, fast_scan

    @classmethod
    def from_relative_positions(cls, gap=4e-3, top_offset=(0.0, 0.0, 0.0),
                                bottom_offset=(0.0, 0.0, 0.0)):
        """Generate a pnCCD geometry from relative module positions.

        The modules are assumed to be separated by the a gap centered
        around the beam (at the origin) in x, y and z = 0, with an
        optional offset applied to each module.

        Parameters
        ----------

        gap: float
          The gap between the detector modules centered around the beam,
          4mm (~50 px) by default.

        top_offset, bottom_offset: array_like of length 3
          Optional offset (x, y, z) for each module relative to the
          centered position.
        """

        top = np.array(top_offset) + np.array([
            -cls.frag_fs_pixels // 2 * cls.pixel_size,
            gap / 2 + cls.frag_ss_pixels * cls.pixel_size,
            0.0
        ])

        bottom = np.array(bottom_offset) + np.array([
            -cls.frag_fs_pixels // 2 * cls.pixel_size,
            -gap / 2,
            0.0
        ])

        return cls.from_absolute_positions(top, bottom)

    @classmethod
    def from_absolute_positions(cls, top, bottom):
        """Generate a pnCCD geometry from absolute module positions.

        Parameters
        ----------

        top, bottom: array_like of length 3
          Absolute position (x, y, z) for the first pixel of each
          module.
        """

        args = (np.array([0, -cls.pixel_size, 0]),
                np.array([cls.pixel_size, 0, 0]),
                cls.frag_ss_pixels, cls.frag_fs_pixels)

        return cls([[GeometryFragment(np.array(top), *args)],
                    [GeometryFragment(np.array(bottom), *args)]])

    @classmethod
    def _ensure_shape(cls, data):
        """Ensure image data has the proper shape.

        As a pnCCD frame is read out and saved as a single array, the
        public interface of this geometry implementation supports
        automatic reshaping to separated modules.
        """

        if data.shape[-2:] == (cls.n_modules * cls.frag_ss_pixels,
                               cls.frag_fs_pixels):
            data = data.reshape(*data.shape[:-2], cls.n_modules,
                                cls.frag_ss_pixels, cls.frag_fs_pixels)

        return data

    def inspect(self, axis_units='px', frontview=True):
        from matplotlib.patches import Rectangle

        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        # Add patches for each rectangular cutout.
        cutout_width = self.cutout_width * scale
        cutout_height = self.cutout_height * scale

        patch_dims = (cutout_width, cutout_height)
        patch_kwargs = dict(hatch='///', ec='red', fc='none')

        top, bottom = self.modules[0][0], self.modules[1][0]

        ax.add_patch(Rectangle(
            (top.centre()[0] * scale - cutout_width/2,
             top.corners()[2][1] * scale),
            *patch_dims, **patch_kwargs))

        ax.add_patch(Rectangle(
            (bottom.centre()[0] * scale - cutout_width/2,
             bottom.corners()[0][1] * scale - cutout_height),
            *patch_dims, **patch_kwargs))

        return ax

    def position_modules_fast(self, data, *args, **kwargs):
        return super().position_modules_fast(self._ensure_shape(data),
                                             *args, **kwargs)

    def plot_data_fast(self, data, *args, **kwargs):
        return super().plot_data_fast(self._ensure_shape(data),
                                      *args, **kwargs)


class EpixGeometryBase(DetectorGeometryBase):
    """Base class for ePix detector geometry. Subclassed for specific detectors.

    The first pixel of the first tile (ASIC 0) corresponds to the first pixel
    in the data `frame [0,0]`, and the first row of this tile corresponds to
    half of the row line in the data `frame[0,:ncol//2]`.

    The four tiles are stacked in 2 rows of 2 columns and numbered along with
    the columns:
    [ASIC 0] [ASIC 1]
    [ASIC 2] [ASIC 3]
    """
    # Geometry layout
    # ---------------
    # oooo-  -oooo
    # oooo-  -oooo
    # ||||#  #||||
    #
    # ||||#  #||||
    # oooo-  -oooo
    # oooo-  -oooo
    #
    # o : normal pixels 50x50 um
    # - : long horizontal pixels 50x175 um
    # | : long vertical pixels 175x50 um
    # # : big square pixels 175x175 um (4 big pixels in the center)
    #
    # pixel sizes given for ePix100, ePix10K has doubled size pixels
    #
    # Electronic layout
    # -----------------
    # detector: 704+4 rows, 768 columns, 4 asics
    # asic: 352+2 rows, 384 columns, 4 banks
    # bank: 352+2 rows, 96 columns
    #     - fast parallel column readout
    #     - sigma-delta ADC per column
    #     - single high speed LVDS link
    #
    # (+2) Each asics has two calibration rows:
    #     - pixel max (the first and last rows in the pixel array)
    #     - baseline (next rows)
    #
    # Schema
    # ------
    #      ||||||  ||||||
    #   =                 =
    #   =    A2      A1   =  352+2
    #
    #   =    A3      A0   =  352+2
    #   =                 =
    #      ||||||  ||||||
    #        384     384
    #
    # Conversions in electronic layout
    # -------------------------------
    # - ASIC position to No.:
    #     tileno = 2*(1 - row*col) + row - col
    #
    # - ASIC No. to position:
    #     hi, lo = tileno // 2, tileno % 2
    #     i, j = (hi - lo + 1) & 1, 1 - hi
    #
    # - ASIC slice:
    #     [(nrow-1)*row : nrow//2-row : 1-2*row] or [-1 : nrow//2-row : 1-2*row]
    #     [(ncol-1)*col : ncol//2-col : 1-2*col] or [-1 : ncol//2-col : 1-2*col]
    #
    # See also
    # --------
    # 1. [EuXFEL ePix documentation]
    #     (https://rtd.xfel.eu/docs/epix-documentation/en/latest/index.html)
    # 2. [A Dragone et al 2014 J. Phys.: Conf. Ser. 493 012012]
    #     (https://doi.org/10.1088/1742-6596/493/1/012012)
    # 3. [SLAC Confluence](https://confluence.slac.stanford.edu/display/PSDM/EPIX10KA)
    #     Describes ePix10KA instead of ePix100, many details are common for these
    #     two versions
    n_modules = 1
    n_tiles_per_module = 4
    fs_tiles = 2
    ss_tiles = 2

    @classmethod
    def from_origin(cls, origin=(0, 0), asic_gap=None, unit=None):
        """Generate an ePix100 geometry from origin position.

        This produces an idealised geometry, assuming all modules are perfectly
        flat, aligned and equally spaced within the detector.

        The default origin (0, 0) of the coordinates is the center of the
        detector. If another coordinate is given as the origin, it is relative
        to the center. Coordinates increase upwards and to the left (looking
        along the beam).

        To give positions in units other than pixels, pass the *unit* parameter
        as the length of the unit in metres. E.g. ``unit=1e-3`` means the
        coordinates are in millimetres.

        .. note::

            ePix100 has 2 different geometry layout:

            - A single monolithic sensor with a 2x2 array of four ASICs bonded to it.
              These would have no dead gaps but would have large pixels in the central
              cross. This is the current default gap implementation.
            - A pair of sensors with each sensor being bonded to two ASICs.
              These would have a dead gap equal to twice the guard ring width (~450-500um)
              plus a mechanical gap of about 200-300 microns. This would result in a total
              dead gap of about 1.25 millimeters.
              For this case see :meth:`from_relative_positions`
        """
        if unit is None:
            unit = cls.pixel_size
        if asic_gap is None:
            asic_gap = cls.asic_gap

        x0, y0 = origin[0] * unit, origin[1] * unit
        tiles = []
        gap = asic_gap * unit
        row_sz = cls.frag_ss_pixels * cls.pixel_size
        col_sz = cls.frag_fs_pixels * cls.pixel_size
        for tileno in range(4):
            row, col = tileno // 2, tileno % 2
            tiles.append(GeometryFragment(
                corner_pos=np.array(
                    [col_sz - col * (col_sz + gap) - x0 + gap / 2,
                     row_sz - row * (row_sz + gap) - y0 + gap / 2,
                     0]),
                ss_vec=np.array([0, -1, 0]) * cls.pixel_size,
                fs_vec=np.array([-1, 0, 0]) * cls.pixel_size,
                ss_pixels=cls.frag_ss_pixels,
                fs_pixels=cls.frag_fs_pixels,
            ))
        return cls([tiles])

    @classmethod
    def _module_coords_to_tile(cls, slow_scan, fast_scan):
        nrow, ncol = cls.frag_ss_pixels, cls.frag_fs_pixels
        row, tile_ss = np.divmod(slow_scan, nrow)
        col, tile_fs = np.divmod(fast_scan, ncol)

        tileno = 2 * row + col
        return tileno.astype(np.int16), tile_ss, tile_fs

    @classmethod
    def _tile_slice(cls, tileno):
        # Which part of the array is this tile?
        # tileno = 0 to 3
        row, col = tileno // 2, tileno % 2
        ss_slice = slice(cls.frag_ss_pixels * row,
                         cls.frag_ss_pixels * (row + 1))
        fs_slice = slice(cls.frag_fs_pixels * col,
                         cls.frag_fs_pixels * (col + 1))
        return ss_slice, fs_slice

    @classmethod
    def split_tiles(cls, module_data):
        # Split into 4 tiles
        return [module_data[
                ...,
                cls.frag_ss_pixels * row:cls.frag_ss_pixels * (row + 1),
                cls.frag_fs_pixels * col:cls.frag_fs_pixels * (col + 1)]
                for row, col in ((i // 2, i % 2) for i in range(4))]

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Axes object.

        Parameters
        ----------
        axis_units : str
            Show the detector scale in pixels ('px') or metres ('m').

        frontview : bool
            If True (the default), x increases to the left, as if you were
            looking along the beam. False gives a 'looking into the beam' view.
        """
        ax = super().inspect(axis_units=axis_units, frontview=frontview)
        scale = self._get_plot_scale_factor(axis_units)

        # Label modules and tiles
        module = self.modules[0]
        for t in range(4):
            cx, cy, _ = module[t].centre() * scale
            ax.text(cx, cy, f'ASIC {t}', fontweight='bold',
                    verticalalignment='center',
                    horizontalalignment='center')

        ax.set_title('{} detector geometry ({})'.format(
            self.detector_type_name, self.filename))
        return ax

    @classmethod
    def asic_seams(cls):
        """Make a boolean array marking the wide pixels

        This returns a full frame array with False for normal pixels, and
        True for the wide pixels at inner edges of ASICs.
        """
        npx_ss = cls.frag_ss_pixels
        npx_fs = cls.frag_fs_pixels
        ss_wides = np.full(cls.ss_tiles * npx_ss, False)
        ss_wides[npx_ss - 1:npx_ss + 1] = True
        fs_wides = np.full(cls.fs_tiles * npx_fs, False)
        fs_wides[npx_fs - 1:npx_fs + 1] = True
        return ss_wides[:, None] + fs_wides[None, :]

    @classmethod
    def pixel_areas(cls):
        """Make an array of pixel areas

        This returns a full frame array with pixel areas. Pixels on inner
        edges of ASICs are bigger.
        """
        npx_ss = cls.frag_ss_pixels
        npx_fs = cls.frag_fs_pixels
        ss_sizes = np.full(cls.ss_tiles * npx_ss, cls.pixel_size)
        ss_sizes[npx_ss - 1:npx_ss + 1] = cls.inner_pixel_size
        fs_sizes = np.full(cls.fs_tiles * npx_fs, cls.pixel_size)
        fs_sizes[npx_fs - 1:npx_fs + 1] = cls.inner_pixel_size
        return np.outer(ss_sizes, fs_sizes)

    @classmethod
    def normalize_data(cls, data):
        """Remove diagnostic pixels from the data

        EuXFEL ePix data can contain extra rows with diagnostic information. this method
        remove these row if they are present.
        """
        if data.shape[-2:] == (2 * cls.frag_ss_pixels + 4, 2 * cls.frag_fs_pixels):
            # EuXFEL data stored extra 2 row per ASIC used for diagnostic
            #   - pixel max (the first and last rows in the pixel array)
            #   - baseline (next rows)
            # we removed them here as they do not contain data
            data = data[..., 2:-2, :]
        return data

    @classmethod
    def _ensure_shape(cls, data):
        """Ensure image data has the proper shape.

        As a ePix frame is read out and saved as a single array, the
        public interface of this geometry implementation supports
        automatic reshaping to adding the modules dimension.
        """
        data = cls.normalize_data(data)

        # add module dimension (ePix data is stored without module dim)
        if data.ndim == 2:
            data = data[None, ...]
        elif data.ndim >= 3 and data.shape[-3] != 1:
            data = data[..., None, :, :]

        return data

    def position_modules_fast(self, data, *args, **kwargs):
        return super().position_modules_fast(self._ensure_shape(data),
                                             *args, **kwargs)

    def plot_data_fast(self, data, *args, **kwargs):
        return super().plot_data_fast(self._ensure_shape(data),
                                      *args, **kwargs)


class Epix100Geometry(EpixGeometryBase):
    """Detector layout for ePix100

    ePix100 detectors have one module, which is built from 4 ASICs
    with wide pixes on inner edges.

    In its default configuration, the complete detector frame is read
    out and written to file as a single image, with the ASICs split
    along the both dimensions.

    There are 4 more rows in raw data. These are calibration pixels.
    They have the same electronics as normal pixels but aren't wired
    to the sensor. They are two first and two last rows in the raw data
    array. This class assumes that calibration rows are cut.
    """
    detector_type_name = 'ePix100'
    pixel_size = 50e-6
    inner_pixel_size = 175e-6
    asic_gap = 2 * (inner_pixel_size - pixel_size) / pixel_size
    frag_ss_pixels = 352  # rows
    frag_fs_pixels = 384  # columns
    expected_data_shape = (
        EpixGeometryBase.n_modules,
        2 * frag_ss_pixels,
        2 * frag_fs_pixels
    )

    @classmethod
    def from_relative_positions(cls, asic_gap=None, unit=None, top=(0., 0., 0.),
                                bottom=(0., 0., 0.)):
        """Generate an ePix100 geometry from relative Asics-pair positions.

        ePix100 has 2 assemblies:

        - a single monolithic sensor with a 2x2 array of four ASICs bonded to it. These
          would have no dead gaps but would have large pixels in the central cross.
          Use :meth:`from_origin` if your detector has this layout.
        - A pair of sensors with each sensor being bonded to two ASICs. These would have
          a dead gap equal to twice the guard ring width (~450-500um) plus a mechanical
          gap of about 200-300 microns. This would result in a total dead gap of about
          1.25 millimeters.

        For the later case, one can determine determine the exact gap existing between
        the 2 (top and bottom) asic pair. A rough estimation of the gap has been seen at
        ~25 pixels. This can be generated with::

            geom = Epix100Geometry.from_relative_positions(
                top=[386.5, 364.5, 0.], bottom=[386.5, -12.5, 0.]
            )

        Parameters
        ----------

        asic_gap: float
            The gap between asics within a pair (default 250um)
        unit: float
            To give positions in units other than pixels, pass the *unit* parameter as
            the length of the unit in metres. E.g. ``unit=1e-3`` means the coordinates
            are in millimetres.
        top, bottom: array_like of length 3 Optional
            offset (x, y, z) for asic pair relative to the centered position.
        """
        unit = unit or cls.pixel_size
        geom = cls.from_origin(asic_gap=asic_gap)
        ref_top = geom.modules[0][0].corner_pos  # asic 0
        ref_bot = geom.modules[0][2].corner_pos  # asic 3

        top, bottom = np.array(top), np.array(bottom)
        position = [top, top, bottom, bottom]
        reference = [ref_top, ref_top, ref_bot, ref_bot]

        return cls([[
            tile.offset(pos * unit - ref)
            for tile, pos, ref in zip(geom.modules[0], position, reference)
        ]])


class Epix10KGeometry(EpixGeometryBase):
    """Detector layout for ePix10K

    ePix10K detectors have one module, which is built from 4 ASICs
    with wide pixes on inner edges.

    In its default configuration, the complete detector frame is read
    out and written to file as a single image, with the ASICs split
    along the both dimensions.

    There are 4 more rows in raw data. These are calibration pixels.
    They have the same electronics as normal pixels but aren't wired
    to the sensor. They are two first and two last rows in the raw data
    array. This class assumes that calibration rows are cut.
    """
    detector_type_name = 'ePix10K'
    pixel_size = 100e-6
    inner_pixel_size = 250e-6
    asic_gap = 2 * (inner_pixel_size - pixel_size) / pixel_size
    frag_ss_pixels = 176  # rows
    frag_fs_pixels = 192  # columns
    expected_data_shape = (
        EpixGeometryBase.n_modules,
        2 * frag_ss_pixels,
        2 * frag_fs_pixels
    )
