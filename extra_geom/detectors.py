"""AGIPD & LPD geometry handling."""
from cfelpyutils.crystfel_utils import load_crystfel_geometry
import h5py
from itertools import chain, product
import numpy as np
from scipy.ndimage import affine_transform
import warnings

from .crystfel_fmt import write_crystfel_geom
from .snapped import GridGeometryFragment, SnappedGeometry

__all__ = ['AGIPD_1MGeometry', 'LPD_1MGeometry']


class GeometryFragment:
    """Holds the 3D position & orientation of one detector tile

    corner_pos refers to the corner of the detector tile where the first pixel
    stored is located. The tile is assumed to be a rectangle of ss_pixels in
    the slow scan dimension and fs_pixels in the fast scan dimension.
    ss_vec and fs_vec are vectors for a step of one pixel in each dimension.

    The coordinates in this class are (x, y, z), in metres.
    """

    def __init__(self, corner_pos, ss_vec, fs_vec, ss_pixels, fs_pixels):
        self.corner_pos = corner_pos
        self.ss_vec = ss_vec
        self.fs_vec = fs_vec
        self.ss_pixels = ss_pixels
        self.fs_pixels = fs_pixels

    @classmethod
    def from_panel_dict(cls, d):
        res = d['res']
        corner_pos = np.array([d['cnx'], d['cny'], d['coffset']]) / res
        ss_vec = np.array([d['ssx'], d['ssy'], d['ssz']]) / res
        fs_vec = np.array([d['fsx'], d['fsy'], d['fsz']]) / res
        ss_pixels = d['max_ss'] - d['min_ss'] + 1
        fs_pixels = d['max_fs'] - d['min_fs'] + 1
        return cls(corner_pos, ss_vec, fs_vec, ss_pixels, fs_pixels)

    def corners(self):
        return np.stack([
            self.corner_pos,
            self.corner_pos + (self.fs_vec * self.fs_pixels),
            self.corner_pos + (self.ss_vec * self.ss_pixels) + (self.fs_vec * self.fs_pixels),
            self.corner_pos + (self.ss_vec * self.ss_pixels),
        ])

    def centre(self):
        return (
            self.corner_pos
            + (0.5 * self.ss_vec * self.ss_pixels)
            + (0.5 * self.fs_vec * self.fs_pixels)
        )

    def offset(self, shift):
        pos = self.corner_pos + shift
        return type(self)(pos, self.ss_vec, self.fs_vec, self.ss_pixels, self.fs_pixels)

    def snap(self, px_shape):
        # Round positions and vectors to integers, drop z dimension
        corner_pos = np.around(self.corner_pos[:2] / px_shape).astype(np.int32)
        ss_vec = np.around(self.ss_vec[:2] / px_shape).astype(np.int32)
        fs_vec = np.around(self.fs_vec[:2] / px_shape).astype(np.int32)

        # We should have one vector in the x direction and one in y, but
        # we don't know which is which.
        assert {tuple(np.abs(ss_vec)), tuple(np.abs(fs_vec))} == {(0, 1), (1, 0)}

        # Convert xy coordinates to yx indexes
        return GridGeometryFragment(
            corner_pos[::-1], ss_vec[::-1], fs_vec[::-1], self.ss_pixels, self.fs_pixels
        )


class DetectorGeometryBase:
    """Base class for detector geometry. Subclassed for specific detectors."""
    # Define in subclasses:
    detector_type_name = ''
    pixel_size = 0.0
    frag_ss_pixels = 0
    frag_fs_pixels = 0
    n_modules = 0
    n_tiles_per_module = 0
    expected_data_shape = (0, 0, 0)
    _pixel_corners = np.array([  # pixel units; overridden for DSSC
        [0, 1, 1, 0],  # slow-scan
        [0, 0, 1, 1]   # fast-scan
    ])
    _draw_first_px_on_tile = 1  # Tile num of 1st pixel - overridden for LPD

    @property
    def _pixel_shape(self):
        """Pixel (x, y) shape. Overridden for DSSC."""
        return np.array([1., 1.], dtype=np.float64) * self.pixel_size

    def __init__(self, modules, filename='No file'):
        # List of lists (1 per module) of fragments (1 per tile)
        self.modules = modules
        # self.filename is metadata for plots, we don't read/write the file.
        # There are separate methods for reading and writing.
        self.filename = filename
        self._snapped_cache = None

    def _get_plot_scale_factor(self, axis_units):
        if axis_units == 'm':
            return 1
        elif axis_units == 'px':
            return 1 / self.pixel_size
        else:
            raise ValueError("axis_units must be 'px' or 'm', not {!r}"
                             .format(axis_units))

    def inspect(self, axis_units='px', frontview=True):
        """Plot the 2D layout of this detector geometry.

        Returns a matplotlib Figure object.
        """
        import matplotlib.pyplot as plt
        from matplotlib.collections import PatchCollection, LineCollection
        from matplotlib.patches import Polygon

        scale = self._get_plot_scale_factor(axis_units)

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)

        rects = []
        first_rows = []
        for module in self.modules:
            for t, fragment in enumerate(module, start=1):
                corners = fragment.corners()[:, :2]  # Drop the Z dimension
                rects.append(Polygon(corners * scale))

                if t == self._draw_first_px_on_tile:
                    # Find the ends of the first row in reading order
                    c1 = fragment.corner_pos * scale
                    c2 = c1 + (fragment.fs_vec * fragment.fs_pixels * scale)
                    first_rows.append((c1[:2], c2[:2]))

        # Add tile shapes
        pc = PatchCollection(rects, facecolor=(0.75, 1.0, 0.75), edgecolor=None)
        ax.add_collection(pc)

        # Add markers for first pixels & lines for first row
        first_rows = np.array(first_rows)
        first_px_x, first_px_y = first_rows[:, 0, 0], first_rows[:, 0, 1]

        ax.scatter(first_px_x, first_px_y, marker='x', label='First pixel')
        ax.add_collection(LineCollection(
            first_rows, linestyles=':', color='k', label='First row'
        ))
        ax.legend()

        cross_size = 0.02 * scale

        # Draw cross in the centre.
        ax.hlines(0, -cross_size, +cross_size, colors='0.75', linewidths=2)
        ax.vlines(0, -cross_size, +cross_size, colors='0.75', linewidths=2)

        if frontview:
            ax.invert_xaxis()

        ax.set_xlabel('metres' if axis_units == 'm' else 'pixels')
        ax.set_ylabel('metres' if axis_units == 'm' else 'pixels')

        return ax

    def compare(self, other, scale=1.0):
        """Show a comparison of this geometry with another in a 2D plot.

        This shows the current geometry like :meth:`inspect`, with the addition
        of arrows showing how each panel is shifted in the other geometry.

        Parameters
        ----------

        other : DetectorGeometryBase
          A second geometry object to compare with this one.
          It should be for the same kind of detector.
        scale : float
          Scale the arrows showing the difference in positions.
          This is useful to show small differences clearly.
        """
        from matplotlib.collections import PatchCollection
        from matplotlib.patches import FancyArrow

        coord_scale = 1 / self.pixel_size
        arrow_scale = scale * coord_scale

        # Draw this geometry first, using pixel units
        ax = self.inspect()

        if len(self.modules) != len(other.modules):
            print("Geometry objects have different numbers of modules!")
        if any(len(mod_a) != len(mod_b) for (mod_a, mod_b) in zip(self.modules, other.modules)):
            print("Geometry objects have different numbers of fragments in a module!")

        arrows = []
        for mod_a, mod_b in zip(self.modules, other.modules):
            for frag_a, frag_b in zip(mod_a, mod_b):
                corners_a = frag_a.corners()[:, :2]  # Drop the Z dimension
                corner_a, corner_a_opp = corners_a[0], corners_a[2]

                corners_b = frag_b.corners()[:, :2]
                corner_b, corner_b_opp = corners_b[0], corners_b[2]

                # Arrow for first corner
                dx, dy = (corner_b - corner_a) * arrow_scale
                if not (dx == dy == 0):
                    sx, sy = corner_a * coord_scale
                    arrows.append(FancyArrow(
                        sx, sy, dx, dy, width=5, head_length=4
                    ))

                # Arrow for third corner
                dx, dy = (corner_b_opp - corner_a_opp) * arrow_scale
                if not (dx == dy == 0):
                    sx, sy = corner_a_opp * coord_scale
                    arrows.append(FancyArrow(
                        sx, sy, dx, dy, width=5, head_length=4
                    ))

        ac = PatchCollection(arrows)
        ax.add_collection(ac)

        ax.set_title('Geometry comparison: {} → {}'
                     .format(self.filename, other.filename))
        ax.text(1, 0, 'Arrows scaled: {}×'.format(scale),
                horizontalalignment="right", verticalalignment="bottom",
                transform=ax.transAxes)
        return ax

    @classmethod
    def from_crystfel_geom(cls, filename):
        """Read a CrystFEL format (.geom) geometry file.

        Returns a new geometry object.
        """
        geom_dict = load_crystfel_geometry(filename)
        modules = []
        for p in range(cls.n_modules):
            tiles = []
            modules.append(tiles)
            for a in range(cls.n_tiles_per_module):
                d = geom_dict['panels']['p{}a{}'.format(p, a)]
                tiles.append(GeometryFragment.from_panel_dict(d))
        return cls(modules, filename=filename)

    def write_crystfel_geom(self, filename, *,
                            data_path='/entry_1/instrument_1/detector_1/data',
                            mask_path=None, dims=('frame', 'modno', 'ss', 'fs'),
                            nquads=4, adu_per_ev=None, clen=None,
                            photon_energy=None):
        """Write this geometry to a CrystFEL format (.geom) geometry file.

        Parameters
        ----------

        filename : str
            Filename of the geometry file to write.
        data_path : str
            Path to the group that contains the data array in the hdf5 file.
            Default: ``'/entry_1/instrument_1/detector_1/data'``.
        mask_path : str
            Path to the group that contains the mask array in the hdf5 file.
        dims : tuple
            Dimensions of the data. Extra dimensions, except for the defaults,
            should be added by their index, e.g.
            ('frame', 'modno', 0, 'ss', 'fs') for raw data.
            Default: ``('frame', 'modno', 'ss', 'fs')``.
            Note: the dimensions must contain frame, ss, fs.
        adu_per_ev : float
            ADU (analog digital units) per electron volt for the considered
            detector.
        clen : float
            Distance between sample and detector in meters
        photon_energy : float
            Beam wave length in eV
        """
        write_crystfel_geom(
            self, filename, data_path=data_path, mask_path=mask_path, dims=dims,
            nquads=nquads, adu_per_ev=adu_per_ev, clen=clen,
            photon_energy=photon_energy,
        )

        if self.filename == 'No file':
            self.filename = filename

    def _snapped(self):
        """Snap geometry to a 2D pixel grid

        This returns a new geometry object. The 'snapped' geometry is
        less accurate, but can assemble data into a 2D array more efficiently,
        because it doesn't do any interpolation.
        """
        if self._snapped_cache is None:
            modules = []
            for module in self.modules:
                tiles = [t.snap(px_shape=self._pixel_shape) for t in module]
                modules.append(tiles)
            centre = -np.min([t.corner_idx for t in chain(*modules)], axis=0)

            # Offset by centre to make all coordinates >= 0
            modules = [
                [t.offset(centre) for t in module]
                for module in modules
            ]
            self._snapped_cache = SnappedGeometry(modules, self, centre)
        return self._snapped_cache

    @staticmethod
    def split_tiles(module_data):
        """Split data from a detector module into tiles.

        Must be implemented in subclasses.
        """
        raise NotImplementedError

    def output_array_for_position_fast(self, extra_shape=(), dtype=np.float32):
        """Make an empty output array to use with position_modules_fast

        You can speed up assembling images by reusing the same output array:
        call this once, and then pass the array as the ``out=`` parameter to
        :meth:`position_modules_fast()`. By default, it allocates a new array on
        each call, which can be slow.

        Parameters
        ----------

        extra_shape : tuple, optional
          By default, a 2D output array is generated, to assemble a single
          detector image. If you are assembling multiple pulses at once, pass
          ``extra_shape=(nframes,)`` to get a 3D output array.
        dtype : optional (Default: np.float32)
        """
        return self._snapped().make_output_array(extra_shape=extra_shape,
                                                 dtype=dtype)

    def position_modules_fast(self, data, out=None, threadpool=None):
        """Assemble data from this detector according to where the pixels are.

        This approximates the geometry to align all pixels to a 2D grid.

        Parameters
        ----------

        data : ndarray
          The last three dimensions should match the modules, then the
          slow scan and fast scan pixel dimensions.
        out : ndarray, optional
          An output array to assemble the image into. By default, a new
          array is allocated. Use :meth:`output_array_for_position_fast` to
          create a suitable array.
          If an array is passed in, it must match the dtype of the data and the
          shape of the array that would have been allocated.
          Parts of the array not covered by detector tiles are not overwritten.
          In general, you can reuse an output array if you are assembling
          similar pulses or pulse trains with the same geometry.
        threadpool : concurrent.futures.ThreadPoolExecutor, optional
          If passed, parallelise copying data into the output image.
          By default, data for different tiles are copied serially.
          For a single 1 MPx image, the default appears to be faster, but for
          assembling a stack of several images at once, multithreading can help.

        Returns
        -------
        out : ndarray
          Array with one dimension fewer than the input.
          The last two dimensions represent pixel y and x in the detector space.
        centre : ndarray
          (y, x) pixel location of the detector centre in this geometry.
        """
        return self._snapped().position_modules(data, out=out, threadpool=threadpool)

    def position_all_modules(self, data, out=None):
        """Deprecated alias for :meth:`position_modules_fast`"""
        return self.position_modules_fast(data, out=out)

    def position_modules_symmetric(self, data, out=None, threadpool=None):
        """Assemble data with the centre in the middle of the output array.

        The assembly process is the same as :meth:`position_modules_fast`,
        aligning each module to a single pixel grid. But this makes the output
        array symmetric, with the centre at (height // 2, width // 2).

        Parameters
        ----------

        data : ndarray
          The last three dimensions should match the modules, then the
          slow scan and fast scan pixel dimensions.
        out : ndarray, optional
          An output array to assemble the image into. By default, a new
          array is created at the minimum size to allow symmetric assembly.
          If an array is passed in, its last two dimensions must be at least
          this size.
        threadpool : concurrent.futures.ThreadPoolExecutor, optional
          If passed, parallelise copying data into the output image.
          See :meth:`position_modules_fast` for details.

        Returns
        -------
        out : ndarray
          Array with one dimension fewer than the input.
          The last two dimensions represent pixel y and x in the detector space.
        """
        return self._snapped().position_modules_symmetric(
            data, out=out, threadpool=threadpool
        )

    def plot_data_fast(self,
                       data, *,
                       axis_units='px',
                       frontview=True,
                       ax=None,
                       figsize=None,
                       colorbar=True,
                       **kwargs):
        """Plot data from the detector using this geometry.

        This approximates the geometry to align all pixels to a 2D grid.

        Returns a matplotlib axes object.

        Parameters
        ----------

        data : ndarray
          Should have exactly 3 dimensions, for the modules, then the
          slow scan and fast scan pixel dimensions.
        axis_units : str
          Show the detector scale in pixels ('px') or metres ('m').
        frontview : bool
          If True (the default), x increases to the left, as if you were looking
          along the beam. False gives a 'looking into the beam' view.
        ax : `~matplotlib.axes.Axes` object, optional
          Axes that will be used to draw the image. If None is given (default)
          a new axes object will be created.
        figsize : tuple
          Size of the figure (width, height) in inches to be drawn
          (default: (10, 10))
        colorbar : bool, dict
          Draw colobar with default values (if boolean is given). Colorbar
          appearance can be controlled by passing a dictionary of properties.
        kwargs :
          Additional keyword arguments passed to `~matplotlib.imshow`
        """
        return self._snapped().plot_data(
            data, axis_units=axis_units, frontview=frontview, figsize=figsize,
            ax=ax, colorbar=colorbar, **kwargs
            )

    @classmethod
    def _distortion_array_slice(cls, m, t):
        """Which part of distortion array each tile is.
        """
        # _tile_slice gives the slice for the tile within its module.
        # The distortion array joins the modules along the slow-scan axis, so
        # we need to offset the slow-scan slice to land in the correct module.
        ss_slice_inmod, fs_slice = cls._tile_slice(t)
        mod_px_ss = cls.expected_data_shape[1]
        mod_offset = m * mod_px_ss
        ss_slice = slice(
            ss_slice_inmod.start + mod_offset, ss_slice_inmod.stop + mod_offset
        )
        return ss_slice, fs_slice

    def to_distortion_array(self, allow_negative_xy=False):
        """Generate a distortion array for pyFAI from this geometry.
        """
        nmods, mod_px_ss, mod_px_fs = self.expected_data_shape
        ncorners = self._pixel_corners.shape[1]
        distortion = np.zeros((nmods * mod_px_ss, mod_px_fs, ncorners, 3),
                              dtype=np.float32)

        pixpos = self.get_pixel_positions(centre=False).reshape(
            (nmods * mod_px_ss, mod_px_fs, 3)
        )
        px, py, pz = np.moveaxis(pixpos, -1, 0)

        corner_ss_offsets = self._pixel_corners[0]
        corner_fs_offsets = self._pixel_corners[1]

        for m, mod in enumerate(self.modules, start=0):
            for t, tile in enumerate(mod, start=0):
                ss_unit_x, ss_unit_y, ss_unit_z = tile.ss_vec
                fs_unit_x, fs_unit_y, fs_unit_z = tile.fs_vec

                # Which part of the array is this tile?
                tile_ss_slice, tile_fs_slice = self._distortion_array_slice(m, t)

                # Get coordinates of each pixel's first corner
                # 2D arrays, shape: (64, 128)
                pixel_corner1_x = px[tile_ss_slice,  tile_fs_slice]
                pixel_corner1_y = py[tile_ss_slice,  tile_fs_slice]
                pixel_corner1_z = pz[tile_ss_slice,  tile_fs_slice]

                # Calculate corner coordinates for each pixel
                # 3D arrays, shape: (64, 128, 4)
                corners_x = (
                        pixel_corner1_x[:, :, np.newaxis]
                        + corner_ss_offsets * ss_unit_x
                        + corner_fs_offsets * fs_unit_x
                )
                corners_y = (
                        pixel_corner1_y[:, :, np.newaxis]
                        + corner_ss_offsets * ss_unit_y
                        + corner_fs_offsets * fs_unit_y
                )
                corners_z = (
                        pixel_corner1_z[:, :, np.newaxis]
                        + corner_ss_offsets * ss_unit_z
                        + corner_fs_offsets * fs_unit_z
                )

                # Insert the data into the array
                distortion[tile_ss_slice, tile_fs_slice, :, 0] = corners_z
                distortion[tile_ss_slice, tile_fs_slice, :, 1] = corners_y
                distortion[tile_ss_slice, tile_fs_slice, :, 2] = corners_x

        if not allow_negative_xy:
            # Shift the x & y origin from the centre to the corner
            min_yx = distortion[..., 1:].min(axis=(0, 1, 2))
            distortion[..., 1:] -= min_yx

        return distortion

    @classmethod
    def _tile_slice(cls, tileno):
        """Implement in subclass: which part of module array each tile is.
        """
        raise NotImplementedError

    def _module_coords_to_tile(self, slow_scan, fast_scan):
        """Implement in subclass: positions in module to tile numbers & pos in tile
        """
        raise NotImplementedError

    @classmethod
    def _adjust_pixel_coords(cls, ss_coords, fs_coords, centre):
        """Called by get_pixel_positions; overridden by DSSC"""
        if centre:
            # A pixel is from n to n+1 in each axis, so centres are at n+0.5.
            ss_coords += 0.5
            fs_coords += 0.5

    def get_pixel_positions(self, centre=True):
        """Get the physical coordinates of each pixel in the detector

        The output is an array with shape like the data, with an extra dimension
        of length 3 to hold (x, y, z) coordinates. Coordinates are in metres.

        If centre=True, the coordinates are calculated for the centre of each
        pixel. If not, the coordinates are for the first corner of the pixel
        (the one nearest the [0, 0] corner of the tile in data space).
        """
        out = np.zeros(self.expected_data_shape + (3,), dtype=np.float64)

        # Prepare some arrays to use inside the loop
        pixel_ss_coord, pixel_fs_coord = np.meshgrid(
            np.arange(0, self.frag_ss_pixels, dtype=np.float64),
            np.arange(0, self.frag_fs_pixels, dtype=np.float64),
            indexing='ij'
        )

        # Shift coordinates from corner to centre if requested.
        # This is also where the DSSC subclass shifts odd rows by half a pixel
        self._adjust_pixel_coords(pixel_ss_coord, pixel_fs_coord, centre)

        for m, mod in enumerate(self.modules, start=0):
            for t, tile in enumerate(mod, start=0):
                corner_x, corner_y, corner_z = tile.corner_pos
                ss_unit_x, ss_unit_y, ss_unit_z = tile.ss_vec
                fs_unit_x, fs_unit_y, fs_unit_z = tile.fs_vec

                # Calculate coordinates of each pixel's first corner
                # 2D arrays, shape: (64, 128)
                pixels_x = (
                        corner_x
                        + pixel_ss_coord * ss_unit_x
                        + pixel_fs_coord * fs_unit_x
                )
                pixels_y = (
                        corner_y
                        + pixel_ss_coord * ss_unit_y
                        + pixel_fs_coord * fs_unit_y
                )
                pixels_z = (
                        corner_z
                        + pixel_ss_coord * ss_unit_z
                        + pixel_fs_coord * fs_unit_z
                )

                # Which part of the array is this tile?
                tile_ss_slice, tile_fs_slice = self._tile_slice(t)

                # Insert the data into the array
                out[m, tile_ss_slice, tile_fs_slice, 0] = pixels_x
                out[m, tile_ss_slice, tile_fs_slice, 1] = pixels_y
                out[m, tile_ss_slice, tile_fs_slice, 2] = pixels_z

        return out

    def data_coords_to_positions(self, module_no, slow_scan, fast_scan):
        """Convert data array coordinates to physical positions

        Data array coordinates are how you might refer to a pixel in an array
        of detector data: module number, and indices in the slow-scan and
        fast-scan directions. But coordinates in the two pixel dimensions aren't
        necessarily integers, e.g. if they refer to the centre of a peak.

        module_no, fast_scan and slow_scan should all be numpy arrays of the
        same shape. module_no should hold integers, starting from 0,
        so 0: Q1M1, 1: Q1M2, etc.

        slow_scan and fast_scan describe positions within that module.
        They may hold floats for sub-pixel positions. In both, 0.5 is the centre
        of the first pixel.

        Returns an array of similar shape with an extra dimension of length 3,
        for (x, y, z) coordinates in metres.

        .. seealso::

           :doc:`agipd_geometry` demonstrates using this method.
        """
        assert module_no.shape == slow_scan.shape == fast_scan.shape

        # We want to avoid iterating over the positions in Python.
        # So we assemble arrays of the corner position and step vectors for all
        # tiles, and then use numpy indexing to select the relevant ones for
        # each set of coordinates.
        tiles_corner_pos = np.stack([
            t.corner_pos for m in self.modules for t in m
        ])
        tiles_ss_vec = np.stack([
            t.ss_vec for m in self.modules for t in m
        ])
        tiles_fs_vec = np.stack([
            t.fs_vec for m in self.modules for t in m
        ])

        # Convert coordinates within each module to coordinates in a tile
        tilenos, tile_ss, tile_fs = self._module_coords_to_tile(slow_scan, fast_scan)

        # The indexes of the relevant tiles in the arrays assembled above
        all_tiles_ix = (module_no * self.n_tiles_per_module) + tilenos

        # Select the relevant tile geometry for each set of coordinates
        coords_tile_corner = tiles_corner_pos[all_tiles_ix]
        coords_ss_vec = tiles_ss_vec[all_tiles_ix]
        coords_fs_vec = tiles_fs_vec[all_tiles_ix]

        # Calculate the physical coordinate for each data coordinate
        return coords_tile_corner \
            + (np.expand_dims(tile_ss, -1) * coords_ss_vec) \
            + (np.expand_dims(tile_fs, -1) * coords_fs_vec)

    def offset(self, shift, *, modules=np.s_[:], tiles=np.s_[:]):
        """Move part or all of the detector, making a new geometry.

        By default, this moves all modules & tiles. To move the centre down in
        the image, move the whole geometry *up* relative to it.

        Returns a new geometry object of the same type.

        ::

            # Move the whole geometry up 2 mm (relative to the beam)
            geom2 = geom.shift((0, 2e-3))

            # Move quadrant 1 (modules 0, 1, 2, 3) up 2 mm
            geom2 = geom.shift((0, 2e-3), modules=np.s_[0:4])

            # Move each module by a separate amount
            shifts = np.zeros((16, 3))
            shifts[5] = (0, 2e-3, 0)    # x, y, z for individual modules
            shifts[10] = (0, -1e-3, 0)
            geom2 = geom.shift(shifts)

        Parameters
        ----------

        shift: numpy.ndarray or tuple
          (x, y) or (x, y, z) shift to apply in metres. Can be a single shift
          for all selected modules, a 2D array with a shift per module, or a
          3D array with a shift per tile (``arr[module, tile, xyz]``).
        modules: slice
          Select modules to move; defaults to all modules.
          Like all Python slicing, the end number is excluded, so ``np.s_[:4]``
          moves modules 0, 1, 2, 3.
        tiles: slice
          Select tiles to move within each module; defaults to all tiles.
        """
        shift = np.asarray(shift)
        if not shift.shape[-1] in (2, 3):
            raise ValueError(
                "Shift must be 2D or 3D coordinate(s). Last dimension "
                f"was {shift.shape[-1]}"
            )

        ntiles = max([len(m) for m in self.modules])
        all_shifts = np.zeros((len(self.modules), ntiles, 3), dtype=shift.dtype)
        sel_shifts = all_shifts[modules, tiles, :shift.shape[-1]]

        if shift.shape[:-1] == sel_shifts.shape[:2]:
            # Per-tile offsets
            sel_shifts[:] = shift
        elif shift.shape[:-1] == sel_shifts.shape[:1]:
            # Per-module offsets - broadcast across tiles
            sel_shifts[:] = shift[:, np.newaxis]
        elif shift.shape[:-1] == ():
            # Single shift - broadcast across modules and tiles
            sel_shifts[:] = shift
        else:
            raise ValueError(
                f"Got {shift.shape[:-1]} coordinates. Expected either a single "
                f"coordinate (), a coordinate per module {sel_shifts.shape[:1]} "
                f"or a coordinate per tile {sel_shifts.shape[:2]}"
            )

        cls = type(self)
        return cls([
            [
                tile.offset(all_shifts[m, t])
                for t, tile in enumerate(module)
            ] for m, module in enumerate(self.modules)
        ])


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
    arr = np.zeros((512, 128), dtype=np.bool)
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

        The quadrant positions are not stored in the file, and must be provided
        separately. By default, both the quadrant positions and the positions
        in the file are measured in millimetres; the unit parameter controls
        this.

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

    def to_h5_file_and_quad_positions(self, path):
        """Write this geometry to an XFEL HDF5 format geometry file

        The quadrant positions are not stored in the file, so they are returned
        separately. These and the numbers in the file are in millimetres.

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
    def from_h5_file_and_quad_positions(cls, path, positions, unit=1e-3):
        """Load a DSSC geometry from an XFEL HDF5 format geometry file

        The quadrant positions are not stored in the file, and must be provided
        separately. The position given should refer to the bottom right (looking
        along the beam) corner of the quadrant.

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

    def to_h5_file_and_quad_positions(self, path):
        """Write this geometry to an XFEL HDF5 format geometry file

        The quadrant positions are not stored in the file, so they are returned
        separately. These and the numbers in the file are in millimetres.

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
    def __init__(self, modules, filename='No file'):
        super().__init__(modules, filename)
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

    def __init__(self, modules, filename='No file'):
        super().__init__(modules, filename)
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

          These offsets are positions for the bottom, beam-left corner of each
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
        module_width = 4 * (cls.frag_fs_pixels + asic_gap)
        module_height = 2 * (cls.frag_ss_pixels + asic_gap)
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

    @classmethod
    def from_crystfel_geom(cls, filename):

        raise NotImplementedError

    def write_crystfel_geom(self, file_name):

        raise NotImplementedError


class PNCCDGeometry(DetectorGeometryBase):
    """Detector layout for pnCCD

    The large-area, pn-junction Charge Coupled Device detector consists
    of two movable modules with a single tile each.

    In its default configuration, the complete detector frame is read
    out and written to file as a single image, with the modules split
    along the slow-scan dimension y. The public methods of this type
    support both the combined image array as well as separated module
    with expected_data_shape.
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
