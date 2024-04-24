from functools import lru_cache
from itertools import chain
from warnings import warn

import numpy as np
from cfelpyutils.geometry import load_crystfel_geometry

from .crystfel_fmt import write_crystfel_geom
from .snapped import GridGeometryFragment, SnappedGeometry
from .motors import read_motors_from_geom


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
        corner_pos = np.array([d['cnx']/res, d['cny']/res, d['coffset']])
        ss_vec = np.array([d['ssx'], d['ssy'], d['ssz']]) / res
        fs_vec = np.array([d['fsx'], d['fsy'], d['fsz']]) / res
        ss_pixels = d['orig_max_ss'] - d['orig_min_ss'] + 1
        fs_pixels = d['orig_max_fs'] - d['orig_min_fs'] + 1
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

    def rotate(self, matrix, center):
        ss_vec = self.ss_vec @ matrix
        fs_vec = self.fs_vec @ matrix
        pos = (self.corner_pos - center) @ matrix + center
        return type(self)(pos, ss_vec, fs_vec, self.ss_pixels, self.fs_pixels)

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
    n_quads = 0
    n_modules = 0
    n_tiles_per_module = 0
    expected_data_shape = (0, 0, 0)
    _pixel_corners = np.array([  # pixel units; overridden for DSSC
        [0, 1, 1, 0],  # slow-scan
        [0, 0, 1, 1]   # fast-scan
    ])
    _draw_first_px_on_tile = 1  # Tile num of 1st pixel - overridden for LPD
    _pyfai_cls_name = None  # Name of class in extra_geom.pyfai

    @property
    def _pixel_shape(self):
        """Pixel (x, y) shape. Overridden for DSSC."""
        return np.array([1., 1.], dtype=np.float64) * self.pixel_size

    def __init__(self, modules, filename='No file', metadata=None):
        # List of lists (1 per module) of fragments (1 per tile)
        self.modules = modules
        # self.filename is metadata for plots, we don't read/write the file.
        # There are separate methods for reading and writing.
        self.filename = filename
        self.metadata = metadata if (metadata is not None) else {}
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
        from matplotlib.collections import LineCollection, PatchCollection
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
    def _cfel_panels_by_data_coord(cls, panels: dict):
        """Arrange panel dicts from CrystFEL geometry by first data coordinate

        Index panels by which part of the data they refer to, rather than
        relying on names like p0a0.
        """
        res = {}
        for pname, info in panels.items():
            dims = info['dim_structure']
            ix_dims = [i for i in dims if isinstance(i, int)]
            if len(ix_dims) > 1:
                raise ValueError(f"Too many index dimensions for {pname}: {dims}")

            min_ss = info['orig_min_ss']
            if ix_dims:
                # Geometry for 3D data, modules stacked along separate axis
                modno = ix_dims[0]
            else:
                # Geometry for 2D data, modules concatenated along slow-scan axis
                modno, min_ss = divmod(min_ss, cls.expected_data_shape[1])

            info['panel_name'] = pname
            res[(modno, min_ss, info['orig_min_fs'])] = info

        return res

    @classmethod
    def from_crystfel_geom(cls, filename):
        """Read a CrystFEL format (.geom) geometry file.

        Returns a new geometry object.
        """
        cfel_geom = load_crystfel_geometry(filename)
        panels_by_data_coord = cls._cfel_panels_by_data_coord(
            cfel_geom.detector['panels']
        )
        n_modules = cls.n_modules
        if n_modules == 0:
            # Detector type with varying number of modules (e.g. JUNGFRAU)
            n_modules = max(c[0] for c in panels_by_data_coord) + 1

        modules = []
        panel_names_to_pNaM = {}
        for p in range(n_modules):
            tiles = []
            modules.append(tiles)
            for a in range(cls.n_tiles_per_module):
                ss_slice, fs_slice = cls._tile_slice(a)
                d = panels_by_data_coord[p, ss_slice.start, fs_slice.start]
                tiles.append(GeometryFragment.from_panel_dict(d))
                panel_names_to_pNaM[d['panel_name']] = f'p{p}a{a}'

        # Store some extra fields to write if we create another .geom file.
        # It's possible for these to have different values for different panels,
        # but it seems to be common to use them like headers, describing all
        # panels, and we're assuming that's the case here.
        cfel_md_keys = ('data', 'mask', 'adu_per_eV', 'clen')
        d1 = panels_by_data_coord[0, 0, 0]
        metadata = {'crystfel': {k: d1.get(k) for k in cfel_md_keys}}
        metadata['crystfel']['photon_energy'] = cfel_geom.beam['photon_energy']

        # Normalise description of bad regions, so we can output it correctly.
        # - Change panel names to uniform pNaM (panel N asic M) format
        # - If the file has a 2D layout (modules arranged along the slow-scan
        #   axis), convert slow-scan coordinates to a 3D layout.
        file_geom_is_2d = not any(isinstance(d, int) for d in d1['dim_structure'])
        adjusted_bad_regions = {}
        for bad_name, bad_d in cfel_geom.detector['bad'].items():
            panel_name = bad_d['panel']
            if panel_name:
                try:
                    bad_d['panel'] = panel_names_to_pNaM[panel_name]
                except KeyError:
                    warn("Discarding {bad_name}, no such panel {panel_name!r}")
                    continue
            if bad_d['is_fsss']:
                if not panel_name:
                    warn("Discarding {bad_name}, ss/fs region without panel name")
                    continue
                if file_geom_is_2d:
                    bad_d['min_ss'] %= cls.expected_data_shape[1]
                    bad_d['max_ss'] %= cls.expected_data_shape[1]
            adjusted_bad_regions[bad_name] = bad_d

        metadata['crystfel']['bad'] = adjusted_bad_regions

        geom = cls(modules, filename=filename, metadata=metadata)

        with open(filename, 'r') as f:
            try:
                geom.motor_positions = read_motors_from_geom(f)
            except ValueError:
                # Motors are not found, ignore
                pass

        return geom

    def write_crystfel_geom(self, filename, *,
                            data_path=None,
                            mask_path=None, dims=('frame', 'modno', 'ss', 'fs'),
                            nquads=None, adu_per_ev=None, clen=None,
                            photon_energy=None):
        """Write this geometry to a CrystFEL format (.geom) geometry file.

        If the geometry was read from a ``.geom`` file by
        :meth:`from_crystfel_geom`, some of the optional fields will be filled
        from metadata if not specified.

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
        if nquads is None:
            nquads = self.n_quads

        # If this geometry came from a .geom file, restore metadata from there
        cfelmeta = self.metadata.get('crystfel', {})
        if data_path is None:
            data_path = cfelmeta.get('data') or '/entry_1/instrument_1/detector_1/data'
        if mask_path is None:
            mask_path = cfelmeta.get('mask')
        if adu_per_ev is None:
            adu_per_ev = cfelmeta.get('adu_per_eV')
        if clen is None:
            clen = cfelmeta.get('clen')
        if photon_energy is None:
            photon_energy = cfelmeta.get('photon_energy')

        write_crystfel_geom(
            self, filename, data_path=data_path, mask_path=mask_path, dims=dims,
            bad_regions=cfelmeta.get('bad', {}),
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

    def output_array_for_position(self, extra_shape=(), dtype=np.float32):
        """Make an empty output array to use with position_modules

        You can speed up assembling images by reusing the same output array:
        call this once, and then pass the array as the ``out=`` parameter to
        :meth:`position_modules`. By default, it allocates a new array on
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

    def output_array_for_position_fast(self, extra_shape=(), dtype=np.float32):
        """Deprecated alias for :meth:`output_array_for_position`"""
        return self.output_array_for_position(extra_shape, dtype)

    def position_modules(self, data, out=None, threadpool=None):
        """Assemble data from this detector according to where the pixels are.

        This approximates the geometry to align all pixels to a 2D grid.

        This method was previously called ``position_modules_fast``, and can
        still be called with that name.

        Parameters
        ----------

        data : ndarray or xarray.DataArray
          The last three dimensions should match the modules, then the
          slow scan and fast scan pixel dimensions. If an xarray labelled array
          is given, it must have a 'module' dimension.
        out : ndarray, optional
          An output array to assemble the image into. By default, a new
          array is allocated. Use :meth:`output_array_for_position` to
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

    def position_modules_fast(self, data, out=None, threadpool=None):
        """Deprecated alias for :meth:`position_modules`"""
        return self.position_modules(data, out=out, threadpool=threadpool)

    def position_all_modules(self, data, out=None):
        """Deprecated alias for :meth:`position_modules`"""
        return self.position_modules(data, out=out)

    def position_modules_symmetric(self, data, out=None, threadpool=None):
        """Assemble data with the centre in the middle of the output array.

        The assembly process is the same as :meth:`position_modules`,
        aligning each module to a single pixel grid. But this makes the output
        array symmetric, with the centre at (height // 2, width // 2).

        Parameters
        ----------

        data : ndarray or xarray.DataArray
          The last three dimensions should match the modules, then the
          slow scan and fast scan pixel dimensions. If an xarray labelled array
          is given, it must have a 'module' dimension.
        out : ndarray, optional
          An output array to assemble the image into. By default, a new
          array is created at the minimum size to allow symmetric assembly.
          If an array is passed in, its last two dimensions must be at least
          this size.
        threadpool : concurrent.futures.ThreadPoolExecutor, optional
          If passed, parallelise copying data into the output image.
          See :meth:`position_modules` for details.

        Returns
        -------
        out : ndarray
          Array with one dimension fewer than the input.
          The last two dimensions represent pixel y and x in the detector space.
        """
        return self._snapped().position_modules_symmetric(
            data, out=out, threadpool=threadpool
        )

    def plot_data(self,
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

        This method was previously called ``plot_data_fast``, and can
        still be called with that name.

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

    def plot_data_fast(self, data, **kwargs):
        """Deprecated alias for :meth:`plot_data`"""
        return self.plot_data(data, **kwargs)

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

    def to_pyfai_detector(self):
        """Make a PyFAI detector object for this detector

        You can use PyFAI to azimuthally integrate detector images around
        the centre point of the geometry. The detector object holds the
        positions of all the pixels. See the examples for how to use this.
        """
        if self._pyfai_cls_name is None:
            raise NotImplementedError

        from . import pyfai
        det = getattr(pyfai, self._pyfai_cls_name)()
        det.set_pixel_corners(self.to_distortion_array(allow_negative_xy=True))
        return det

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

           :doc:`convert_coords` demonstrates using this method.
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
            geom2 = geom.offset((0, 2e-3))

            # Move quadrant 1 (modules 0, 1, 2, 3) up 2 mm
            geom2 = geom.offset((0, 2e-3), modules=np.s_[0:4])

            # Move each module by a separate amount
            shifts = np.zeros((16, 3))
            shifts[5] = (0, 2e-3, 0)    # x, y, z for individual modules
            shifts[10] = (0, -1e-3, 0)
            geom2 = geom.offset(shifts)

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

    def rotate(self, angles, center=None, modules=np.s_[:], tiles=np.s_[:],
               degrees=True):
        """Rotate part or all of the detector, making a new geometry.

        The rotation is defined by composition of rotations about the axes of
        the coordinate system (https://en.wikipedia.org/wiki/Euler_angles),
        i.e. a rotation around the z axis rotates the xy (detector) plan.
        We use the right-hand rule, xy being the detector plan with z
        increasing looking toward the detector front plan, x increasing to
        the left, y increasing to the top. Positive rotations are clockwise.

        In other words:

        - Positive x angles tilt the top edge of the detector backwards, away from the source
        - Positive y angles tilt the right-hand edge (looking along the beam) away from the source
        - Positive z angles turn the detector clockwise (looking along the beam)

        By default, this rotates all modules & tiles.
        Returns a new geometry object of the same type.

        ::

            # Rotate the whole geometry by 90 degree in the xy plan
            geom2 = geom.rotate((0, 0, 90))

            # Move the tile 0 in the module 0 around its center by 90 degrees
            geom2 = geom.rotate((0, 0, 90), modules=np.s_[:1], tiles=np.s_[:1])

            # Rotate each module by a separate amount
            rotate = np.zeros((16, 3))
            rotate[5] = (3, 5, 1)  # x, y, z for individual modules
            rotate[10] = (0, -2, 1)
            geom2 = geom.rotate(rotate)

        Parameters
        ----------

        angles: np.array or tuple
            (x, y, z) rotations to apply in degree. Can be a single rotation
            for all selected modules, a 2D array with a rotation per module, or
            a 3D array with a rotation per tile (``arr[module, tile, xyz]``).
        center: np.array or tuple
            center of rotation. Shape must match angles.shape.
            If set to None (default), the rotation center is set to:
            * all modules: center of the detector
            * selected modules: centers of modules
            * selected tiles: centers of tiles
        modules: slice
            Select modules to rotate; defaults to all modules.
        tiles: slice
            Select tiles to move within each module; defaults to all tiles.
        degrees: bool
            If True (default), angles are in degrees. If False, angles are in
            radians.
        """
        if degrees:
            rot = np.deg2rad(angles)
        else:
            rot = np.asarray(angles)

        if rot.shape[-1] != 3:
            raise ValueError('Rotation angles must be 3d (x, y, z). '
                             f'Last dimension was {rot.shape[-1]}')
        if center is not None and np.asarray(center).shape != rot.shape:
            raise ValueError('Rotation center must have same shape as angles')

        ntiles = len(max((m for m in self.modules), key=len))
        all_rot = np.zeros((len(self.modules), ntiles, 3), dtype=rot.dtype)
        sel_rot = all_rot[modules, tiles]

        refs = np.zeros((len(self.modules), ntiles, 3), dtype=np.float64)
        sel_refs = refs[modules, tiles]

        if center is None:
            if rot.shape[:-1] == sel_rot.shape[:2] or tiles != np.s_[:]:
                # Per-tile rotation
                refs[modules, tiles] = np.array(
                    [[t.centre() for t in m[tiles]]
                    for m in self.modules[modules]]
                )
            elif rot.shape[:-1] == sel_rot.shape[:1] or modules != np.s_[:]:
                # Per-module rotation
                mods = []
                for module in self.modules[modules]:
                    center = np.mean([t.centre() for t in module], axis=0)
                    mods.append(np.array([center for _ in module[tiles]]))
                refs[modules, tiles] = np.array(mods)
            else:
                # Single rotation, the reference is the detector center
                # no shift required (refs is already initialized to 0)
                pass
        else:
            center = np.asarray(center)
            if center.shape[:-1] == sel_refs.shape[:2]:
                # Per-tile offsets
                sel_refs[:] = center
            elif center.shape[:-1] == sel_refs.shape[:1]:
                # Per-module offsets - broadcast across tiles
                sel_refs[:] = center[:, np.newaxis]
            elif center.shape[:-1] == ():
                # Single shift - broadcast across modules and tiles
                sel_refs[:] = center
            else:
                raise ValueError(
                    f"Got {center.shape[:-1]} center coordinates. Expected"
                    f"either a single coordinate (), a coordinate per module "
                    f"{sel_refs.shape[:1]} or a coordinate per tile "
                    f"{sel_refs.shape[:2]}"
                )

        if rot.shape[:-1] == sel_rot.shape[:2]:
            # per-tile rotations
            sel_rot[:] = rot
        elif rot.shape[:-1] == sel_rot.shape[:1]:
            # per-module rotations - broadcast across tiles
            sel_rot[:] = rot[:, np.newaxis]
        elif rot.shape[:-1] == ():
            # single rotation - broadcast across modules and tiles
            sel_rot[:] = rot
        else:
            raise ValueError(
                f"Got {rot.shape[:-1]} rotation sequences. Expected either a single "
                f"sequence (), a sequence per module {sel_rot.shape[:1]} "
                f"or a sequence per tile {sel_rot.shape[:2]}"
            )

        @lru_cache()
        def _rot(rotations):
            # generate rotation matrix for the given angles
            # https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
            w, p, k = rotations

            rot_x = np.array([[1, 0, 0],
                              [0, np.cos(w), -np.sin(w)],
                              [0, np.sin(w), np.cos(w)]])
            rot_y = np.array([[np.cos(p), 0, np.sin(p)],
                              [0, 1, 0],
                              [-np.sin(p), 0, np.cos(p)]])
            rot_z = np.array([[np.cos(k), -np.sin(k), 0],
                              [np.sin(k), np.cos(k), 0],
                              [0, 0, 1]])

            return rot_x @ rot_y @ rot_z

        return type(self)([
            [
                tile.rotate(_rot(tuple(all_rot[m, t])), refs[m, t])
                for t, tile in enumerate(module)
            ] for m, module in enumerate(self.modules)
        ])
