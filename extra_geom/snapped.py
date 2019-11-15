"""This deals with geometry approximated to fit a regular 2D pixel grid

This is useful to copy data from multiple modules into a single image array.
This module is not a public API: it's used internally by the classes in
extra_geom.detectors.
"""

from copy import copy
import numpy as np


class GridGeometryFragment:
    """Holds the 2D axis-aligned position and orientation of one detector tile.

    This is used in 'snapped' geometry which efficiently assembles a detector
    image into a 2D array.

    These coordinates are all (y, x), suitable for indexing a numpy array.

    ss_vec and fs_vec must be length 1 vectors in either positive or negative
    x or y direction. In the output array, the fast scan dimension is always x.
    So if the input data is oriented with fast-scan vertical, we need to
    transpose it first.

    Regardless of transposition, we may also need to flip the data on one or
    both axes; the fs_order and ss_order variables handle this.
    """
    def __init__(self, corner_pos, ss_vec, fs_vec, ss_pixels, fs_pixels):
        self.ss_vec = ss_vec
        self.fs_vec = fs_vec
        self.ss_pixels = ss_pixels
        self.fs_pixels = fs_pixels

        if fs_vec[0] == 0:
            # Fast scan is x dimension: Flip without transposing
            fs_order = fs_vec[1]
            ss_order = ss_vec[0]
            self.transform = lambda arr: arr[..., ::ss_order, ::fs_order]
            corner_shift = np.array([
                min(ss_order, 0) * self.ss_pixels,
                min(fs_order, 0) * self.fs_pixels
            ])
            self.pixel_dims = np.array([self.ss_pixels, self.fs_pixels])
        else:
            # Fast scan is y : Transpose so fast scan -> x and then flip
            fs_order = fs_vec[0]
            ss_order = ss_vec[1]
            self.transform = lambda arr: arr.swapaxes(-1, -2)[..., ::fs_order, ::ss_order]
            corner_shift = np.array([
                min(fs_order, 0) * self.fs_pixels,
                min(ss_order, 0) * self.ss_pixels
            ])
            self.pixel_dims = np.array([self.fs_pixels, self.ss_pixels])
        self.corner_idx = corner_pos + corner_shift
        self.opp_corner_idx = self.corner_idx + self.pixel_dims


class SnappedGeometry:
    """Detector geometry approximated to align modules to a 2D grid

    The coordinates used in this class are (y, x) suitable for indexing a
    Numpy array; this does not match the (x, y, z) coordinates in the more
    precise geometry above.
    """
    def __init__(self, modules, geom):
        self.modules = modules
        self.geom = geom
        self.size_yx, self.centre = self._get_dimensions()

    def make_output_array(self, extra_shape=(), dtype=np.float32):
        """Make an output array for self.position_modules()
        """
        shape = extra_shape + self.size_yx
        return np.full(shape, np.nan, dtype=dtype)

    def position_modules(self, data, out=None):
        """Implementation for position_modules_fast
        """
        assert data.shape[-3:] == self.geom.expected_data_shape
        if out is None:
            out = self.make_output_array(data.shape[:-3], data.dtype)
        else:
            assert out.shape == data.shape[:-3] + self.size_yx
            if not np.can_cast(data.dtype, out.dtype, casting='safe'):
                raise TypeError("{} cannot be safely cast to {}".
                                format(data.dtype, out.dtype))

        for i, module in enumerate(self.modules):
            mod_data = data[..., i, :, :]
            tiles_data = self.geom.split_tiles(mod_data)
            for j, tile in enumerate(module):
                tile_data = tiles_data[j]
                # Offset by centre to make all coordinates positive
                y, x = tile.corner_idx + self.centre
                h, w = tile.pixel_dims
                out[..., y : y + h, x : x + w] = tile.transform(tile_data)

        return out, self.centre

    def _get_dimensions(self):
        """Calculate appropriate array dimensions for assembling data.

        Returns (size_y, size_x), (centre_y, centre_x)
        """
        corners = []
        for module in self.modules:
            for tile in module:
                corners.append(tile.corner_idx)
                corners.append(tile.opp_corner_idx)
        corners = np.stack(corners)

        # Find extremes
        min_yx = corners.min(axis=0)
        max_yx = corners.max(axis=0)

        size = max_yx - min_yx
        centre = -min_yx
        return tuple(size), centre

    def plot_data(self,
                  modules_data, *,
                  axis_units='px',
                  frontview=True,
                  ax=None,
                  figsize=None,
                  colorbar=False,
                  **kwargs):
        """Implementation for plot_data_fast
        """
        from matplotlib.cm import viridis
        import matplotlib.pyplot as plt

        if axis_units not in {'px', 'm'}:
            raise ValueError("axis_units must be 'px' or 'm', not {!r}"
                             .format(axis_units))

        res, centre = self.position_modules(modules_data)
        min_y, min_x = -centre
        max_y, max_x = np.array(res.shape) - centre

        _extent = np.array((min_x - 0.5, max_x + 0.5, min_y - 0.5, max_y + 0.5))
        cross_size = 20
        if axis_units == 'm':
            _extent *= self.geom.pixel_size
            cross_size *= self.geom.pixel_size

        # Use a dark grey for missing data
        _cmap = copy(viridis)
        _cmap.set_bad('0.25', 1.0)

        kwargs.setdefault('cmap', _cmap)
        kwargs.setdefault('extent', _extent)
        kwargs.setdefault('origin', 'lower')

        if ax is None:
            fig = plt.figure(figsize=figsize or (10, 10))
            ax = fig.add_subplot(1, 1, 1)

        im = ax.imshow(res, **kwargs)
        if isinstance(colorbar, dict) or colorbar is True:
            if isinstance(colorbar, bool):
                colorbar = {}
            colorbar = plt.colorbar(im, ax=ax, **colorbar)

        ax.set_xlabel('metres' if axis_units == 'm' else 'pixels')
        ax.set_ylabel('metres' if axis_units == 'm' else 'pixels')

        if frontview:
            ax.invert_xaxis()

        # Draw a cross at the centre
        ax.hlines(0, -cross_size, +cross_size, colors='w', linewidths=1)
        ax.vlines(0, -cross_size, +cross_size, colors='w', linewidths=1)
        return ax
