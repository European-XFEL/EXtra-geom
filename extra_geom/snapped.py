"""This deals with geometry approximated to fit a regular 2D pixel grid

This is useful to copy data from multiple modules into a single image array.
This module is not a public API: it's used internally by the classes in
extra_geom.detectors.
"""

from copy import copy
from itertools import chain
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
            self.pixel_dims = (self.ss_pixels, self.fs_pixels)
        else:
            # Fast scan is y : Transpose so fast scan -> x and then flip
            fs_order = fs_vec[0]
            ss_order = ss_vec[1]
            self.transform = lambda arr: arr.swapaxes(-1, -2)[..., ::fs_order, ::ss_order]
            corner_shift = np.array([
                min(fs_order, 0) * self.fs_pixels,
                min(ss_order, 0) * self.ss_pixels
            ])
            self.pixel_dims = (self.fs_pixels, self.ss_pixels)
        self.corner_idx = tuple(corner_pos + corner_shift)

    def offset(self, y_x) -> 'GridGeometryFragment':
        new = copy(self)
        new.corner_idx = tuple(np.array(self.corner_idx) + y_x)
        return new


class SnappedGeometry:
    """Detector geometry approximated to align modules to a 2D grid

    The coordinates used in this class are (y, x) suitable for indexing a
    Numpy array; this does not match the (x, y, z) coordinates in the more
    precise geometry above.
    """
    def __init__(self, modules, geom, centre):
        self.modules = modules
        self.geom = geom
        self.centre = centre

        # The fragments here are already shifted so corner_idx starts from 0 in
        # each dim, so the max outer edges define the output image size.
        self.size_yx = tuple(np.max([
            np.array(frag.corner_idx) + np.array(frag.pixel_dims)
            for frag in chain(*modules)
        ], axis=0))

    def make_output_array(self, extra_shape=(), dtype=np.float32):
        """Make an output array for self.position_modules()
        """
        shape = extra_shape + self.size_yx
        return np.full(shape, np.nan, dtype=dtype)

    def position_modules(self, data, out=None, threadpool=None):
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

        copy_pairs = []
        for i, module in enumerate(self.modules):
            mod_data = data[..., i, :, :]
            tiles_data = self.geom.split_tiles(mod_data)
            for tile, tile_data in zip(module, tiles_data):
                y, x = tile.corner_idx
                h, w = tile.pixel_dims

                copy_pairs.append((
                    out[..., y : y + h, x : x + w], tile.transform(tile_data)
                ))

        if threadpool is not None:
            def copy_data(pair):
                dst, src = pair
                dst[:] = src
            # concurrent.futures map() is async, so call list() to wait for it
            list(threadpool.map(copy_data, copy_pairs))
        else:
            for dst, src in copy_pairs:
                dst[:] = src

        return out, self.centre

    def position_modules_symmetric(self, data, out=None, threadpool=None):
        """Assemble data so the centre is in the middle of the output array"""
        assert data.shape[-3:] == self.geom.expected_data_shape
        min_shape = np.stack([self.centre * 2, self.size_yx]).max(axis=0)
        if out is None:
            img_shape = min_shape
            out = np.full(data.shape[:-3] + tuple(img_shape), np.nan, dtype=data.dtype)
        else:
            assert out.shape[:-2] == data.shape[:-3]
            img_shape = np.array(out.shape[-2:])
            if (img_shape < np.array(min_shape)).any():
                raise ValueError(
                    f"Output shape {img_shape} less than required {min_shape}"
                )

        y, x = (img_shape // 2) - self.centre  # Find offset
        h, w = self.size_yx
        self.position_modules(
            data, out[..., y:y+h, x:x+w], threadpool=threadpool
        )
        return out

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
            _extent[:2] *= self.geom._pixel_shape[0]  # x
            _extent[2:] *= self.geom._pixel_shape[1]  # y
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
