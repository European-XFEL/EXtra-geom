from extra_geom import GenericGeometry

import extra_geom
import matplotlib.pyplot as plt
import numpy as np
import pytest

pixel_size = 0.01

simple_config = {'pixel_size': pixel_size,
                 'fast_pixels': 128,
                 'slow_pixels': 64,
                 'corner_coordinates': [pixel_size * np.array([5, 1, 0]),
                                        pixel_size * np.array([-133, 5, 0]),
                                        pixel_size * np.array([-63, -130, 0])],
                 'n_tiles_per_module': 2,
                 'fs_vec': np.array([0, 1, 0]),
                 'ss_vec': np.array([1, 0, 0]),
                 'tile_gap': 2 * pixel_size,
                 'tile_vec': [0, -1, 0]
                 }


@pytest.fixture
def simple():
    return GenericGeometry(**simple_config)


def test_simple_geometry(simple):
    """ Tests that an object was created """
    assert isinstance(simple, extra_geom.detectors.GenericGeometry)
    assert simple.expected_data_shape == (3, 64, 256)


def test__tile_silce(simple):
    _slice = simple._tile_slice(1)
    assert isinstance(_slice, tuple)
    assert isinstance(_slice[0], slice)
    assert isinstance(_slice[1], slice)
    assert _slice == (slice(0, 64), slice(128, 256))


def test_write_crystfel_geom(simple):
    simple.write_crystfel_geom('test.geom')


def test_compare(simple):
    ax = simple.compare(simple)
    assert isinstance(ax, plt.Axes)


def test_inspect(simple):
    assert isinstance(simple.inspect(), plt.Axes)


def test_get_pixel_positions(simple):
    out = simple.get_pixel_positions()
    assert isinstance(out, np.ndarray)
    assert out.shape == simple.expected_data_shape + tuple([3])


def test_split_tiles(simple):
    data = np.zeros(simple.expected_data_shape)
    split_tile = simple.split_tiles(data)
    assert isinstance(split_tile, list)
    assert len(split_tile) == simple.n_tiles_per_module
    s1, s2, s3 = split_tile[0].shape
    assert (s1, s2 * simple.ss_tiles, s3 * simple.fs_tiles) == simple.expected_data_shape


def test_module_coords_to_tile(simple):
    """ expected data shape is (3, 64 [slow], 256 [fast: 2 tiles * 128 pixels])
    The points are:
    (5, 20),    <- t=0
    (10, 50),   <- t=0
    (30, 130),  <- t=1
    (60, 200),  <- t=1
    """
    slow_scan = np.array([5, 10, 30, 60])
    fast_scan = np.array([10, 50, 130, 200])
    tileno, tile_ss, tile_fs = simple._module_coords_to_tile(slow_scan, fast_scan)

    assert all(tileno == [0, 0, 1, 1])
    assert all(tile_ss == slow_scan)
    assert all(tile_fs == np.mod(fast_scan, 128))

