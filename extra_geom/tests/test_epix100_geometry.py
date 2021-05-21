from extra_geom import Epix100Geometry

import extra_geom
import matplotlib.pyplot as plt
import numpy as np
import pytest


@pytest.fixture
def epix():
    return Epix100Geometry.from_origin()


def test_epix_geometry(epix):
    """ Tests that an object was created """
    assert isinstance(epix, extra_geom.detectors.Epix100Geometry)
    assert epix.expected_data_shape == (1, 704, 768)
    assert epix.n_modules == 1
    assert epix.n_quads == 0
    assert epix.n_tiles_per_module == 4
    assert isinstance(epix.modules, list)
    assert len(epix.modules) == 1
    assert isinstance(epix.modules[0], list)
    assert len(epix.modules[0]) == 4
    for i in range(epix.n_tiles_per_module):
        assert isinstance(epix.modules[0][0],
                          extra_geom.detectors.GeometryFragment)


def test__tile_silce(epix):
    _slice = epix._tile_slice(0)
    assert isinstance(_slice, tuple)
    assert isinstance(_slice[0], slice)
    assert isinstance(_slice[1], slice)
    assert _slice == (slice(352, 704), slice(384, 768))


def test_write_crystfel_geom(epix, tmp_path):
    epix.write_crystfel_geom(tmp_path / 'test.geom')


def test_compare(epix):
    ax = epix.compare(epix)
    assert isinstance(ax, plt.Axes)


def test_inspect(epix):
    assert isinstance(epix.inspect(), plt.Axes)


def test_get_pixel_positions(epix):
    out = epix.get_pixel_positions()
    assert isinstance(out, np.ndarray)
    assert out.shape == epix.expected_data_shape + (3,)


def test_split_tiles(epix):
    data = np.zeros(epix.expected_data_shape)
    split_tile = epix.split_tiles(data)
    assert isinstance(split_tile, list)
    assert len(split_tile) == epix.n_tiles_per_module
    s1, s2, s3 = split_tile[0].shape
    assert ((s1, s2 * epix.ss_tiles, s3 * epix.fs_tiles) ==
            epix.expected_data_shape)


def test_module_coords_to_tile(epix):
    """ expected data shape is
        (3, 704 [slow: 2 tiles * 352 px], 768 [fast: 2 tiles * 384 px])

    The points are:
    (5, 20),    <- t=2
    (400, 430), <- t=0
    (30, 450),  <- t=1
    (510, 60),  <- t=3
    """
    slow_scan = np.array([ 5, 400,  30, 510])
    fast_scan = np.array([20, 430, 450,  60])
    tileno, tile_ss, tile_fs = epix._module_coords_to_tile(
        slow_scan, fast_scan)

    np.testing.assert_array_equal(tileno, [2, 0, 1, 3])
    np.testing.assert_array_equal(tile_ss, np.mod(slow_scan, 352))
    np.testing.assert_allclose(tile_fs, np.mod(fast_scan, 384))
