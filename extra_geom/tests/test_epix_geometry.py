import matplotlib.pyplot as plt
import numpy as np
import pytest
from cfelpyutils.crystfel_utils import load_crystfel_geometry

import extra_geom
from extra_geom import Epix10KGeometry, Epix100Geometry

from .utils import assert_geom_close


@pytest.fixture
def args(request):
    return request.getfixturevalue(request.param)


@pytest.fixture
def epix100():
    return (
        Epix100Geometry.from_origin(),
        (Epix100Geometry.frag_ss_pixels, Epix100Geometry.frag_fs_pixels),
        Epix100Geometry.pixel_size,
        extra_geom.detectors.Epix100Geometry)


@pytest.fixture
def epix10K():
    return (
        Epix10KGeometry.from_origin(),
        (Epix10KGeometry.frag_ss_pixels, Epix10KGeometry.frag_fs_pixels),
        Epix10KGeometry.pixel_size,
        extra_geom.detectors.Epix10KGeometry)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_epix_geometry(args):
    """ Tests that an object was created """
    epix, (nrow, ncol), pxsz, cls = args
    assert isinstance(epix, cls)
    assert epix.expected_data_shape == (1, 2*nrow, 2*ncol)
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


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test__tile_silce(args):
    epix, (nrow, ncol), pxsz, cls = args
    _slice = epix._tile_slice(0)
    assert isinstance(_slice, tuple)
    assert isinstance(_slice[0], slice)
    assert isinstance(_slice[1], slice)
    assert _slice == (slice(0, nrow), slice(0, ncol))


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_compare(args):
    epix, (nrow, ncol), pxsz, cls = args
    ax = epix.compare(epix)
    assert isinstance(ax, plt.Axes)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_inspect(args):
    epix, (nrow, ncol), pxsz, cls = args
    assert isinstance(epix.inspect(), plt.Axes)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_get_pixel_positions(args):
    epix, (nrow, ncol), pxsz, cls = args
    out = epix.get_pixel_positions()
    assert isinstance(out, np.ndarray)
    assert out.shape == epix.expected_data_shape + (3,)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_split_tiles(args):
    epix, (nrow, ncol), pxsz, cls = args
    data = np.zeros(epix.expected_data_shape)
    split_tile = epix.split_tiles(data)
    assert isinstance(split_tile, list)
    assert len(split_tile) == epix.n_tiles_per_module
    for tileno in range(4):
        s1, s2, s3 = split_tile[0].shape
        assert ((s1, s2 * epix.ss_tiles, s3 * epix.fs_tiles) ==
                epix.expected_data_shape)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_module_coords_to_tile(args):
    """ expected data shape is
        ePix100: (3, 704 [slow: 2 tiles * 352 px], 768 [fast: 2 tiles * 384 px])
        ePix10K: (3, 352 [slow: 2 tiles * 176 px], 384 [fast: 2 tiles * 192 px])

    The points are:
    (  5,  10),  <- t=0, (5, 10) 
    (200, 215),  <- t=3, (24, 23)
    (15,  225),  <- t=1, (15, 33)
    (255,  30),  <- t=2, (79, 30)

    doubled for ePix100
    """
    epix, (nrow, ncol), pxsz, cls = args
    css, cfs = nrow // 176, ncol // 192

    slow_scan = np.array([ 5, 200,  15, 255])*css
    fast_scan = np.array([10, 215, 225,  30])*cfs
    tileno, tile_ss, tile_fs = epix._module_coords_to_tile(
        slow_scan, fast_scan)

    np.testing.assert_array_equal(tileno, [0, 3, 1, 2])
    np.testing.assert_allclose(tile_ss, np.array([ 5, 24, 15, 79])*css)
    np.testing.assert_allclose(tile_fs, np.array([10, 23, 33, 30])*cfs)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
def test_write_read_crystfel_file(args, tmpdir):
    epix, (nrow, ncol), pxsz, cls = args
    path = str(tmpdir / 'test.geom')
    epix.write_crystfel_geom(filename=path, photon_energy=9000,
                             adu_per_ev=0.0042, clen=0.101)

    loaded = cls.from_crystfel_geom(path)
    assert_geom_close(loaded, epix)

    # Load the geometry file with cfelpyutils and test the rigid groups
    geom_dict = load_crystfel_geometry(path)
    assert geom_dict['panels']['p0a0']['res'] == 1 / pxsz
    assert len(geom_dict['panels']) == 4
    p0a0 = geom_dict['panels']['p0a0']
    assert p0a0['max_ss'] == nrow - 1
    assert p0a0['min_ss'] == 0
    assert p0a0['max_fs'] == ncol - 1
    assert p0a0['min_fs'] == 0


@pytest.mark.parametrize('args', ['epix100'], indirect=True)
def test_from_relative_position(args):
    epix, (nrow, ncol), pxsz, cls = args

    geom = cls.from_relative_positions(
        top=(ncol+cls.asic_gap/2, nrow+25/2, 0),
        bottom=(ncol+cls.asic_gap/2, -25/2, 0)
    )

    np.testing.assert_almost_equal(geom.modules[0][0].corner_pos, [0.019325, 0.018225, 0.0])

    pos = geom.get_pixel_positions(centre=False)
    print(pos.shape)
    assert pos.shape == (1, 704, 768, 3)
    px = pos[..., 0]
    py = pos[..., 1]

    np.testing.assert_almost_equal(px.max(), (ncol+cls.asic_gap/2)*pxsz)
    np.testing.assert_almost_equal(px.min(), -(ncol+cls.asic_gap/2-1)*pxsz)
    np.testing.assert_almost_equal(py.max(), (nrow+12.5)*pxsz)
    np.testing.assert_almost_equal(py.min(), -(nrow+12.5-1)*pxsz)


@pytest.mark.parametrize('args', ['epix100', 'epix10K'], indirect=True)
@pytest.mark.parametrize('shape', [(1,), (33, 1), (33,), tuple()])
def test_ensure_shape(args, shape):
    epix, (nrow, ncol), pxsz, cls = args

    data = np.zeros(shape + epix.expected_data_shape[1:])
    img, centre = epix.position_modules_fast(data)
    expected_img_shape = (
        np.empty(shape).squeeze().shape
        + (epix.frag_ss_pixels * 2 + epix.asic_gap,
           epix.frag_fs_pixels * 2 + epix.asic_gap)
    )
    np.testing.assert_allclose(img.shape, expected_img_shape, atol=1)

    # with 4 extra diagnostic rows
    data = np.zeros(shape + (epix.expected_data_shape[-2] + 4, epix.expected_data_shape[-1]))
    img, centre = epix.position_modules_fast(data)
    expected_img_shape = (
        np.empty(shape).squeeze().shape
        + (epix.frag_ss_pixels * 2 + epix.asic_gap,
           epix.frag_fs_pixels * 2 + epix.asic_gap)
    )
    np.testing.assert_allclose(img.shape, expected_img_shape, atol=1)
