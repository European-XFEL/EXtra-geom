from cfelpyutils.crystfel_utils import load_crystfel_geometry
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from matplotlib.axes import Axes
import numpy as np
import pytest

from extra_geom import AGIPD_500K2GGeometry


def test_snap_assemble_data():

    def check_result(img, centre):
        assert img.shape == (602, 1068)
        assert tuple(centre) == (-100, -100)
        assert np.isnan(img[0, 535])
        assert img[50, 50] == 0

    geom = AGIPD_500K2GGeometry.from_origin((-100, -100))

    stacked_data = np.zeros((8, 512, 128))
    img, centre = geom.position_modules_fast(stacked_data)
    check_result(img, centre)

    # test unsafe cast with output array
    stacked_data = np.zeros((8, 512, 128), dtype=np.float64)
    out = geom.output_array_for_position_fast(dtype=np.float32)
    with pytest.raises(TypeError):
        img, centre = geom.position_modules_fast(stacked_data, out=out)

    # test safe cast with output array
    stacked_data = np.zeros((8, 512, 128), dtype=np.uint16)
    out = geom.output_array_for_position_fast(dtype=np.float32)
    img, centre = geom.position_modules_fast(stacked_data, out=out)
    check_result(img, centre)
    check_result(out, centre)
    assert img.dtype == out.dtype == np.float32

    # Assemble in parallel
    stacked_data = np.zeros((8, 512, 128))
    with ThreadPoolExecutor(max_workers=2) as tpool:
        img, centre = geom.position_modules_fast(stacked_data, threadpool=tpool)
        check_result(img, centre)

def test_write_read_crystfel_file(tmpdir):
    geom = AGIPD_500K2GGeometry.from_origin()
    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path, photon_energy=9000,
                             adu_per_ev=0.0075, clen=0.2)

    loaded = AGIPD_500K2GGeometry.from_crystfel_geom(path)
    np.testing.assert_allclose(
        loaded.modules[0][0].corner_pos, geom.modules[0][0].corner_pos
    )
    np.testing.assert_allclose(loaded.modules[0][0].fs_vec, geom.modules[0][0].fs_vec)

    # Load the geometry file with cfelpyutils and test the rigid groups
    geom_dict = load_crystfel_geometry(path)
    quad_gr0 = [  # quadrant: p0a0 ... p7a7
        'p{}a{}'.format(p, a) for p, a in product(range(8), range(8))
    ]
    assert geom_dict['rigid_groups']['p0'] == quad_gr0[:8]
    assert geom_dict['rigid_groups']['p7'] == quad_gr0[-8:]
    assert geom_dict['rigid_groups']['q0'] == quad_gr0
    assert geom_dict['panels']['p0a0']['res'] == 5000  # 5000 pixels/metre
    p3a7 = geom_dict['panels']['p3a7']
    assert p3a7['min_ss'] == 448
    assert p3a7['max_ss'] == 511
    assert p3a7['min_fs'] == 0
    assert p3a7['max_fs'] == 127


def test_write_read_crystfel_file_2d(tmpdir):
    geom = AGIPD_500K2GGeometry.from_origin()
    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path, dims=('frame', 'ss', 'fs'),
                             adu_per_ev=0.0075, clen=0.2)

    loaded = AGIPD_500K2GGeometry.from_crystfel_geom(path)
    np.testing.assert_allclose(
        loaded.modules[0][0].corner_pos, geom.modules[0][0].corner_pos
    )
    np.testing.assert_allclose(loaded.modules[0][0].fs_vec, geom.modules[0][0].fs_vec)

    # Load the geometry file with cfelpyutils and check some values
    geom_dict = load_crystfel_geometry(path)

    p3a7 = geom_dict['panels']['p3a7']
    assert p3a7['dim_structure'] == ['%', 'ss', 'fs']
    assert p3a7['min_ss'] == (3 * 512) + 448
    assert p3a7['max_ss'] == (3 * 512) + 511
    assert p3a7['min_fs'] == 0
    assert p3a7['max_fs'] == 127


def test_inspect():
    geom = AGIPD_500K2GGeometry.from_origin()
    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)


def test_compare():
    geom1 = AGIPD_500K2GGeometry.from_origin((5, 10))
    geom2 = AGIPD_500K2GGeometry.from_origin((7, -5))
    # Smoketest
    ax = geom1.compare(geom2)
    assert isinstance(ax, Axes)


def test_to_distortion_array():
    geom = AGIPD_500K2GGeometry.from_origin()
    # Smoketest
    distortion = geom.to_distortion_array()
    assert isinstance(distortion, np.ndarray)
    assert distortion.shape == (4096, 128, 4, 3)

    # Coordinates in m, origin at corner; max x ~ 21cm & y ~ 12cm
    assert 0.12 < distortion[..., 1].max() < 0.121
    assert 0.21 < distortion[..., 2].max() < 0.22
    assert 0.0 <= distortion[..., 1].min() < 0.01
    assert 0.0 <= distortion[..., 2].min() < 0.01


def test_get_pixel_positions():
    geom = AGIPD_500K2GGeometry.from_origin()

    pixelpos = geom.get_pixel_positions()
    assert pixelpos.shape == (8, 512, 128, 3)
    px = pixelpos[..., 0]
    py = pixelpos[..., 1]

    assert 0. < px.min() < .0002
    assert 0.218 > px.max() > 0.21
    assert 0.13 > py.max() > 0.
    assert  0. < py.min() < 0.0002

def test_data_coords_to_positions():
    geom = AGIPD_500K2GGeometry.from_origin((0, -100))

    module_no = np.full(8, fill_value=6, dtype=np.int16)
    slow_scan = np.linspace(0, 511, num=8, dtype=np.float32)
    fast_scan = np.zeros(8, dtype=np.float32)

    res = geom.data_coords_to_positions(module_no, slow_scan, fast_scan)

    assert res.shape == (8, 3)
    print(res)

    resx, resy, resz = res.T

    np.testing.assert_allclose(resz, 0)
    np.testing.assert_allclose(resy, 100 * geom.pixel_size)

    assert (np.diff(resx) < 0).all()   # Monotonically decreasing
    np.testing.assert_allclose(resx[0], 526 * geom.pixel_size)
    assert -0.01 < resx[-1] < 0.01
