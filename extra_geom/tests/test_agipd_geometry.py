from concurrent.futures import ThreadPoolExecutor
from itertools import product

import numpy as np
import pyFAI.detectors
import pytest
from cfelpyutils.geometry import load_crystfel_geometry
from matplotlib.axes import Axes

from extra_geom import AGIPD_1MGeometry, agipd_asic_seams
from .utils import assert_geom_close

def test_snap_assemble_data():

    def check_result(img, centre):
        assert img.shape == (1256, 1092)
        assert tuple(centre) == (631, 550)
        assert np.isnan(img[0, 0])
        assert img[50, 50] == 0

    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )

    stacked_data = np.zeros((16, 512, 128))
    img, centre = geom.position_modules_fast(stacked_data)
    check_result(img, centre)

    # test unsafe cast with output array
    stacked_data = np.zeros((16, 512, 128), dtype=np.float64)
    out = geom.output_array_for_position_fast(dtype=np.float32)
    with pytest.raises(TypeError):
        img, centre = geom.position_modules_fast(stacked_data, out=out)

    # test safe cast with output array
    stacked_data = np.zeros((16, 512, 128), dtype=np.uint16)
    out = geom.output_array_for_position_fast(dtype=np.float32)
    img, centre = geom.position_modules_fast(stacked_data, out=out)
    check_result(img, centre)
    check_result(out, centre)
    assert img.dtype == out.dtype == np.float32

    # Assemble in parallel
    stacked_data = np.zeros((16, 512, 128))
    with ThreadPoolExecutor(max_workers=2) as tpool:
        img, centre = geom.position_modules_fast(stacked_data, threadpool=tpool)
        check_result(img, centre)


def test_assemble_symmetric():
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )
    print("2 centre", geom._snapped().centre * 2)

    stacked_data = np.zeros((16, 512, 128))
    img = geom.position_modules_symmetric(stacked_data)

    assert img.shape == (1262, 1100)
    assert np.isnan(img[0, 0])
    assert np.isnan(img[img.shape[0] // 2, img.shape[1] // 2])
    assert img[50, 50] == 0

    # Smoketest assembling into suitable output array
    geom.position_modules_symmetric(stacked_data, out=img)

    with pytest.raises(ValueError):
        # Output array not big enough
        geom.position_modules_symmetric(stacked_data, out=img[:-1, :-1])


bad_xy = {
    'is_fsss': False,
    'min_x': -20., 'max_x': 20., 'min_y': -100., 'max_y': 100,
    'min_fs': 0, 'max_fs': 0, 'min_ss': 0, 'max_ss': 0, 'panel': '',
}
bad_fsss = {
    'is_fsss': True,
    'min_x': None, 'max_x': None, 'min_y': None, 'max_y': None,
    'min_fs': 10, 'max_fs': 100, 'min_ss': 450, 'max_ss': 500, 'panel': 'p3a7',
}
def assert_bad_region_like(actual, expected):
    actual = {k: None if isinstance(v, float) and np.isnan(v) else v
              for (k, v) in actual.items()}
    assert actual == expected


def test_write_read_crystfel_file(tmpdir):
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )
    # Use the z dimension (coffset in .geom)
    geom = geom.offset((0, 0, 0.001), modules=np.s_[8:12])

    # Add some bad regions in CrystFEL file
    geom.metadata['crystfel'] = {'bad': {'bad_xy': bad_xy, 'bad_fsss': bad_fsss}}

    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path, photon_energy=9000,
                             adu_per_ev=0.0075, clen=0.2)

    loaded = AGIPD_1MGeometry.from_crystfel_geom(path)
    assert_geom_close(loaded, geom)
    assert_bad_region_like(loaded.metadata['crystfel']['bad']['bad_xy'], bad_xy)
    assert_bad_region_like(loaded.metadata['crystfel']['bad']['bad_fsss'], bad_fsss)

    # Load the geometry file with cfelpyutils and test the rigid groups
    geom_dict = load_crystfel_geometry(path).detector
    quad_gr0 = [  # 1st quadrant: p0a0 ... p3a7
        'p{}a{}'.format(p, a) for p, a in product(range(4), range(8))
    ]
    assert geom_dict['rigid_groups']['p0'] == quad_gr0[:8]
    assert geom_dict['rigid_groups']['p3'] == quad_gr0[-8:]
    assert geom_dict['rigid_groups']['q0'] == quad_gr0
    assert geom_dict['panels']['p0a0']['res'] == 5000  # 5000 pixels/metre
    p3a7 = geom_dict['panels']['p3a7']
    assert p3a7['orig_min_ss'] == 448
    assert p3a7['orig_max_ss'] == 511
    assert p3a7['orig_min_fs'] == 0
    assert p3a7['orig_max_fs'] == 127

    print(geom_dict['bad'])
    assert_bad_region_like(geom_dict['bad']['bad_xy'], bad_xy)
    assert_bad_region_like(geom_dict['bad']['bad_fsss'], bad_fsss)


def test_write_read_crystfel_file_2d(tmpdir):
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )
    path = str(tmpdir / 'test.geom')

    # Add some bad regions in CrystFEL file
    geom.metadata['crystfel'] = {'bad': {'bad_xy': bad_xy, 'bad_fsss': bad_fsss}}

    geom.write_crystfel_geom(filename=path, dims=('frame', 'ss', 'fs'),
                             adu_per_ev=0.0075, clen=0.2)

    loaded = AGIPD_1MGeometry.from_crystfel_geom(path)
    assert_geom_close(loaded, geom)
    assert_bad_region_like(loaded.metadata['crystfel']['bad']['bad_xy'], bad_xy)
    assert_bad_region_like(loaded.metadata['crystfel']['bad']['bad_fsss'], bad_fsss)

    # Load the geometry file with cfelpyutils and check some values
    geom_dict = load_crystfel_geometry(path).detector

    p3a7 = geom_dict['panels']['p3a7']
    assert p3a7['dim_structure'] == ['%', 'ss', 'fs']
    assert p3a7['orig_min_ss'] == (3 * 512) + 448
    assert p3a7['orig_max_ss'] == (3 * 512) + 511
    assert p3a7['orig_min_fs'] == 0
    assert p3a7['orig_max_fs'] == 127

    assert geom_dict['bad']['bad_fsss']['min_ss'] == (3 * 512) + 450
    assert geom_dict['bad']['bad_fsss']['max_ss'] == (3 * 512) + 500


def test_quad_positions():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)

    np.testing.assert_allclose(geom.quad_positions(), quad_pos)


def test_inspect():
    geom = AGIPD_1MGeometry.example()
    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)


def test_compare():
    geom1 = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )
    geom2 = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-527, 625), (-548, -10), (520, -162), (542.5, 473)]
    )
    # Smoketest
    ax = geom1.compare(geom2)
    assert isinstance(ax, Axes)


def test_to_distortion_array():
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )
    # Smoketest
    distortion = geom.to_distortion_array()
    assert isinstance(distortion, np.ndarray)
    assert distortion.shape == (8192, 128, 4, 3)

    # Coordinates in m, origin at corner; max x & y should be ~ 25cm
    assert 0.20 < distortion[..., 1].max() < 0.30
    assert 0.20 < distortion[..., 2].max() < 0.30
    assert 0.0 <= distortion[..., 1].min() < 0.01
    assert 0.0 <= distortion[..., 2].min() < 0.01

def test_get_pixel_positions():
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )

    pixelpos = geom.get_pixel_positions()
    assert pixelpos.shape == (16, 512, 128, 3)
    px = pixelpos[..., 0]
    py = pixelpos[..., 1]

    assert -0.12 < px.min() < -0.1
    assert  0.12 > px.max() > 0.1
    assert -0.14 < py.min() < -0.12
    assert  0.14 > py.max() >  0.12

def test_data_coords_to_positions():
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )

    module_no = np.zeros(16, dtype=np.int16)
    slow_scan = np.linspace(0, 500, num=16, dtype=np.float32)
    fast_scan = np.zeros(16, dtype=np.float32)

    res = geom.data_coords_to_positions(module_no, slow_scan, fast_scan)

    assert res.shape == (16, 3)

    resx, resy, resz = res.T

    np.testing.assert_allclose(resz, 0)
    np.testing.assert_allclose(resy, 625 * geom.pixel_size)

    assert (np.diff(resx) > 0).all()   # Monotonically increasing
    np.testing.assert_allclose(resx[0], -525 * geom.pixel_size)
    assert -0.01 < resx[-1] < 0.01


def test_asic_seams():
    arr = agipd_asic_seams()
    assert arr.shape == (512, 128)
    assert not arr[30, 0]
    assert arr[63, 0]


def test_to_pyfai_detector():
    geom = AGIPD_1MGeometry.from_quad_positions(
        quad_pos=[(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    )
    agipd_pyfai = geom.to_pyfai_detector()
    assert isinstance(agipd_pyfai, pyFAI.detectors.Detector)
