from os.path import abspath, dirname

import numpy as np
import pyFAI.detectors
import pytest
from cfelpyutils.geometry import load_crystfel_geometry
from matplotlib.axes import Axes

from extra_geom import LPD_MiniGeometry


tests_dir = dirname(abspath(__file__))


def test_write_read_crystfel_file(tmpdir):
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])
    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path, clen=0.119, adu_per_ev=0.0075)

    # Load the geometry file with cfelpyutils and test the ridget groups
    loaded = LPD_MiniGeometry.from_crystfel_geom(path)
    np.testing.assert_allclose(
        loaded.modules[0][0].corner_pos, geom.modules[0][0].corner_pos
    )
    np.testing.assert_allclose(loaded.modules[0][0].fs_vec, geom.modules[0][0].fs_vec)

    geom_dict = load_crystfel_geometry(path).detector
    assert geom_dict['rigid_groups']['p0'] == ['p0a0', 'p0a1']


def test_offset():
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])
    y_orig = np.array([m[0].corner_pos[1] for m in geom.modules])

    # Shift module, all tiles
    all_shifted = geom.offset((0, 1e-3))
    y1 = np.array([m[0].corner_pos[1] for m in all_shifted.modules])
    np.testing.assert_allclose(y1, y_orig + 1e-3)

    # Per-tile shift
    shift = np.zeros((1, 2, 3), dtype=np.float64)
    shift[0, 1, 1] = 3e-3  # Shift T1 in y dimension
    t1_shifted = geom.offset(shift)
    y_t0 = np.array([m[0].corner_pos[1] for m in t1_shifted.modules])
    np.testing.assert_allclose(y_t0, y_orig)
    y_t1 = np.array([m[1].corner_pos[1] for m in t1_shifted.modules])
    y_t1_orig = np.array([m[1].corner_pos[1] for m in geom.modules])
    np.testing.assert_allclose(y_t1, y_t1_orig + 3e-3)

    # Wrong number of modules
    with pytest.raises(ValueError):
        geom.offset(np.zeros((15, 2)))

    # Coordinates must be 2D or 3D
    with pytest.raises(ValueError):
        geom.offset(np.zeros(4))


def test_rotate():
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])

    # Uniform rotation for all modules, all tiles
    all_rotated = geom.rotate((90, 0, 0))
    y = np.array([m[0].corner_pos[1] for m in all_rotated.modules])
    np.testing.assert_allclose(y, 0., atol=1e-7)

    all_rotated = geom.rotate((0, 90, 0))
    x = np.array([m[0].corner_pos[0] for m in all_rotated.modules])
    np.testing.assert_allclose(x, 0., atol=1e-7)

    all_rotated = geom.rotate((0, 0, 90))
    assert all_rotated.output_array_for_position_fast().T.shape == \
        geom.output_array_for_position_fast().shape

    all_rotated = geom.rotate((90, 0, 0), center=(0, 1, 0))
    y = np.array([m[0].corner_pos[1] for m in all_rotated.modules])
    np.testing.assert_allclose(y, 1)

    # Per-tile rotation
    rot = np.zeros((1, 2, 3), dtype=np.float64)
    rot[:, 1, :] = (0, 0, 180)  # rotate T1 of each module
    rotated = geom.rotate(rot)
    for mod_rot, mod_ref in zip(rotated.modules, geom.modules):
        np.testing.assert_array_almost_equal(
            mod_rot[1].corners(),
            np.roll(mod_ref[1].corners(), 2, 0)
        )

    # Wrong number of modules
    with pytest.raises(ValueError):
        geom.rotate(np.ones((15, 3)))

    # angles must be 3D
    with pytest.raises(ValueError):
        geom.rotate((1, 1))

    # angles in radian
    deg = geom.rotate((90, 0, 0), modules=np.s_[:1], tiles=np.s_[:1])
    rad = geom.rotate((np.pi / 2, 0, 0), modules=np.s_[:1], tiles=np.s_[:1], degrees=False)
    np.testing.assert_array_almost_equal(
        deg.modules[0][1].corners(),
        rad.modules[0][1].corners()
    )


def test_inspect():
    geom = LPD_MiniGeometry.example(4)
    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)


def test_snap_assemble_data():
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])

    stacked_data = np.zeros((1, 32, 256))
    img, centre = geom.position_modules_fast(stacked_data)
    assert img.shape == (68, 128)
    assert tuple(centre) == (0, 128)
    assert img[50, 50] == 0
    assert np.isnan(img[34, 50])  # Gap between modules


def test_to_distortion_array():
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])
    # Smoketest
    distortion = geom.to_distortion_array()
    assert isinstance(distortion, np.ndarray)
    assert distortion.shape == (32, 256, 4, 3)

    # Coordinates in m, origin at corner; max x & y should be ~ 63mm and 32mm
    assert 0.032 < distortion[..., 1].max() < 0.036
    assert 0.06 < distortion[..., 2].max() < 0.065
    assert 0.0 <= distortion[..., 1].min() < 0.01
    assert 0.0 <= distortion[..., 2].min() < 0.01


def test_data_coords_to_positions():
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])

    module_no = np.zeros(4, dtype=np.int16)
    slow_scan = np.array([12, 24, 12, 24])
    fast_scan = np.array([30, 100, 158, 228])

    tileno, tile_ss, tile_fs = geom._module_coords_to_tile(slow_scan, fast_scan)
    np.testing.assert_allclose(tileno, [0, 0, 1, 1])
    np.testing.assert_allclose(tile_ss, [12, 24, 12, 24])
    np.testing.assert_allclose(tile_fs, [30, 100, 30, 100])

    res = geom.data_coords_to_positions(module_no, slow_scan, fast_scan)

    assert res.shape == (4, 3)

    resx, resy, resz = res.T

    np.testing.assert_allclose(resz, 0)

    # T0 is read in the opposite direction to our coordinate axes, i.e.
    # top to bottom and left to right (viewed from the front)
    np.testing.assert_allclose(resx[1] - resx[0], -35e-3)
    np.testing.assert_allclose(resy[1] - resy[0], -6e-3)
    # And T1 is the opposite, read matching our coordinate directions.
    np.testing.assert_allclose(resx[3] - resx[2], 35e-3)
    np.testing.assert_allclose(resy[3] - resy[2], 6e-3)
    assert 0.001 > resx.max() > resx.min() > -0.065


def test_to_pyfai_detector():
    geom = LPD_MiniGeometry.from_module_positions([(0, 0)])
    agipd_pyfai = geom.to_pyfai_detector()
    assert isinstance(agipd_pyfai, pyFAI.detectors.Detector)
