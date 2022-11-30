from os.path import abspath, dirname
from os.path import join as pjoin

import h5py
import numpy as np
import pyFAI.detectors
import pytest
from cfelpyutils.geometry import load_crystfel_geometry
from matplotlib.axes import Axes

from extra_geom import LPD_MiniGeometry
from extra_geom.detectors import invert_xfel_lpd_geom


tests_dir = dirname(abspath(__file__))


def test_write_read_crystfel_file(tmpdir):
    geom = LPD_MiniGeometry.from_origin()
    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path)

    with open(path, 'r') as f:
        contents = f.read()
    with open(path, 'w') as f:
        f.write('clen = 0.119\n')
        f.write('adu_per_eV = 0.0075\n')

        f.write(contents)
    # Load the geometry file with cfelpyutils and test the ridget groups
    loaded = LPD_MiniGeometry.from_crystfel_geom(path)
    np.testing.assert_allclose(
        loaded.modules[0][0].corner_pos, geom.modules[0][0].corner_pos
    )
    np.testing.assert_allclose(loaded.modules[0][0].fs_vec, geom.modules[0][0].fs_vec)


    geom_dict = load_crystfel_geometry(path).detector
    print(geom_dict['rigid_groups'].keys())
    quad_gr0 = ['p0a0', 'p0a1', 'p0a2', 'p0a3', 'p0a4', 'p0a5', 'p0a6', 'p0a7',
                'p0a8', 'p0a9', 'p0a10','p0a11', 'p0a12', 'p0a13', 'p0a14',
                'p0a15']
    assert geom_dict['rigid_groups']['p0'] == quad_gr0[:16]
    assert geom_dict['rigid_groups']['q0'] == quad_gr0


def test_offset():
    geom = LPD_MiniGeometry.from_origin()
    y_orig = np.array([m[0].corner_pos[1] for m in geom.modules])

    # Shift module, all tiles
    all_shifted = geom.offset((0, 1e-3))
    y1 = np.array([m[0].corner_pos[1] for m in all_shifted.modules])
    np.testing.assert_allclose(y1, y_orig + 1e-3)

    # Per-tile shift
    shift = np.zeros((1, 16, 3), dtype=np.float64)
    shift[0, 5, 1] = 3e-3  # Shift T6 in y dimension
    t2_shifted = geom.offset(shift)
    y_t1 = np.array([m[0].corner_pos[1] for m in t2_shifted.modules])
    np.testing.assert_allclose(y_t1, y_orig)
    y_t6 = np.array([m[5].corner_pos[1] for m in t2_shifted.modules])
    y_t6_orig = np.array([m[5].corner_pos[1] for m in geom.modules])
    np.testing.assert_allclose(y_t6[4:8], y_t6_orig[4:8] + 3e-3)

    # Wrong number of modules
    with pytest.raises(ValueError):
        geom.offset(np.zeros((15, 2)))

    # Offsets for 16 modules, but only 4 selected
    with pytest.raises(ValueError):
        geom.offset(np.zeros((16, 16, 3)), modules=np.s_[:4])

    # Coordinates must be 2D or 3D
    with pytest.raises(ValueError):
        geom.offset(np.zeros(4))


def test_rotate():
    geom = LPD_MiniGeometry.from_origin()

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
    rot = np.zeros((1, 16, 3), dtype=np.float64)
    rot[:, 5, :] = (0, 0, 180)  # rotate T6 of each module
    rotated = geom.rotate(rot)
    for mod_rot, mod_ref in zip(rotated.modules[4:8], geom.modules[4:8]):
        np.testing.assert_array_almost_equal(
            mod_rot[5].corners(),
            np.roll(mod_ref[5].corners(), 2, 0)
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
    geom = LPD_MiniGeometry.from_origin((11.4, 299))
    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)


def test_snap_assemble_data():
    geom = LPD_MiniGeometry.from_origin()

    stacked_data = np.zeros((1, 256, 256))
    img, centre = geom.position_modules_fast(stacked_data)
    assert img.shape == (284, 260)
    assert tuple(centre) == (284, 260)
    assert img[50, 50] == 0


def test_to_distortion_array():
    geom = LPD_MiniGeometry.from_origin()
    # Smoketest
    distortion = geom.to_distortion_array()
    assert isinstance(distortion, np.ndarray)
    assert distortion.shape == (256, 256, 4, 3)

    # Coordinates in m, origin at corner; max x & y should be ~ 50cm
    assert 0.12 < distortion[..., 1].max() < 0.70
    assert 0.12 < distortion[..., 2].max() < 0.70
    assert 0.0 <= distortion[..., 1].min() < 0.01
    assert 0.0 <= distortion[..., 2].min() < 0.01


def test_data_coords_to_positions():
    geom = LPD_MiniGeometry.from_origin()

    module_no = np.zeros(16, dtype=np.int16)
    # Points near the centre of each tile
    slow_scan = np.tile(np.linspace(16, 240, num=8, dtype=np.float32), 2)
    fast_scan = np.array([64, 192], dtype=np.float32).repeat(8)

    tileno, tile_ss, tile_fs = geom._module_coords_to_tile(slow_scan, fast_scan)
    np.testing.assert_allclose(tileno,
                       [7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 10, 11, 12, 13, 14, 15])
    np.testing.assert_allclose(tile_ss, 16)
    np.testing.assert_allclose(tile_fs, 64)

    res = geom.data_coords_to_positions(module_no, slow_scan, fast_scan)

    assert res.shape == (16, 3)

    resx, resy, resz = res.T

    np.testing.assert_allclose(resz, 0)

    assert (np.diff(resy[:8]) > 0).all()  # T1-T8 Monotonically increasing
    assert (np.diff(resy[9:]) > 0).all()  # T9-T16 Monotonically increasing
    assert -0.031 > resx.max() > resx.min() > -0.099


def test_invert_xfel_lpd_geom(tmpdir):
    src_file = pjoin(tests_dir, 'lpd_mar_18.h5')
    dst_file = pjoin(str(tmpdir), 'lpd_inverted.h5')
    invert_xfel_lpd_geom(src_file, dst_file)
    with h5py.File(src_file, 'r') as fsrc, h5py.File(dst_file, 'r') as fdst:
        np.testing.assert_array_equal(
            fsrc['Q1/M1/Position'][:], -1 * fdst['Q1/M1/Position'][:]
        )
        np.testing.assert_array_equal(
            fsrc['Q1/M1/T07/Position'][:], -1 * fdst['Q1/M1/T07/Position'][:]
        )


def test_to_pyfai_detector():
    geom = LPD_MiniGeometry.from_origin()
    agipd_pyfai = geom.to_pyfai_detector()
    assert isinstance(agipd_pyfai, pyFAI.detectors.Detector)
