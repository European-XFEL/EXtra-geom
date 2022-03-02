from os.path import abspath, dirname
from os.path import join as pjoin

import h5py
import numpy as np
import pyFAI.detectors
import pytest
from cfelpyutils.geometry import load_crystfel_geometry
from matplotlib.axes import Axes
from testpath import assert_isfile

from extra_geom import LPD_1MGeometry
from extra_geom.detectors import invert_xfel_lpd_geom

from .utils import assert_geom_close

tests_dir = dirname(abspath(__file__))


def test_write_read_crystfel_file(tmpdir):
    geom = LPD_1MGeometry.from_quad_positions(
        [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    )
    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path)

    with open(path, 'r') as f:
        contents = f.read()
    with open(path, 'w') as f:
        f.write('clen = 0.119\n')
        f.write('adu_per_eV = 0.0075\n')

        f.write(contents)
    # Load the geometry file with cfelpyutils and test the ridget groups
    loaded = LPD_1MGeometry.from_crystfel_geom(path)
    np.testing.assert_allclose(
        loaded.modules[0][0].corner_pos, geom.modules[0][0].corner_pos
    )
    np.testing.assert_allclose(loaded.modules[0][0].fs_vec, geom.modules[0][0].fs_vec)


    geom_dict = load_crystfel_geometry(path).detector
    quad_gr0 = ['p0a0', 'p0a1', 'p0a2', 'p0a3', 'p0a4', 'p0a5', 'p0a6', 'p0a7',
                'p0a8', 'p0a9', 'p0a10','p0a11', 'p0a12', 'p0a13', 'p0a14',
                'p0a15', 'p1a0', 'p1a1','p1a2', 'p1a3','p1a4','p1a5','p1a6',
                'p1a7', 'p1a8', 'p1a9', 'p1a10', 'p1a11', 'p1a12', 'p1a13',
                'p1a14', 'p1a15', 'p2a0', 'p2a1', 'p2a2', 'p2a3', 'p2a4', 'p2a5',
                'p2a6', 'p2a7', 'p2a8', 'p2a9', 'p2a10', 'p2a11', 'p2a12', 'p2a13',
                'p2a14','p2a15', 'p3a0', 'p3a1','p3a2', 'p3a3', 'p3a4', 'p3a5',
                'p3a6', 'p3a7', 'p3a8','p3a9', 'p3a10', 'p3a11', 'p3a12', 'p3a13',
                'p3a14', 'p3a15']
    assert geom_dict['rigid_groups']['p0'] == quad_gr0[:16]
    assert geom_dict['rigid_groups']['p3'] == quad_gr0[-16:]
    assert geom_dict['rigid_groups']['q0'] == quad_gr0


def test_read_write_xfel_file_quadpos(tmpdir):
    quad_pos = [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    geom = LPD_1MGeometry.from_quad_positions(quad_pos)
    path = str(tmpdir / 'lpd_geom.h5')
    quad_pos_out = geom.to_h5_file_and_quad_positions(path)

    np.testing.assert_allclose(quad_pos_out, quad_pos)
    assert_isfile(path)

    loaded = LPD_1MGeometry.from_h5_file_and_quad_positions(path, quad_pos_out)
    assert_geom_close(loaded, geom)


def test_read_write_xfel_file(tmpdir):
    quad_pos = [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    geom = LPD_1MGeometry.from_quad_positions(quad_pos)
    path = str(tmpdir / 'lpd_geom.h5')
    geom.to_h5_file_and_quad_positions(path)

    assert_isfile(path)

    loaded = LPD_1MGeometry.from_h5_file(path)
    assert_geom_close(loaded, geom)


def test_quad_positions_with_file():
    path = pjoin(tests_dir, 'lpd_mar_18.h5')
    quad_pos = [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    geom = LPD_1MGeometry.from_h5_file_and_quad_positions(path, quad_pos)

    np.testing.assert_allclose(geom.quad_positions(path), quad_pos)


def test_quad_positions_no_file():
    quad_pos = [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    geom = LPD_1MGeometry.from_quad_positions(quad_pos)

    np.testing.assert_allclose(geom.quad_positions(), quad_pos)


def test_offset():
    quad_pos = [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    geom = LPD_1MGeometry.from_quad_positions(quad_pos)
    y_orig = np.array([m[0].corner_pos[1] for m in geom.modules])

    # Uniform shift for all modules, all tiles
    all_shifted = geom.offset((0, 1e-3))
    y1 = np.array([m[0].corner_pos[1] for m in all_shifted.modules])
    np.testing.assert_allclose(y1, y_orig + 1e-3)

    # Select some modules
    q4_shifted = geom.offset((0, 2e-3), modules=np.s_[12:])
    y2 = np.array([m[0].corner_pos[1] for m in q4_shifted.modules])
    np.testing.assert_allclose(y2[:12], y_orig[:12])
    np.testing.assert_allclose(y2[12:], y_orig[12:] + 2e-3)

    quad_pos_modified = q4_shifted.quad_positions()
    np.testing.assert_allclose(quad_pos_modified[:3], quad_pos[:3])
    np.testing.assert_allclose(
        quad_pos_modified[3], np.array(quad_pos[3]) + [0, 2]  # quad positions in mm
    )

    # Per-module shift
    q3_shifted = geom.offset(np.repeat([
        (0, 0), (0, 0), (0, 2e-3), (0, 0),
    ], repeats=4, axis=0))
    y3 = np.array([m[0].corner_pos[1] for m in q3_shifted.modules])
    np.testing.assert_allclose(y3[:8], y_orig[:8])
    np.testing.assert_allclose(y3[8:12], y_orig[8:12] + 2e-3)
    np.testing.assert_allclose(y3[12:], y_orig[12:])

    # Per-tile shift
    shift = np.zeros((4, 16, 3), dtype=np.float64)
    shift[:, 5, 1] = 3e-3  # Shift T6 of each module in y dimension
    q2_t2_shifted = geom.offset(shift, modules=np.s_[4:8])
    y_t1 = np.array([m[0].corner_pos[1] for m in q2_t2_shifted.modules])
    np.testing.assert_allclose(y_t1, y_orig)
    y_t6 = np.array([m[5].corner_pos[1] for m in q2_t2_shifted.modules])
    y_t6_orig = np.array([m[5].corner_pos[1] for m in geom.modules])
    np.testing.assert_allclose(y_t6[:4], y_t6_orig[:4])
    np.testing.assert_allclose(y_t6[4:8], y_t6_orig[4:8] + 3e-3)
    np.testing.assert_allclose(y_t6[8:], y_t6_orig[8:])

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
    quad_pos = [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    geom = LPD_1MGeometry.from_quad_positions(quad_pos)

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

    # Select some modules
    q4_rotated = geom.rotate((0, 0, 180), modules=np.s_[12:])
    assert q4_rotated.output_array_for_position_fast().shape == \
        geom.output_array_for_position_fast().shape
    np.testing.assert_array_almost_equal(
        q4_rotated.modules[-1][7].corners(),
        np.roll(geom.modules[-1][15].corners(), 2, 0)
    )

    # Per-module rotation
    rotations = [(0, 0, 90), (0, 0, 180), (90, 0, 0), (0, 90, 0)]
    q1_rotated = geom.rotate(rotations, modules=np.s_[:4])
    mod3_center = np.mean([t.centre() for t in geom.modules[2]], axis=0)
    mod4_center = np.mean([t.centre() for t in geom.modules[3]], axis=0)
    np.testing.assert_allclose([t.corner_pos[1] for t in q1_rotated.modules[2]], mod3_center[1])
    np.testing.assert_allclose([t.corner_pos[0] for t in q1_rotated.modules[3]], mod4_center[0])

    # Per-tile rotation
    rot = np.zeros((4, 16, 3), dtype=np.float64)
    rot[:, 5, :] = (0, 0, 180)  # rotate T6 of each module
    q2_t2_rotated = geom.rotate(rot, modules=np.s_[4:8])
    for mod_rot, mod_ref in zip(q2_t2_rotated.modules[4:8], geom.modules[4:8]):
        np.testing.assert_array_almost_equal(
            mod_rot[5].corners(),
            np.roll(mod_ref[5].corners(), 2, 0)
        )
    for mod_rot, mod_ref in zip(q2_t2_rotated.modules[:4], geom.modules[:4]):
        np.testing.assert_array_almost_equal(
            mod_rot[5].corners(),
            mod_ref[5].corners(),
        )

    # set rotation center, per-det/mod/tile
    rot = np.zeros((4, 16, 3), dtype=np.float64)
    rot[:, 2, :] = (90, 0, 0)
    center = np.zeros((4, 16, 3), dtype=np.float64)
    center[0, 2, :] = geom.modules[8][2].centre()
    center[1, 2, :] = geom.modules[9][2].centre()
    center[2, 2, :] = geom.modules[10][2].centre()
    center[3, 2, :] = geom.modules[11][2].centre()
    q3_t2_rotated = geom.rotate(rot, center=center, modules=np.s_[8:12])
    y = np.array([m[2].corner_pos[1] for m in q3_t2_rotated.modules[8:12]])
    np.testing.assert_array_almost_equal(
        y,
        [
            geom.modules[8][2].centre()[1],
            geom.modules[9][2].centre()[1],
            geom.modules[10][2].centre()[1],
            geom.modules[11][2].centre()[1]
        ]
    )

    # Wrong number of modules
    with pytest.raises(ValueError):
        geom.rotate(np.ones((15, 3)))

    # Rotation for 16 modules, but only 4 selected
    with pytest.raises(ValueError):
        geom.rotate(np.ones((16, 16, 3)), modules=np.s_[:4])

    # angles must be 3D
    with pytest.raises(ValueError):
        geom.rotate((1, 1))

    # angles in radian
    deg = geom.rotate((90, 0, 0), modules=np.s_[:1], tiles=np.s_[:1])
    rad = geom.rotate((np.pi / 2, 0, 0), modules=np.s_[:1], tiles=np.s_[:1], degrees=False)
    np.testing.assert_array_almost_equal(
        deg.modules[1][1].corners(),
        rad.modules[1][1].corners()
    )


def test_inspect():
    geom = LPD_1MGeometry.from_quad_positions(
        [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    )
    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)


def test_snap_assemble_data():
    geom = LPD_1MGeometry.from_quad_positions(
        [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    )

    stacked_data = np.zeros((16, 256, 256))
    img, centre = geom.position_modules_fast(stacked_data)
    assert img.shape == (1202, 1104)
    assert tuple(centre) == (604, 547)
    assert np.isnan(img[0, 0])
    assert img[50, 50] == 0

def test_to_distortion_array():
    geom = LPD_1MGeometry.from_quad_positions(
        [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    )
    # Smoketest
    distortion = geom.to_distortion_array()
    assert isinstance(distortion, np.ndarray)
    assert distortion.shape == (4096, 256, 4, 3)

    # Coordinates in m, origin at corner; max x & y should be ~ 50cm
    assert 0.40 < distortion[..., 1].max() < 0.70
    assert 0.40 < distortion[..., 2].max() < 0.70
    assert 0.0 <= distortion[..., 1].min() < 0.01
    assert 0.0 <= distortion[..., 2].min() < 0.01

def test_data_coords_to_positions():
    geom = LPD_1MGeometry.from_quad_positions(
        [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    )

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
    assert (np.diff(resy[8:]) > 0).all()  # T9-T16 Monotonically increasing
    assert -0.128 > resx.max() > resx.min() > -0.280

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
    geom = LPD_1MGeometry.from_quad_positions(
        [(11.4, 299), (-11.5, 8), (254.5, -16), (278.5, 275)]
    )
    agipd_pyfai = geom.to_pyfai_detector()
    assert isinstance(agipd_pyfai, pyFAI.detectors.Detector)
