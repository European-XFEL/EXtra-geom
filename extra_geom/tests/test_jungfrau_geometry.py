
from tempfile import TemporaryDirectory

import numpy as np
import pyFAI.detectors
from matplotlib.axes import Axes
from cfelpyutils.geometry import load_crystfel_geometry
from extra_data import RunDirectory
from extra_data.components import JUNGFRAU
from extra_data.tests.make_examples import make_jungfrau_run
from extra_geom import JUNGFRAUGeometry

from .utils import assert_geom_close


def jf4m_geometry():
    x_start, y_start = 1125, 1078
    mod_width = (256 * 4) + (2 * 3)
    mod_height = (256 * 2) + 2

    module_pos = [  #
        (x_start - mod_width, y_start - mod_height - (i * (mod_height + 33)))
        for i in range(4)
    ] + [
        (-x_start, -y_start + (i * (mod_height + 33))) for i in range(4)
    ]
    orientations = [(-1, -1) for _ in range(4)] + [(1, 1) for _ in range(4)]

    return JUNGFRAUGeometry.from_module_positions(module_pos, orientations=orientations)


def test_write_read_crystfel_file(tmpdir):
    geom = jf4m_geometry()
    path = str(tmpdir / 'test.geom')
    geom.write_crystfel_geom(filename=path, photon_energy=9000,
                             adu_per_ev=0.0042, clen=0.101)

    loaded = JUNGFRAUGeometry.from_crystfel_geom(path)
    assert_geom_close(loaded, geom)
    assert loaded.metadata['crystfel']['adu_per_eV'] == 0.0042
    assert loaded.metadata['crystfel']['clen'] == 0.101
    assert loaded.metadata['crystfel']['photon_energy'] == 9000

    # Load the geometry file with cfelpyutils and test the rigid groups
    geom_dict = load_crystfel_geometry(path).detector
    assert geom_dict['panels']['p0a0']['res'] == 1 / 75e-6
    p3a7 = geom_dict['panels']['p3a7']
    assert p3a7['orig_min_ss'] == 256
    assert p3a7['orig_max_ss'] == 511
    assert p3a7['orig_min_fs'] == 768
    assert p3a7['orig_max_fs'] == 1023

    # Check that metadata is written back to .geom file OK
    path2 = str(tmpdir / 'test2.geom')
    loaded.write_crystfel_geom(filename=path2, photon_energy=9000)
    re_loaded = JUNGFRAUGeometry.from_crystfel_geom(path)
    assert re_loaded.metadata == loaded.metadata


def test_get_pixel_positions():
    geom = jf4m_geometry()

    pixelpos = geom.get_pixel_positions()
    assert pixelpos.shape == (8, 512, 1024, 3)
    px = pixelpos[..., 0]
    py = pixelpos[..., 1]

    assert -0.09 < px.min() < -0.07
    assert 0.09 > px.max() > 0.07
    assert -0.09 < py.min() < -0.07
    assert 0.09 > py.max() > 0.07


def test_position_modules_with_labelled_array():
    geom = jf4m_geometry()

    with TemporaryDirectory() as td:
        make_jungfrau_run(td)
        run = RunDirectory(td).select_trains(np.s_[:5])
        jf = JUNGFRAU(run)
        data = jf.get_array('data.adc')
        positioned, centre = geom.position_modules(data)


def test_to_pyfai_detector():
    geom = jf4m_geometry()
    jf4m_pyfai = geom.to_pyfai_detector()
    assert isinstance(jf4m_pyfai, pyFAI.detectors.Detector)
    assert jf4m_pyfai.MAX_SHAPE == (8*512, 1024)


def test_inspect():
    # Smoke test
    geom = JUNGFRAUGeometry.example(n_modules=8)

    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)

def test_rotate_z():
    geom = JUNGFRAUGeometry.example(n_modules=1)
    px_size = geom.pixel_size

    # Start with the slow-scan axis pointing up (+y), fast-scan pointing beam-left (+x)
    np.testing.assert_allclose(geom.modules[0][0].ss_vec / px_size, [0, 1, 0])
    np.testing.assert_allclose(geom.modules[0][0].fs_vec / px_size, [1, 0, 0])

    # Step from ASIC 0 to ASIC 4, vertically above it
    tile_step = geom.modules[0][4].corner_pos - geom.modules[0][0].corner_pos
    np.testing.assert_almost_equal(tile_step[0], 0)
    assert .01 < tile_step[1] < .03  # ASICS are about 2cm
    np.testing.assert_almost_equal(tile_step[2], 0)
    tile_step_magnitude = tile_step[1]

    # Rotate clockwise around the z axis (looking along the axis)
    rotated = geom.rotate((0, 0, 10), degrees=True)

    # Slow scan: [0, 1, 0] -> [-?, 1-?, 0]
    ss_x, ss_y, ss_z = rotated.modules[0][0].ss_vec / px_size
    assert -0.2 < ss_x < 0
    assert 0.8 < ss_y < 1
    np.testing.assert_almost_equal(ss_z, 0)

    # Fast scan: [1, 0, 0] -> [1-?, ?, 0]
    fs_x, fs_y, fs_z = rotated.modules[0][0].fs_vec / px_size
    assert 0.8 < fs_x < 1
    assert 0 < fs_y < 0.2
    np.testing.assert_almost_equal(fs_z, 0)

    # Tile step: [0, 1, 0] -> [-?, 1-?, 0]
    tile_step = rotated.modules[0][4].corner_pos - rotated.modules[0][0].corner_pos
    assert -0.2 < tile_step[0] / tile_step_magnitude < 0
    assert 0.8 < tile_step[1] / tile_step_magnitude < 1
    np.testing.assert_almost_equal(tile_step[2], 0)

    # Rotate anticlockwise (looking along the axis)
    rotated2 = geom.rotate((0, 0, -5), degrees=True)

    # Slow scan: [0, 1, 0] -> [+?, 1-?, 0]
    ss_x, ss_y, ss_z = rotated2.modules[0][0].ss_vec / px_size
    assert 0 < ss_x < 0.2
    assert 0.8 < ss_y < 1
    np.testing.assert_almost_equal(ss_z, 0)

    # Fast scan: [1, 0, 0] -> [-?, 1-?, 0]
    fs_x, fs_y, fs_z = rotated2.modules[0][0].fs_vec / px_size
    assert 0.8 < fs_x < 1
    assert -0.2 < fs_y < 0
    np.testing.assert_almost_equal(fs_z, 0)

    # Tile step: [0, 1, 0] -> [+?, 1-?, 0]
    tile_step = rotated2.modules[0][4].corner_pos - rotated2.modules[0][0].corner_pos
    assert 0 < tile_step[0] / tile_step_magnitude < 0.2
    assert 0.8 < tile_step[1] / tile_step_magnitude < 1
    np.testing.assert_almost_equal(tile_step[2], 0)

def test_rotate_x():
    geom = JUNGFRAUGeometry.example(n_modules=1)
    px_size = geom.pixel_size

    # Start with the slow-scan axis pointing up (+y), fast-scan pointing beam-left (+x)
    np.testing.assert_allclose(geom.modules[0][0].ss_vec / px_size, [0, 1, 0])
    np.testing.assert_allclose(geom.modules[0][0].fs_vec / px_size, [1, 0, 0])

    # Step from ASIC 0 to ASIC 4, vertically above it
    tile_step = geom.modules[0][4].corner_pos - geom.modules[0][0].corner_pos
    np.testing.assert_almost_equal(tile_step[0], 0)
    assert .01 < tile_step[1] < .03  # ASICS are about 2cm
    np.testing.assert_almost_equal(tile_step[2], 0)
    tile_step_magnitude = tile_step[1]

    # Postive rotation around x tilts the top away from the source (+z)
    rotated = geom.rotate((5, 0, 0), degrees=True)

    # Slow scan: [0, 1, 0] -> [0, 1-?, +?]
    ss_x, ss_y, ss_z = rotated.modules[0][0].ss_vec / px_size
    np.testing.assert_almost_equal(ss_x, 0)
    assert 0.8 < ss_y < 1  # y decreases
    assert 0 < ss_z < 0.2  # z increases

    # Fast scan: no change
    np.testing.assert_allclose(geom.modules[0][0].fs_vec / px_size, [1, 0, 0])

    # Tile step: [0, 1, 0] -> [0, 1-?, +?]
    tile_step = rotated.modules[0][4].corner_pos - rotated.modules[0][0].corner_pos
    np.testing.assert_almost_equal(tile_step[0], 0)
    assert 0.8 < tile_step[1] / tile_step_magnitude < 1
    assert 0 < tile_step[2] < 0.2
