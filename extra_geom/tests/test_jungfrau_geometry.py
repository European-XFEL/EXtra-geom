
import numpy as np
from cfelpyutils.crystfel_utils import load_crystfel_geometry

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

    # Load the geometry file with cfelpyutils and test the rigid groups
    geom_dict = load_crystfel_geometry(path)
    assert geom_dict['panels']['p0a0']['res'] == 1 / 75e-6
    p3a7 = geom_dict['panels']['p3a7']
    assert p3a7['min_ss'] == 256
    assert p3a7['max_ss'] == 511
    assert p3a7['min_fs'] == 768
    assert p3a7['max_fs'] == 1023

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
