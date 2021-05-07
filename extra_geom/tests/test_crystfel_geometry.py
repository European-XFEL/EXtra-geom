from extra_geom import CrystFEL_Geometry, AGIPD_1MGeometry
from cfelpyutils.crystfel_utils import load_crystfel_geometry

import extra_geom
import pytest

agipd_geom_file = 'agipd_simple_2d.geom'
jungfrau_geom_file = 'jungfrau.geom'

simple_config = {'pixel_size': 0.001,
                 'fast_pixels': 128,
                 'slow_pixels': 64,}

@pytest.fixture
def simple():
    return CrystFEL_Geometry(**simple_config)

@pytest.fixture
def geom():
    return CrystFEL_Geometry.from_crystfel_geom(agipd_geom_file)

@pytest.fixture
def agipd():
    return AGIPD_1MGeometry.from_crystfel_geom(agipd_geom_file)

@pytest.fixture()
def geom_dict():
    return load_crystfel_geometry(agipd_geom_file)


def test_simple_geometry(simple):
    """ Tests that an object was created """
    assert type(simple) == extra_geom.detectors.CrystFEL_Geometry


def test_from_crystfel_geometry():
    """ Tests that an object was created """
    assert type(geom) == extra_geom.detectors.CrystFEL_Geometry


def test_content(geom, agipd):
    """ Test that the mandatory properties are populated:

    detector_type_name
    pixel_size
    frag_ss_pixels
    frag_fs_pixels
    n_quads
    n_modules
    n_tiles_per_module
    expected_data_shap
    _pixel_corners = np.array([  # pixel units; overridden for DSSC
        [0, 1, 1, 0],  # slow-scan
        [0, 0, 1, 1]   # fast-scan
    ])
    _draw_first_px_on_tile = 1
    """
    assert geom.pixel_size == agipd.pixel_size
    assert geom.frag_fs_pixels == agipd.frag_fs_pixels
    assert geom.frag_ss_pixels == agipd.frag_ss_pixels

    assert geom.n_quads == agipd.n_quads
    assert geom.n_modules == agipd.n_modules


def test_detect_shapes(geom_dict):
    CrystFEL_Geometry.detect_shapes(geom_dict)
