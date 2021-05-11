from extra_geom import GenericGeometry

import extra_geom
import numpy as np
import pytest

pixel_size = 0.01

simple_config = {'pixel_size': pixel_size,
                 'fast_pixels': 128,
                 'slow_pixels': 64,
                 'corner_coordinates': [pixel_size * np.array([5, 1, 0]),
                                        pixel_size * np.array([-133, 5, 0]),
                                        pixel_size * np.array([2, -133, 0]),
                                        pixel_size * np.array([-131, -130, 0])],
                 'n_tiles_per_module': 2,
                 'fs_vec': np.array([0, 1, 0]),
                 'ss_vec': np.array([1, 0, 0]),
                 'tile_gap': 2 * pixel_size,
                 'tile_vec': [0, -1, 0]
                 }


@pytest.fixture
def simple():
    return GenericGeometry(**simple_config)


def test_simple_geometry(simple):
    """ Tests that an object was created """
    assert type(simple) == extra_geom.detectors.GenericGeometry

