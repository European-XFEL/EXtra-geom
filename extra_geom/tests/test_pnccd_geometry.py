import numpy as np
import pytest
from matplotlib.axes import Axes

from extra_geom import PNCCDGeometry

from .utils import assert_geom_close


def test_positions():
    # Positional constants
    gap_pixels = 50
    offset_pixels = 10

    # Calculate parameters for the relative/absolute methods
    top_offset = (-offset_pixels * PNCCDGeometry.pixel_size, 0, 0)
    bottom_offset = (offset_pixels * PNCCDGeometry.pixel_size, 0, 0)
    absolute_pos = np.array(
        [(-PNCCDGeometry.frag_fs_pixels // 2 - offset_pixels,
          PNCCDGeometry.frag_ss_pixels + gap_pixels // 2,
          0),
         (-PNCCDGeometry.frag_fs_pixels // 2 + offset_pixels,
          -gap_pixels // 2,
          0)]
    ) * PNCCDGeometry.pixel_size

    # Create the geometries
    absolute_geom = PNCCDGeometry.from_absolute_positions(top=absolute_pos[0], bottom=absolute_pos[1])
    relative_geom = PNCCDGeometry.from_relative_positions(gap=gap_pixels * PNCCDGeometry.pixel_size,
                                                          top_offset=top_offset,
                                                          bottom_offset=bottom_offset)

    # They should be the same
    assert_geom_close(absolute_geom, relative_geom)

# No re-shaping is required for the first shape, the class should automatically
# reshape the data with the second shape.
@pytest.mark.parametrize("shape", [PNCCDGeometry.expected_data_shape,
                                   (PNCCDGeometry.frag_ss_pixels * PNCCDGeometry.n_modules, PNCCDGeometry.frag_fs_pixels)])
def test_snap_assemble_data(shape):
    gap_pixels = 100
    geom = PNCCDGeometry.from_relative_positions(gap=gap_pixels * PNCCDGeometry.pixel_size)

    stacked_data = np.zeros(shape)
    img, centre = geom.position_modules_fast(stacked_data)
    expected_img_shape = (PNCCDGeometry.frag_ss_pixels * PNCCDGeometry.n_modules + gap_pixels,
                          PNCCDGeometry.frag_fs_pixels)

    assert img.shape == expected_img_shape
    assert tuple(centre) == (expected_img_shape[0] // 2, expected_img_shape[1] // 2)
    # Everywhere in the gap should be NaN
    assert np.isnan(img[img.shape[0] // 2, 0])
    assert img[50, 50] == 0

def test_inspect():
    geom = PNCCDGeometry.from_relative_positions()

    # Smoketest
    ax = geom.inspect()
    assert isinstance(ax, Axes)
