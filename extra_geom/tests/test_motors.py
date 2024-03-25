import numpy as np
import pytest
from extra_geom import AGIPD_1MGeometry
from extra_geom.crystfel_fmt import motors_to_geom
from extra_geom.motors import AGIPD_1MMotors, read_motors_from_geom

motor_text_reference = """\
;XGEOM MOTORS=4,2
;XGEOM MOTOR_Q1=1,-2
;XGEOM MOTOR_Q2=1,-2
;XGEOM MOTOR_Q3=1,-2
;XGEOM MOTOR_Q4=1,-2
"""


def test_move_geom_by():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    quad_pos_new = [(-530, 635), (-555., -20.), (525., -170.), (547.5, 485.)]
    motor_pos = [[1, -2], [1, -2], [1, -2], [1, -2]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)
    tracker = AGIPD_1MMotors(geom)
    geom2 = tracker.move_geom_by(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)


def test_geom_at():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    quad_pos_new = [(-530, 635), (-555., -20.), (525., -170.), (547.5, 485.)]
    motor_ref = [[0, 0], [0, 0], [0, 0], [0, 0]]
    motor_pos = [[1, -2], [1, -2], [1, -2], [1, -2]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)
    with pytest.raises(ValueError):
        AGIPD_1MMotors.with_reference_positions(geom)

    tracker = AGIPD_1MMotors(geom)
    with pytest.raises(ValueError):
        tracker.geom_at(motor_pos)

    tracker = AGIPD_1MMotors.with_reference_positions(geom, motor_ref)
    geom2 = tracker.geom_at(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)

    geom.motor_positions = np.array(motor_ref)
    tracker = AGIPD_1MMotors.with_reference_positions(geom)
    geom2 = tracker.geom_at(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)


def test_read_write_motors():
    motor_pos = [[1, -2], [1, -2], [1, -2], [1, -2]]

    motor_text = motors_to_geom(motor_pos)
    assert motor_text_reference == motor_text

    motor_pos_read = read_motors_from_geom(motor_text)
    np.testing.assert_allclose(motor_pos, motor_pos_read)
