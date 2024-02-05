import numpy as np
from extra_geom import AGIPD_1MGeometry
from extra_geom.motors import read_motors_from_geom

motor_text_reference = """\
;XGEOM MOTORS=4,2
;XGEOM MOTOR_Q1=1,-2
;XGEOM MOTOR_Q2=1,-2
;XGEOM MOTOR_Q3=1,-2
;XGEOM MOTOR_Q4=1,-2
"""


def test_move_by_motors():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    quad_pos_new = [(-530, 635), (-555., -20.), (525., -170.), (547.5, 485.)]
    motor_pos = [[1, -2], [1, -2], [1, -2], [1, -2]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)
    geom2 = geom.move_by_motors(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)


def test_read_write_motors():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    motor_pos = [[1, -2], [1, -2], [1, -2], [1, -2]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)
    geom2 = geom.move_by_motors(motor_pos)

    motor_text = geom2.motors_to_geom()
    assert motor_text_reference == motor_text

    motor_pos_read = read_motors_from_geom(motor_text)
    np.testing.assert_allclose(motor_pos, motor_pos_read)

    geom.set_motor_positions(motor_pos_read)
    np.testing.assert_allclose(geom.motor_positions, motor_pos)
