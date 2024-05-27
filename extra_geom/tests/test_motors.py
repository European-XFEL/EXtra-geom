import numpy as np
import pytest
from extra_geom import AGIPD_1MGeometry
from extra_geom.crystfel_fmt import motors_to_geom
from extra_geom.motors import (AGIPD_1MMotors, BaseMotorTracker,
                               read_motors_from_geom)

motor_text_template = """\
;XGEOM MOTORS={num_groups},{num_motors}
;XGEOM MOTOR_Q1={q1m1},{q1m2}
;XGEOM MOTOR_Q2={q2m1},{q2m2}
;XGEOM MOTOR_Q3={q3m1},{q3m2}
;XGEOM MOTOR_Q4={q4m1},{q4m2}
"""


class InvalidMotorTracker(BaseMotorTracker):
    groups = [
        np.s_[0:4], np.s_[4:8],
    ]
    default_motor_axes = np.array([
        [[0, -1], [-1, 0]],  # Q1
    ])


def test_move_geom_by():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    quad_pos_new = [(-530, 635), (-555., -20.), (525., -170.), (547.5, 485.)]
    motor_pos = [[-2, 1], [-2, 1], [-2, 1], [-2, 1]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)
    tracker = AGIPD_1MMotors(geom)
    geom2 = tracker.move_geom_by(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)

    with pytest.raises(ValueError):
        tracker.move_geom_by([[1, -2]] * 3)


def test_geom_at():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    quad_pos_new = [(-530, 635), (-555., -20.), (525., -170.), (547.5, 485.)]
    motor_ref = [[0, 0], [0, 0], [0, 0], [0, 0]]
    motor_pos = [[-2, 1], [-2, 1], [-2, 1], [-2, 1]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)

    tracker = AGIPD_1MMotors(geom)
    with pytest.raises(ValueError):
        tracker.geom_at(motor_pos)

    tracker = AGIPD_1MMotors(geom, motor_ref)
    geom2 = tracker.geom_at(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)

    geom.motor_positions = np.array(motor_ref)
    tracker = AGIPD_1MMotors(geom)
    geom2 = tracker.geom_at(motor_pos)

    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)

    with pytest.raises(ValueError):
        tracker.geom_at([[1, -2]] * 3)


def test_read_write_motors():
    motor_pos = [[-2, 1], [-2, 1], [-2, 1], [-2, 1]]
    motor_param = {
        f"q{i // 2 + 1}m{i % 2 + 1}": motor_pos[i // 2][i % 2]
        for i in range(8)
    }
    motor_param.update({
        "num_groups": len(motor_pos),
        "num_motors": len(motor_pos[0]),
    })
    motor_text_reference = motor_text_template.format(**motor_param)

    motor_text = motors_to_geom(motor_pos)
    assert motor_text_reference == motor_text

    motor_pos_read = read_motors_from_geom(motor_text)
    np.testing.assert_allclose(motor_pos, motor_pos_read)

    motor_invalid_param = motor_param.copy()
    motor_invalid_param["num_groups"] = 'a'
    motor_text = motor_text_template.format(**motor_invalid_param)
    with pytest.raises(ValueError):
        read_motors_from_geom(motor_text)

    motor_invalid_param = motor_param.copy()
    motor_invalid_param["q2m2"] = 'a'
    motor_text = motor_text_template.format(**motor_invalid_param)
    with pytest.raises(ValueError):
        read_motors_from_geom(motor_text)

    motor_invalid_param = motor_param.copy()
    motor_invalid_param["q2m2"] = "0,0"
    motor_text = motor_text_template.format(**motor_invalid_param)
    with pytest.raises(ValueError):
        read_motors_from_geom(motor_text)

    motor_text = motor_text_template.format(**motor_param)
    motor_text = '\n'.join(motor_text.split('\n')[:-2] + [""])
    with pytest.raises(ValueError):
        read_motors_from_geom(motor_text)


def test_other_methods():
    quad_pos = [(-525, 625), (-550, -10), (520, -160), (542.5, 475)]
    quad_pos_new = [(-520, 615), (-545., 0.), (515., -150.), (537.5, 465.)]
    motor_pos = [[-2, 1], [-2, 1], [-2, 1], [-2, 1]]

    geom = AGIPD_1MGeometry.from_quad_positions(quad_pos)
    with pytest.raises(ValueError):
        InvalidMotorTracker(geom)

    with pytest.raises(ValueError):
        AGIPD_1MMotors(geom, [[1, -2]] * 3)

    tracker = AGIPD_1MMotors(geom)
    with pytest.raises(ValueError):
        tracker.with_motor_axes(InvalidMotorTracker.default_motor_axes)

    tracker2 = tracker.with_motor_axes(AGIPD_1MMotors.default_motor_axes * -1)
    geom2 = tracker2.move_geom_by(motor_pos)
    np.testing.assert_allclose(geom2.quad_positions(), quad_pos_new)
