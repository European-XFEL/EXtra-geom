import io

import numpy as np


def read_motors_from_data(dc, n_groups, n_motors, position_key, atol=0.001):
    """ Read motors from experiment data.

    The function reads motor position of ``n_motors`` motors for
    ``n_groups`` movable groups (e.g. quadrants) from experimental
    data in EXDF format.

    ::
        # open run
        run = open_run(propno, runno)

        # create function returning (source, key)
        # to read motor `m` for group `q`
        position_key = lambda q, m: (
            f"SPB_IRU_AGIPD1M/MOTOR/Q{q+1}M{m+1}",
            "actualPosition"
        )

        # read motor positions
        motors = read_motors_from_data(run, 4, 2, position_key)

        # add aliases
        run2 = run.with_aliases({
            f"motor_q{q+1}m{m+1}": (
                f"SPB_IRU_AGIPD1M/MOTOR/Q{q+1}M{m+1}",
                "actualPosition"
            )
            for q, m in product(range(4), range(2))
        })

        # read motor positions using aliases
        motors = read_motors_from_data(
            run2.alias, 4, 2, lambda q, m: f"motor_q{q+1}m{m+1}")

    Parameters
    ----------
    dc: extra_data.DataCollection or extra_data.AliasIndexer
      Experimental data
    n_groups: int
      The number of movable groups
    n_motors: int
      The number of motors per group
    position_key: callable
      A function `position_key(q, m)` returning an indentificator of
      motor position property for motor `m` in group `q`. If the `dc`
      is DataCollection, then the identificator is a tuple of source
      and key. If the `dc` is AliasIndexer, then - alias

    Returns
    -------
    numpy.ndarray:
      an array with shape (n_groups, n_motors)
    """
    positions = [
        [
            dc[position_key(q, m)].as_single_value(atol=atol)
            for m in range(n_motors)
        ] for q in range(n_groups)
    ]
    return np.array(positions)


def motors_to_geom(positions):
    """Prints the motor positions in the text format."""
    n_groups, n_motors = positions.shape
    meta_lines = [f";XGEOM MOTORS={n_groups},{n_motors}"]
    meta_lines += [
        f";XGEOM MOTOR_Q{q+1}=" + ",".join(
            (str(positions[q, m]) for m in range(n_motors))
        ) for q in range(n_groups)
    ]
    return "\n".join(meta_lines) + "\n"


def read_motors_from_geom(text):
    """Reads the motor positions from the text format."""
    if isinstance(text, str):
        file = io.StringIO(text)
    else:
        file = text

    meta = {
        line[0]: line[2]
        for line in (
            line[7:].partition('=')
            for line in file if line.startswith(";XGEOM ")
        )
    }
    try:
        n_groups, n_motors = (int(n) for n in meta["MOTORS"].split(','))
    except KeyError:
        raise ValueError("MOTORS record is not found")
    except ValueError:
        raise ValueError(
            "Invalid MOTORS format, expected two comma-separated integers")

    positions = []
    for q in range(n_groups):
        try:
            key = f"MOTOR_Q{q+1}"
            pos = [float(n) for n in meta[key].split(',')]
        except KeyError:
            raise ValueError(key + " record is not found")
        except ValueError:
            raise ValueError(
                f"Invalid {key} format, expected {n_motors} "
                "comma-separated floats")
        if len(pos) != n_motors:
            raise ValueError(
                f"Wrong length of {key}, expected {n_motors} floats")

        positions.append(pos)

    return np.array(positions)


class MotorMixin:
    n_movable_groups = 4
    n_motor_per_group = 2
    motor_position_shape = (4, 2)
    motor_axis_shape = (4, 2, 2)

    # groups of modules driven by motors together
    # Q1, Q2, Q3, Q4
    movable_groups = [
        np.s_[0:4], np.s_[4:8], np.s_[8:12], np.s_[12:16],
    ]

    # transformation matrix (h,v) -> (x,y), where
    #    (h, v) - local motor coordinates
    #    (x, y) - laboratory cooridnates (looking downstream)
    #  | hx vx | |h|
    #  | hy vy | |v|
    motor_axes = np.array([
        [[-1, 0], [0, -1]],  # Q1
        [[-1, 0], [0, +1]],  # Q2
        [[+1, 0], [0, +1]],  # Q3
        [[+1, 0], [0, -1]],  # Q4
    ])

    # motor positions in local motor coordinates (h, v)
    # for each movable group of modules
    # [[Q1M1, Q1M2], ..., [Q4M1, Q4M2]]
    # motor_positions = np.array([
    #    [0, 0], [0, 0], [0, 0], [0, 0],
    # ])

    def __init__(self, modules, filename='No file', metadata=None):
        super().__init__(modules, filename, metadata)

        self.motor_axes_shape = (
            self.n_movable_groups,
            2,
            self.n_motor_per_group
        )
        self.motor_position_shape = (
            self.n_movable_groups,
            self.n_motor_per_group
        )
        try:
            with open(filename) as f:
                self.motor_positions = read_motors_from_geom(f)
        except (ValueError, FileNotFoundError):
            pass

    def set_motor_positions(self, new_motor_positions):
        """Set motor positions for the geometry.

        Parameters
        ----------
        new_motor_positions: array or list
          New motor positions as array of the number of movable groups
          (quadrants) by the number of motor per group.
        """
        new_motor_positions = np.array(new_motor_positions, copy=True)
        if new_motor_positions.shape != self.motor_position_shape:
            raise ValueError(f"Expects array{self.motor_position_shape}: "
                             f"{self.n_movable_groups} groups moving by "
                             f"{self.n_motor_per_group} motor each.")
        self.motor_positions = new_motor_positions

    def set_motor_axes(self, new_motor_axes):
        """Set the matrices of transformation motor positions in
        the positions of detector panels.

        ::
            (h, v) - local motor coordinates
            (x, y) - laboratory cooridnates (looking downstream)
            (hx, hy) - the axis of horizontal motor in laboratory coordinates
            (vx, vy) - the axis of vertical motor in laboratory coordinates

            |x| _ | hx vx | |h|
            |y| ‾ | hy vy | |v|

        Parameters
        ----------
        new_motor_axes: array or list
          New matrices of motor axes (transmation matrices). The matrices
          are expected as three dimention array of the number of movable
          groups (quadrants) by the number of panel coordinates (two) by
          the number of motors per group.
        """
        new_motor_axes = np.array(new_motor_axes, copy=True)
        if new_motor_axes.shape != self.motor_axes_shape:
            raise ValueError(f"Expects array{self.motor_axes_shape}: "
                             f"{self.n_movable_groups} groups moving by "
                             f"{self.n_motor_per_group} motor each.")
        self.motor_axes = new_motor_axes

    def move_by_motors(self, new_motor_positions):
        """Move the geometry according to the given motor positions.

        This changes the geometry according to the given motor positions
        with respect the current motor position. If the geometry does not
        have current motor positions, then this assumes that all motors
        are in zero positions.

        Parameters
        ----------
        new_motor_positions: array or list
          New motor positions as array of the number of movable groups
          (quadrants) by the number of motor per group.

        Returns
        -------
        geometry: the same class as self
          a new geometry
        """
        new_motor_positions = np.array(new_motor_positions, copy=True)
        if new_motor_positions.shape != self.motor_position_shape:
            raise ValueError(f"Expects array{self.motor_position_shape}: "
                             f"{self.n_movable_groups} groups moving by "
                             f"{self.n_motor_per_group} motor each.")

        new_geom = self.offset((0, 0))
        if hasattr(self, "motor_positions"):
            motor_diff = (new_motor_positions - self.motor_positions) * 1e-3
        else:
            motor_diff = new_motor_positions * 1e-3

        for i in range(self.n_movable_groups):
            det_diff = self.motor_axes[i] @ motor_diff[i]
            new_geom = new_geom.offset(
                det_diff, modules=self.movable_groups[i])

        new_geom.motor_positions = new_motor_positions
        return new_geom

    def motors_to_geom(self):
        """Format the current motor position as text ready to store in Crystfel
        geometry file.

        Returns
        -------
        str:
          text with motor positions.
        """
        if hasattr(self, "motor_positions"):
            return motors_to_geom(self.motor_positions)
        else:
            return ""
