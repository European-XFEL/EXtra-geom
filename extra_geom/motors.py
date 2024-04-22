import numpy as np


def read_motors_from_geom(text):
    """Reads the motor positions from the text format."""
    if isinstance(text, str):
        import io
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


class BaseMotorTracker:
    """Detector motor tracker updates geometry according to motor positions.
    """
    def __init__(self, ref_geom):
        """
        Parameters
        ----------
        ref_geom: one of `extra_geom.DetectorGeometryBase` implementation
            Geometry
        ref_motor_positions: array or sequence
            Reference motor positions
        """
        self.motor_axes = self.default_motor_axes
        self.num_groups, _, self.num_motors = self.motor_axes.shape
        if self.num_groups != len(self.groups):
            raise ValueError(
                "The len of groups does not match the len of axes")

        self.motor_position_shape = (self.num_groups, self.num_motors)
        self.motor_axes_shape = self.motor_axes.shape

        self.ref_geom = ref_geom
        self.ref_motor_positions = None

    @classmethod
    def with_reference_positions(cls, ref_geom, ref_motor_positions=None):
        tracker = cls(ref_geom)
        if ref_motor_positions is None:
            if hasattr(tracker.ref_geom, "motor_positions"):
                ref_motor_positions = tracker.ref_geom.motor_positions
            else:
                raise ValueError(
                    "There is no motor positions in reference geometry")

        tracker.set_reference_positions(ref_motor_positions)
        return tracker

    def set_reference_positions(self, ref_motor_positions):
        """Set reference motor positions.

        Parameters
        ----------
        ref_motor_positions: array or list
            New reference motor positions. The positions are expected
            as array or list of the number of movable groups (quadrants)
            by the number of motors per group.
        """
        ref_motor_positions = np.array(ref_motor_positions, copy=True)
        if ref_motor_positions.shape != self.motor_position_shape:
            raise ValueError(f"Expects array{self.motor_position_shape}: "
                             f"{self.num_groups} groups moving by "
                             f"{self.num_motors} motor each.")

        self.ref_motor_positions = np.array(ref_motor_positions)

    def set_motor_axes(self, new_motor_axes):
        """Set the matrices of transformation motor positions in the positions
        of detector panels.

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
            New matrices of motor axes (transformation matrices). The matrices
            are expected as three dimention array of the number of movable
            groups (quadrants) by the number of panel coordinates (two) by
            the number of motors per group.
        """
        new_motor_axes = np.array(new_motor_axes, copy=True)
        if new_motor_axes.shape != self.motor_axes_shape:
            raise ValueError(f"Expects array{self.motor_axes_shape}: "
                             f"{self.num_groups} groups moving by "
                             f"{self.num_motors} motor each.")
        self.motor_axes = new_motor_axes

    def geom_at(self, motor_positions):
        """Get geometry for absolute motor positions

        This returns a new geometry according to the given
        motor positions with respect to the reference motor positions.

        Parameters
        ----------
        motor_positions: array or list
            New motor positions as array of the number of movable groups
            (quadrants) by the number of motor per group.

        Returns
        -------
        geometry: the same class `ref_geom`
            a new geometry
        """
        if self.ref_motor_positions is None:
            raise ValueError(
                "Define `ref_motor_position` to use this method")

        motor_positions = np.array(motor_positions, copy=True)

        if motor_positions.shape != self.motor_position_shape:
            raise ValueError(f"Expects array{self.motor_position_shape}: "
                             f"{self.num_groups} groups moving by "
                             f"{self.num_motors} motor each.")

        motor_diff = motor_positions - self.ref_motor_positions
        new_geom = self._move_by(1e-3 * motor_diff)
        new_geom.motor_positions = motor_positions
        return new_geom

    def move_geom_by(self, motor_diff):
        """Get geometry for changes of motor positions.

        This returns a new geometry according to the relative changes
        of motor positions from the reference geometry.

        Parameters
        ----------
        motor_diff: array or list
            Changes of motor positions as array of the number of movable groups
            (quadrants) by the number of motor per group.

        Returns
        -------
        geometry: the same class `ref_geom`
            a new geometry
        """
        motor_diff = np.array(motor_diff)
        if motor_diff.shape != self.motor_position_shape:
            raise ValueError(f"Expects array{self.motor_position_shape}: "
                             f"{self.num_groups} groups moving by "
                             f"{self.num_motors} motor each.")
        return self._move_by(1e-3 * motor_diff)

    def _move_by(self, motor_diff):
        """Moves reference geometry relative current position."""
        new_geom = self.ref_geom.offset((0, 0))
        for i in range(self.num_groups):
            det_diff = self.motor_axes[i] @ motor_diff[i]
            new_geom = new_geom.offset(
                det_diff, modules=self.groups[i])
        return new_geom


class AGIPD_1MMotors(BaseMotorTracker):
    # groups of modules driven by motors together
    # Q1, Q2, Q3, Q4
    groups = [
        np.s_[0:4], np.s_[4:8], np.s_[8:12], np.s_[12:16],
    ]

    # transformation matrix (h,v) -> (x,y), where
    #    (h, v) - local motor coordinates
    #    (x, y) - laboratory cooridnates (looking downstream)
    #  | hx vx | |h|
    #  | hy vy | |v|
    default_motor_axes = np.array([
        [[-1, 0], [0, -1]],  # Q1
        [[-1, 0], [0, +1]],  # Q2
        [[+1, 0], [0, +1]],  # Q3
        [[+1, 0], [0, -1]],  # Q4
    ])
