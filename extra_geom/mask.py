"""Convert between masks as arrays and as a set of rectangular regions
"""

import os
import re
import warnings

from cfelpyutils.crystfel_utils import load_crystfel_geometry
import numpy as np

from . import crystfel_fmt


def delta_method(mask):
    """
    Generalized delta-method algorithm to decompose two-dimensional mask
    (boolean numpy array) into rectangles.

    Based on paper:
    http://library.utia.cas.cz/separaty/2012/ZOI
    /suk-rectangular%20decomposition%20of%20binary%20images.pdf

    Args:
        mask (2D np.array, dtype=bool): input mask.

    Raises:
        ValueError: mask is not a 2D boolean numpy array.

    Returns:
        list: list of rectangles that represent masked regionsin in the form:
            [((x_min, x_max),(y_min, y_max)), ...]
    """

    # Check input:
    if (mask.ndim != 2
            or mask.dtype != bool):
        raise ValueError(f"Expected input - 2D boolean numpy array.")

    A = np.copy(mask)
    res = []

    i_y = 0
    while not np.array_equiv(A, False):

        assert i_y < A.shape[0], (
            "Decomposition GDM: algorithm did not convert the whole matrix.")

        x_min = x_max = y_min = y_max = -1
        for i_x in range(A.shape[1]):
            if x_min < 0 and A[i_y, i_x]:
                x_min = i_x
                y_min = i_y
                if i_x == A.shape[1] - 1:
                    x_max = A.shape[1]
            elif x_min >= 0 and not A[i_y, i_x]:
                x_max = i_x
            elif x_min >= 0 and i_x == A.shape[1] - 1:
                x_max = A.shape[1]

            if x_max >= 0:
                sx = slice(x_min, x_max)
                for i_y2 in range(i_y + 1, A.shape[0]):
                    if not np.array_equal(A[i_y2, sx], A[i_y, sx]):
                        y_max = i_y2
                        break
                else:
                    y_max = A.shape[0]
                res.append(((x_min, x_max), (y_min, y_max)))
                A[y_min:y_max, sx] = False
                x_min = x_max = y_min = y_max = -1

        i_y += 1

    return res

class RegionRect:
    """One masked region within a 3D array of data (modules stacked)

    This can apply to a single module (modno=5) or all modules (modno=None).
    """
    def __init__(self, modno, start_ss, stop_ss, start_fs, stop_fs):
        self.modno = modno  # None -> applies to all modules
        self.start_ss = start_ss
        self.stop_ss = stop_ss
        self.start_fs = start_fs
        self.stop_fs = stop_fs

    def _tuple(self):
        return self.modno, self.start_ss, self.stop_ss, self.start_fs, self.stop_fs

    def __repr__(self):
        return f'RegionRect{self._tuple()}'

    def __hash__(self):
        return hash(self._tuple())

    def __eq__(self, other):
        return isinstance(other, RegionRect) and self._tuple() == other._tuple()

    @property
    def array_slice(self):
        """Get a tuple to use for slicing an array"""
        mod = np.s_[:] if (self.modno is None) else self.modno
        return np.s_[mod, self.start_ss:self.stop_ss, self.start_fs:self.stop_fs]

    def intersection(self, other: 'RegionRect'):
        """Find the intersection of this and another RegionRect

        Returns (overlap, intersection). The first value is True if the two
        regions overlap, and the second is a RegionRect for their intersection.
        """
        modno = self.modno if (self.modno is not None) else other.modno
        if other.modno in (modno, None):
            # On the same module, or one of self/other is for all modules
            r = RegionRect(
                modno,
                max(self.start_ss, other.start_ss),
                min(self.stop_ss, other.stop_ss),
                max(self.start_fs, other.start_fs),
                min(self.stop_fs, other.stop_fs),
            )
            overlap = (r.start_ss < r.stop_ss) and (r.start_fs < r.stop_fs)
            return overlap, r

        return False, RegionRect(0, 0, 0, 0, 0)


class MaskRegions:
    """A set of rectangular regions masked in a stack of module images"""
    def __init__(self, shape, regions=()):
        assert len(shape) == 3
        self.shape = shape
        self.regions = regions  # [RegionRect]

    def __repr__(self):
        npx = self.to_mask_array().sum()
        return (f"<MaskRegions for {self.shape} array: "
                f"{len(self.regions)} covering {npx} pixels>")

    @classmethod
    def from_crystfel_geom(cls, filename):
        """Read mask from a CrystFEL format geometry file.
        """
        geom_dict = load_crystfel_geometry(filename)
        arr_shape = crystfel_fmt.data_shape(geom_dict['panels'])
        if len(arr_shape) == 2:
            arr_shape = (1,) + arr_shape

        regions = []
        for name, info in geom_dict['bad'].items():
            if not info['is_fsss']:
                # Regions can also be defined with x & y coordinates, but these
                # are not handled yet.
                continue

            min_ss, max_ss = info['min_ss'], info['max_ss']
            min_fs, max_fs = info['min_fs'], info['max_fs']
            if (max_ss < min_ss) or (max_fs < min_fs):
                warnings.warn(
                    f"Geometry file {filename} - incomplete information for {name}."
                )
            if 'panel' in info:
                panel_info = geom_dict['panels'][info['panel']]
                modno = crystfel_fmt.panel_modno(panel_info, info['panel'])
                modno = modno or 0  # None -> 0

                # Clip the region to the panel limits
                min_ss = max(min_ss, panel_info['min_ss'])
                max_ss = min(max_ss, panel_info['max_ss'])
                min_fs = max(min_fs, panel_info['min_fs'])
                max_fs = min(max_fs, panel_info['max_fs'])
            else:
                modno = None  # Applies to all panels

            regions.append(RegionRect(
                modno, min_ss, max_ss + 1, min_fs, max_fs + 1,
            ))

        return cls(arr_shape, regions)

    @classmethod
    def from_mask_array(cls, arr):
        """Convert mask from a 2D or 3D boolean numpy array to rectangular regions.

        The masked/selected regions will be where the array contains True.
        """
        if arr.ndim not in {2, 3}:
            raise TypeError(f"Array must be 2D or 3D (got {arr.ndim} dims)")

        def find_regions(arr_2d, modno):
            return [RegionRect(
                modno, min_ss, max_ss, min_fs, max_fs
            ) for ((min_fs, max_fs), (min_ss, max_ss)) in delta_method(arr_2d)]

        if arr.ndim == 2:
            return cls((1,) + arr.shape, find_regions(arr, modno=None))

        # 3D array (modno, slow_scan, fast_scan)

        # First check for regions to be excluded in all panels:
        panel_all = np.logical_and.reduce(arr, axis=0)
        regions = find_regions(panel_all, modno=None)
        is_panel_all_empty = not np.any(panel_all)

        # Loop over all panels:
        for i, panel_arr in enumerate(arr):
            if not is_panel_all_empty:
                panel_arr = np.copy(panel_arr) & ~panel_all

            regions.extend(find_regions(panel_arr, modno=i))

        return cls(arr.shape, regions)

    def to_mask_array(self):
        """Convert the mask rectangles to a 3D boolean numpy array.

        The selected regions will be True in the array.
        """
        res_mask = np.zeros(self.shape, dtype=np.bool_)

        for region in self.regions:
            res_mask[region.array_slice] = True

        return res_mask

    def make_crystfel_bad_regions(self, panels_dict):
        modno_to_panels = {}
        for pname, pinfo in panels_dict.items():
            modno = crystfel_fmt.panel_modno(pinfo, pname)
            modno_to_panels.setdefault(modno, []).append((pname, pinfo))

        def to_dict(region: RegionRect):
            return {
                'min_ss': region.start_ss, 'max_ss': region.stop_ss - 1,
                'min_fs': region.start_fs, 'max_fs': region.stop_fs - 1,
            }

        res = []
        for mask_region in self.regions:
            if mask_region.modno is None:
                res.append(to_dict(mask_region))  # Mask for all panels

            else:
                for pname, pinfo in modno_to_panels[mask_region.modno]:
                    panel_region = RegionRect(
                        mask_region.modno,
                        pinfo['min_ss'], pinfo['max_ss'] + 1,
                        pinfo['min_fs'], pinfo['max_fs'] + 1,
                    )
                    overlap, matching = mask_region.intersection(panel_region)
                    if overlap:
                        d = to_dict(matching)
                        d['panel'] = pname
                        res.append(d)

        return res

    def write_crystfel_geom(self, filename, write_mode):
        """Write the mask regions to a CrystFEL format geometry file.

        *filename* needs to exist already; it will be read to identify the
        panels which the mask applies to.
        """

        text_before = text_after = []
        n_area_start = 0

        # Store and process content of the existing geometry file
        if os.path.exists(filename):
            with open(filename, 'r') as f_geo:
                contents = f_geo.readlines()

            idx_write = len(contents)
            for i, line in enumerate(contents):

                # Search for the mask information
                if "bad_" in line:
                    idx_write = i + 1

                    # Comment existing mask for the 'replace' mode
                    if all([
                        write_mode == 'replace',
                        "bad_" in line.partition(';')[0]
                    ]):
                        contents[i] = "; " + line

                    # Check existing numbered 'bad_area's
                    sa_obj = re.search(r"bad_area(\d+)/", line)
                    if sa_obj is not None:
                        sa_num = int(sa_obj.group(1))
                        n_area_start = max(n_area_start, sa_num + 1)

                # Search for the 'rigid_group'
                # (to put mask before it in case no 'bad_' area found)
                if all([
                    idx_write == len(contents),
                    "rigid_group" in line
                ]):
                    idx_write = i

            text_before = contents[:idx_write]
            text_after = contents[idx_write:]

        # Format mask as a list of text lines
        text_mask = []
        if text_before and text_before[-1].strip() != "":
            text_mask.append("\n")

        geom_dict = load_crystfel_geometry(filename)

        new_bad_regions = self.make_crystfel_bad_regions(geom_dict['panels'])
        for i, bad_d in enumerate(new_bad_regions, start=n_area_start):
            text_mask.extend([
                f'bad_area{i}/{k} = {v}\n' for (k, v) in bad_d.items()
            ] + ['\n'])

        text_write = "".join(text_before + text_mask + text_after)

        with open(filename, 'w') as f_geo:
            f_geo.write(text_write)

    def __eq__(self, other):
        if not (isinstance(other, MaskRegions) and (self.shape == other.shape)):
            return False

        if set(self.regions) == set(other.regions):
            return True

        # Two equivalent masks could have regions described in different ways.
        # Converting them both to arrays is the easiest way to normalise them,
        # though probably not the most efficient.
        return np.array_equal(self.to_mask_array(), other.to_mask_array())

