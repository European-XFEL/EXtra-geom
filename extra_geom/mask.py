"""Convert between
"""

import os
import re
import warnings

import h5py
import numpy as np


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


class MaskConverter:

    modes_avail = ('hd52geom', 'geom2hd5')
    write_modes = ('replace', 'add')

    def __init__(self, hd5file, geofile, run_mode, write_mode,
                 hd5path, hd5entry, detector, data_type, invert):
        """
        Construct a mask converter with provided parameters.

        Args:
            hd5file (string): relative or absolute path to the HD5 file.
            geofile (string): relative or absolute path to the geometry file.
            run_mode (string): mode of the mask converter operation.
            write_mode (string): mode of writing to the existing file.
            hd5path (string): path to the mask in HD5 file.
            hd5entry (int): entry number of the mask in HD5 file.
            detector (string): name of the detector.
            data_type (string): type of the detector data.
            invert (bool): invert the mask after reading from of before
                writing to the HD5 file.

        Raises:
            KeyError: unexpected run mode.
            KeyError: unexpected write mode.
            KeyError: detector name missing in the detector_info keys.
            KeyError: data type missing in detector_info for specified
                detector.
        """

        self._hd5file = hd5file
        self._hd5path = hd5path
        self._hd5entry = hd5entry
        self._geofile = geofile
        self._run_mode = run_mode
        self._write_mode = write_mode
        self._detector = detector
        self._data_type = data_type
        self._invert = invert

        if self._run_mode not in self.modes_avail:
            raise KeyError(
                f"Unknown MaskConverter run mode: {self._run_mode}.")
        if self._write_mode not in self.write_modes:
            raise KeyError(
                f"Unknown MaskConverter write mode: {self._write_mode}.")
        if self._detector not in di.detector_info.keys():
            raise KeyError(
                f"Unknown detector: {self._detector}.")
        if self._data_type not in di.detector_info[self._detector].keys():
            raise KeyError(
                f"Missing info on '{self._data_type}' data type for "
                f"detector {self._detector}.")

        self._det_info = di.detector_info[self._detector][self._data_type]

        self._read_mask()

    @property
    def mask(self):
        """
        Provides a copy of the detector mask.

        Returns:
            np.array: Detector mask as a 2D boolean numpy array.
        """
        return np.copy(self.__mask)

    def convert(self):
        """
        Convert detector mask and write to the output file.
        """
        self._convert_mask()
        self._write_mask()

    def _read_mask_hd5(self):
        """
        Read mask from the HD5 file to boolean numpy array.

        Returns:
            np.array: Detector mask as a 2D boolean numpy array.
        """

        res_mask = None
        mask_shape = self._det_info['shape']

        with h5py.File(self._hd5file, 'r') as f_hd5:
            hd5mask = f_hd5.get(self._hd5path)

            # Check shape of the HDF5 mask
            mu.check_hd5mask(hd5mask, mask_shape, self._hd5entry)

            # Remove 'n_data' dimension, convert to np.array
            if hd5mask.ndim > len(mask_shape):
                hd5mask = hd5mask[self._hd5entry]
            else:
                hd5mask = np.array(hd5mask)

        # Convert mask values to True (masked) and False
        res_mask = self._det_info['read_mask'](hd5mask)

        if self._invert:
            res_mask = np.logical_not(res_mask)

        return res_mask

    def _read_mask_geo(self):
        """
        Read mask from the geometry file to the dictionary of rectangles
        (in the same format as in the geometry file).

        Raises:
            ValueError: mask panel description not expected for the
                specified data type.
            ValueError: mask panel description does not suit specified
                detector and data type.

        Returns:
            dict: dictionary of rectangles (in the same format as in the
            geometry file) representing masked regions.
        """

        res_dict = {}
        bad_dict = {}

        with open(self._geofile, 'r') as f_geo:
            for line in f_geo:
                # To ignore comments:
                line = line.partition(';')[0].rstrip()

                m_obj = re.search(r"bad_(.+)/(\S+)\s*=\s*(\S+)", line)
                if m_obj is not None:
                    area, var, val = m_obj.groups()
                    if area not in bad_dict.keys():
                        bad_dict[area] = {
                            'min_fs': -1,
                            'max_fs': -1,
                            'min_ss': -1,
                            'max_ss': -1,
                            'panel': 'all'
                        }
                    if var == 'panel':
                        bad_dict[area][var] = val

                        # Check whether panel value is expected
                        if len(self._det_info['shape']) != 3:
                            raise ValueError(
                                f"Panel description ({val}) not expected "
                                f"for {self._data_type} data.")

                        # Check <val> to be suitable description
                        # of a panel and asic
                        m_obj = re.match(r"(p\d+)(a\d+)", val)
                        if (m_obj is None
                            or m_obj.group(1) not in
                                self._det_info['panel_names']
                            or m_obj.group(2) not in
                                self._det_info['asic_names']):
                            raise ValueError(
                                f"Not suitable panel description: {val}.")
                    elif var in bad_dict[area].keys():
                        bad_dict[area][var] = int(val)
                    else:
                        warnings.warn(
                            f"Geometry file - unsupported mask variable: "
                            f"{var} in bad_{area}.")

        # Check if rectangle information is complete
        i = 0
        for area in bad_dict.keys():
            if all([
                bad_dict[area]['min_fs'] >= 0,
                bad_dict[area]['max_fs'] >= 0,
                bad_dict[area]['min_ss'] >= 0,
                bad_dict[area]['max_ss'] >= 0
            ]):
                res_dict[i] = bad_dict[area]
                i += 1
            else:
                warnings.warn(
                    f"Geometry file - incomplete information for bad_{area}.")

        return res_dict

    def _read_mask(self):
        """
        Read mask from the HD5 or geometry file, depending on mode.
        """

        self.__mask = np.zeros(self._det_info['shape'], dtype=bool)
        self.__rect = {}

        if self._run_mode == "hd52geom":
            self.__mask = self._read_mask_hd5()

            # Reduce mask in case of write option 'add'
            if self._write_mode == 'add' and os.path.exists(self._geofile):
                self.__rect = self._read_mask_geo()
                reduce_mask = self._convert_rectd2nparr()
                self.__mask = np.logical_and(self.__mask,
                                             np.logical_not(reduce_mask))
                self.__rect = {}

        elif self._run_mode == "geom2hd5":
            self.__rect = self._read_mask_geo()

            # Read also mask from HD5 in case of write option 'add'
            if (self._write_mode == 'add'
                    and os.path.exists(self._hd5file)):
                self.__mask = self._read_mask_hd5()

    def _convert_nparr2rectd(self):
        """
        Convert mask from the 2D boolean numpy array to the dictionary
        of rectangles (in the same format as in the geometry file).

        Returns:
            dict: Dictionary of rectangles (in the same format as in the
            geometry file) representing masked regions.
        """

        res_dict = {}

        # First check for regions to be excluded in all panels:
        if self.__mask.ndim == 3:
            panel_all = self.mask[0]
            for i in range(1, self.__mask.shape[0]):
                panel_all = np.logical_and(panel_all, self.__mask[i])
        else:
            panel_all = self.mask

        is_panel_all_empty = np.array_equiv(panel_all, False)
        if not is_panel_all_empty:
            res_dict.update(mu.rect2dict(delta_method(panel_all), 'all'))

        if self.__mask.ndim == 3:
            panels = self._det_info['panel_names']
            asics = self._det_info['asic_names']
            asic_range = self._det_info['asic_range']

            # Loop over all panels:
            for i in range(len(panels)):
                if is_panel_all_empty:
                    panel_i = self.__mask[i]
                else:
                    panel_i = np.copy(self.__mask[i])
                    panel_i = np.logical_and(panel_i,
                                             np.logical_not(panel_all))

                # Loop over all asics in the panel:
                for j in range(len(asics)):
                    asic_j = np.zeros(panel_i.shape, dtype=bool)
                    asic_j[asic_range[i][j]] = True
                    panel_i_asic_j = np.logical_and(panel_i, asic_j)
                    res_dict.update(
                        mu.rect2dict(delta_method(panel_i_asic_j),
                                     f"{panels[i]}{asics[j]}"))

        return res_dict

    def _convert_rectd2nparr(self):
        """
        Convert mask from the dictionary of rectangles (same format as
        in the geometry file) to the 2D boolean numpy array.

        Returns:
            np.array: Detector mask as a 2D boolean numpy array.
        """

        shape = self._det_info['shape']
        res_mask = np.zeros(shape, dtype=bool)

        for area in self.__rect.keys():
            slice_ss = slice(self.__rect[area]['min_ss'],
                             self.__rect[area]['max_ss'] + 1)
            slice_fs = slice(self.__rect[area]['min_fs'],
                             self.__rect[area]['max_fs'] + 1)
            if self.__rect[area]['panel'] == 'all':
                if len(shape) == 3:
                    res_mask[:, slice_ss, slice_fs] = True
                else:
                    res_mask[slice_ss, slice_fs] = True
            else:
                assert len(shape) == 3, (
                    "Convert rectd2nparr: mask has to be dimensions 3 to "
                    "apply rectangles per panel.")

                panel, asic = re.match(
                    r"(p\d+)(a\d+)", self.__rect[area]['panel']).groups()
                panel_n = self._det_info['panel_names'].index(panel)
                asic_n = self._det_info['asic_names'].index(asic)
                asic_range = self._det_info['asic_range'][panel_n][asic_n]

                slice_ss_in_asic = slice(
                    max(slice_ss.start, asic_range[0].start),
                    min(slice_ss.stop, asic_range[0].stop))
                slice_fs_in_asic = slice(
                    max(slice_fs.start, asic_range[1].start),
                    min(slice_fs.stop, asic_range[1].stop))

                res_mask[panel_n, slice_ss_in_asic, slice_fs_in_asic] = True

        return res_mask

    def _convert_mask(self):
        """
        Convert the mask, depending on mode.
        """
        if self._run_mode == "hd52geom":
            self.__rect = self._convert_nparr2rectd()
        elif self._run_mode == "geom2hd5":
            rect_mask = self._convert_rectd2nparr()
            self.__mask = np.logical_or(self.__mask, rect_mask)

    def _write_mask_hd5(self):
        """
        Write converted mask to the HD5 file.
        """

        mask_shape = self.__mask.shape
        mask_tmp = self.__mask
        if self._invert:
            mask_tmp = np.logical_not(mask_tmp)
        mask_to_write = self._det_info['write_mask'](mask_tmp)

        with h5py.File(self._hd5file, 'a') as f_hd5:
            if self._hd5path in f_hd5:
                hd5mask = f_hd5[self._hd5path]

                # Check shape of the existing HDF5 mask
                mu.check_hd5mask(hd5mask, mask_shape, self._hd5entry)

                if hd5mask.ndim > len(mask_shape):
                    hd5mask[self._hd5entry] = mask_to_write
                else:
                    hd5mask[...] = mask_to_write
            else:
                f_hd5.create_dataset(self._hd5path, data=mask_to_write)

    def _write_mask_geo(self):
        """
        Write converted mask to the geometry file.
        """

        text_before = text_after = []
        n_area_start = 0

        # Store and process content of the existing geometry file
        if os.path.exists(self._geofile):
            with open(self._geofile, 'r') as f_geo:
                contents = f_geo.readlines()

            idx_write = len(contents)
            for i, line in enumerate(contents):

                # Search for the mask information
                if "bad_" in line:
                    idx_write = i + 1

                    # Comment existing mask for the 'replace' mode
                    if all([
                        self._write_mode == 'replace',
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
        if (text_before
                and text_before[-1].strip() != ""):
            text_mask.append("\n")

        for area in self.__rect:
            for dim in ['min_fs', 'max_fs', 'min_ss', 'max_ss']:
                text_mask.append(
                    f"bad_area{n_area_start + area}/{dim} = "
                    f"{self.__rect[area][dim]}\n")
            if self.__rect[area]['panel'] != 'all':
                text_mask.append(
                    f"bad_area{n_area_start + area}/panel = "
                    f"{self.__rect[area]['panel']}\n")
            text_mask.append("\n")

        text_write = "".join(text_before + text_mask + text_after)

        with open(self._geofile, 'w') as f_geo:
            f_geo.write(text_write)

    def _write_mask(self):
        """
        Write converted mask to the HD5 or geometry file, depending on mode.
        """
        # Write the mask depending on mode.
        if self._run_mode == "hd52geom":
            self._write_mask_geo()
        elif self._run_mode == "geom2hd5":
            self._write_mask_hd5()

