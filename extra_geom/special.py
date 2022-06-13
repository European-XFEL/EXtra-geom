from warnings import warn

import numpy as np

from .detectors import JUNGFRAUGeometry as JF


def jungfrau_module_data(data, use_mask=True, edge_pixels=1, trains=np.s_[:]):
    """Get data for the given JUNGFRAU detector

    The returned data has inter ASIC gaps added and masked applied if requested

    ::

        from extra_data import open_run

        run = open_run(1234, 1)
        jungfrau = run['HED_IA1_JF500K1/DET/JNGFR01:daqOutput']
        data = jungfrau_module_data(jungfrau)

    data: extra_data.SourceData
        A SourceData object for a JUNGFRAU detector
    use_mask: bool
        if True, use the mask present in the data (if in use with CORR data)
    edge_pixels: int
        number of pixels to mask at the edges of detector ASICS
    trains: slice
        Select a subset of trains from this data collection.
            Slice trains by position within this data::
                np.s_[:5]  # first 5 trains
            Or select trains by train ID, with a slice or a list::
                from extra_data import by_id
                by_id[142844490 : 142844495]
                by_id[[142844490, 142844493, 142844494]]
    """
    jf = data.select_trains(trains)
    data = jf['data.adc'].ndarray()

    if use_mask:
        if 'data.mask' in jf.keys():
            mask = jf['data.mask'].ndarray().astype(np.bool_)
            data = data * ~mask
        else:
            warn("No mask present in data")

    if edge_pixels > 0:
        asic = np.zeros((JF.frag_ss_pixels, JF.frag_fs_pixels), dtype=np.bool_)
        asic[edge_pixels:-edge_pixels, edge_pixels:-edge_pixels] = True
        edge_mask = np.tile(
            asic, [JF.expected_data_shape[-2] // JF.frag_ss_pixels,
                   JF.expected_data_shape[-1] // JF.frag_fs_pixels])
        data = data * edge_mask

    assembled_data, centre = JF.from_module_positions().position_modules(data)
    return assembled_data, centre
