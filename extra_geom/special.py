import numpy as np

from .detectors import JUNGFRAUGeometry as JF


def jungfrau_data(dc, name, use_mask=True, mask_edges=True, pix_edge=1,
                  trains=np.s_[:]):
    """Get Jungfrau data for the given detector in the data collection

    The returned data has inter ASIC gaps added and masked applied if requested

    dc: extra_data.DataCollection
        an extra_data data collection
    name: str
        name of the detector
    use_mask: bool
        if True, use the mask present in the data (if in use with CORR data)
    mask_edges: bool
        if True, mask the edges of the detector
    pix_edge: int
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
    jf = dc[name]
    jf = jf.select_trains(trains)
    data = jf['data.adc'].ndarray()

    if use_mask and 'data.mask' in jf.keys():
        mask = jf['data.mask'].ndarray().astype(np.bool_)
        data = data * ~mask

    if mask_edges:
        asic = np.zeros((JF.frag_ss_pixels, JF.frag_fs_pixels), dtype=np.bool_)
        asic[pix_edge:-pix_edge, pix_edge:-pix_edge] = True
        edge_mask = np.tile(
            asic, [JF.expected_data_shape[-2] // JF.frag_ss_pixels,
                   JF.expected_data_shape[-1] // JF.frag_fs_pixels])
        data = data * edge_mask

    assembled_data, centre = JF.from_module_positions().position_modules(data)
    return assembled_data, centre
