"""Write geometry in CrystFEL format.
"""
from itertools import product

import numpy as np

HEADER_TEMPLATE = """\
; {detector} geometry file written by EXtra-geom {version}
; You may need to edit this file to add:
; - data and mask locations in the file
; - mask_good & mask_bad values to interpret the mask
; - adu_per_eV & photon_energy
; - clen (detector distance)
;
; See: http://www.desy.de/~twhite/crystfel/manual-crystfel_geometry.html

{paths}
{frame_dim}
res = {resolution} ; pixels per metre

; Beam energy in eV
{photon_energy}

; Camera length, aka detector distance
{clen}

; Analogue Digital Units per eV
{adu_per_ev}
"""


PANEL_TEMPLATE = """
{dims}
{name}/min_fs = {min_fs}
{name}/min_ss = {min_ss}
{name}/max_fs = {max_fs}
{name}/max_ss = {max_ss}
{name}/fs = {fs_vec}
{name}/ss = {ss_vec}
{name}/corner_x = {corner_x}
{name}/corner_y = {corner_y}
{name}/coffset = {coffset}
"""


def _crystfel_format_vec(vec):
    """Convert an array of 3 numbers to CrystFEL format like "+1.0x -0.1y"
    """
    s = '{:+}x {:+}y'.format(*vec[:2])
    if vec[2] != 0:
        s += ' {:+}z'.format(vec[2])
    return s


def frag_to_crystfel(fragment, p, a, ss_slice, fs_slice, dims, pixel_size):
    tile_name = 'p{}a{}'.format(p, a)
    c = fragment.corner_pos[:2] / pixel_size
    dim_list = []
    for num, value in dims.items():
        if value == 'modno':
            key = p
        else:
            key = value
        dim_list.append('{}/dim{} = {}'.format(tile_name, num, key))

    return PANEL_TEMPLATE.format(
        dims='\n'.join(dim_list),
        name=tile_name,
        min_ss=ss_slice.start,
        max_ss=ss_slice.stop - 1,
        min_fs=fs_slice.start,
        max_fs=fs_slice.stop - 1,
        ss_vec=_crystfel_format_vec(fragment.ss_vec / pixel_size),
        fs_vec=_crystfel_format_vec(fragment.fs_vec/ pixel_size),
        corner_x=c[0],
        corner_y=c[1],
        coffset=fragment.corner_pos[2],
    )

def write_crystfel_geom(self, filename, *,
                        data_path='/entry_1/instrument_1/detector_1/data',
                        mask_path=None, dims=('frame', 'modno', 'ss', 'fs'),
                        nquads=4, adu_per_ev=None, clen=None,
                        photon_energy=None):
    """Write this geometry to a CrystFEL format (.geom) geometry file.
    """
    from . import __version__

    if adu_per_ev is None:
        adu_per_ev_str = '; adu_per_eV = SET ME'
        # TODO: adu_per_ev should be fixed for each detector, we should
        #       find out the values and set them.
    else:
        adu_per_ev_str = 'adu_per_eV = {}'.format(adu_per_ev)

    if clen is None:
        clen_str = '; clen = SET ME'
    else:
        clen_str = 'clen = {}'.format(clen)

    if photon_energy is None:
        photon_energy_str = '; photon_energy = SET ME'
    else:
        photon_energy_str = 'photon_energy = {}'.format(photon_energy)

    # Get the frame dimension
    tile_dims = {}

    frame_dim = None
    for nn, dim_name in enumerate(dims):
        if dim_name == 'frame':
            frame_dim = 'dim{} = %'.format(nn)
        else:
            tile_dims[nn] = dim_name
    if frame_dim is None:
        raise ValueError('No frame dimension given')

    panel_chunks = []
    for p, module in enumerate(self.modules):
        for a, fragment in enumerate(module):
            ss_slice, fs_slice = self._tile_slice(a)
            if 'modno' not in dims:
                # If we don't have a modno dimension, assume modules are
                # concatenated along the slow-scan dim, e.g. AGIPD (8192, 128)
                module_offset = p * self.expected_data_shape[1]
                ss_slice = slice(
                    ss_slice.start + module_offset,
                    ss_slice.stop + module_offset
                )

            panel_chunks.append(frag_to_crystfel(
                fragment, p, a, ss_slice, fs_slice, tile_dims, self.pixel_size
            ))

    resolution = 1.0 / self.pixel_size  # Pixels per metre

    paths = dict(data=data_path)
    if mask_path:
        paths['mask'] = mask_path
    path_str = '\n'.join('{} = {} ;'.format(i, j) for i, j in paths.items())

    with open(filename, 'w') as f:
        f.write(HEADER_TEMPLATE.format(
            detector=self.detector_type_name,
            version=__version__,
            paths=path_str,
            frame_dim=frame_dim,
            resolution=resolution,
            adu_per_ev=adu_per_ev_str,
            clen=clen_str,
            photon_energy=photon_energy_str
        ))
        rigid_groups = get_rigid_groups(self, nquads=nquads)
        f.write(rigid_groups)
        for chunk in panel_chunks:
            f.write(chunk)

def get_rigid_groups(geom, nquads=4):
    """Create string for rigid groups definition."""

    lines = []

    module_groups = {
        f'p{p}': [f"p{p}a{a}" for a in range(len(tiles))]
        for p, tiles in enumerate(geom.modules)
    }
    lines += [
        f'rigid_group_{k} = ' + ','.join(v) for k, v in module_groups.items()
    ] + [
        'rigid_group_collection_modules = ' + ','.join(module_groups),
        ''
    ]

    if nquads:
        if len(geom.modules) % nquads:
            raise ValueError(
                f"Can't divide {len(geom.modules)} modules into {nquads} groups"
            )
        prod = product(range(len(geom.modules)), range(geom.n_tiles_per_module))
        all_panels = [f"p{p}a{a}" for (p, a) in prod]
        quad_groups = {
            f'q{q}': panels
            for q, panels in enumerate(np.array_split(all_panels, nquads))
        }
        lines += [
            f'rigid_group_{k} = ' + ','.join(v) for k, v in quad_groups.items()
        ] + [
            'rigid_group_collection_quads = ' + ','.join(quad_groups),
            ''
        ]

    return '\n'.join(lines)
