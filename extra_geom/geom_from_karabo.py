# coding: utf-8
import argparse
import asyncio
import base64
import gzip
import pickle

from karabo_proxy import AsyncKaraboProxy

FRM_DIMS = {
    "crystfel": ('frame', 'modno', 'ss', 'fs'),
    "cheetah": ('frame', 'ss', 'fs'),
}

presets = {
    "SPB/AGIPD1M": {
        "proxy_url": "http://spb-rr-sys-con-gui1:8484",
        "geometry_device": "SPB_DET_AGIPD1M-1/CAL/GEOMETRY",
        "z_motor": "SPB_IRU_AGIPD1M/MOTOR/Z_STEPPER",
        "energy_device": "SPB_XTD2_UND/DOOCS/ENERGY",
        "clen_min": 0.119,
    },
    "SPB/JF4M": {
        "proxy_url": "http://spb-rr-sys-con-gui1:8484",
        "geometry_device": "SPB_IRDA_JF4M/CAL/GEOMETRY",
        "z_motor": "SPB_IRDA_JF4M/MOTOR/Z",
        "energy_device": "SPB_XTD2_UND/DOOCS/ENERGY",
        "clen_min": 0.020,
    },
    "MID": {
        "proxy_url": "http://mid-rr-sys-con-gui1:8484",
        "geometry_device": "MID_DET_AGIPD1M-1/CAL/GEOMETRY",
        "z_motor": "MID_EXP_UPP/MOTOR/T4",
        "energy_device": "MID_XTD1_UND/DOOCS/ENERGY",
        "clen_min": 0.0,
    },
}


def deserialize_geometry(serialized: str):
    return pickle.loads(gzip.decompress(base64.b64decode(serialized.encode())))


async def get_geometry(proxy_url: str, geometry_device: str,
                       z_motor: str, energy_device: str, clen_min: float):
    client = AsyncKaraboProxy(proxy_url)
    data = await asyncio.gather(
        client.get_device_config_path(geometry_device, "serializedGeometry"),
        client.get_device_config_path(z_motor, "actualPosition"),
        client.get_device_config_path(energy_device, "actualPosition"),
    )
    geom = deserialize_geometry(data[0].value)
    meta = type("GeometryMetadata", (), dict(
        clen=data[1].value * 0.001 + clen_min,
        photon_energy=int(1000 * data[2].value + 0.5)
    ))
    return geom, meta


def main():
    parser = argparse.ArgumentParser(
        description="Dump a crystfel geometry from Karabo to a file"
    )
    parser.add_argument('detector', choices=list(presets.keys()),
                        help="Detector: " + ", ".join(presets.keys()))
    parser.add_argument('output_file',
                        help="The name of output file")
    parser.add_argument('-f', '--format',
                        choices=["crystfel", "cheetah"], default="crystfel",
                        help="The output format dialect, default: crystfel")
    parser.add_argument('--adu_per_eV', type=float, default=0.001,
                        help="Unit factor ADU/eV, default=0.001")
    parser.add_argument('--tube_len', type=float, default=0.0,
                        help="Lenght of additional spacing between detector "
                        "and sample, e.g. tubes length")

    args = parser.parse_args()
    prm = presets[args.detector].copy()
    prm["clen_min"] += args.tube_len

    geom, meta = asyncio.run(get_geometry(**prm))

    geom.write_crystfel_geom(
        args.output_file,
        dims=FRM_DIMS[args.format],
        adu_per_ev=args.adu_per_eV,
        clen=meta.clen,
        photon_energy=meta.photon_energy
    )


if __name__ == "__main__":
    main()
