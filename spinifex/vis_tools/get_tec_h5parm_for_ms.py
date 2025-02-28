#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from spinifex import h5parm_tools
from spinifex.vis_tools import ms_tools

logger = logging.getLogger(__name__)


def get_tecs_for_ms(
    msname,
    ionex_dir="./ionex_data/",
    ionex_prefix="uqr",
    ionex_server="chapman",
):
    if not Path(msname).is_dir():
        message = f"Measurement set {msname} is not valid."
        raise RuntimeError(message)

    iono_kwargs = {
        "output_directory": Path(ionex_dir),
        "prefix": ionex_prefix,
        "server": ionex_server,
    }

    ms_metadata = ms_tools.get_metadata_from_ms(Path(msname))
    return ms_tools.get_dtec_from_ms(
        Path(msname), iono_kwargs=iono_kwargs, use_stations=ms_metadata.station_names
    )


def main():
    # Initialize parser
    parser = argparse.ArgumentParser(
        description="Calculate tec values using spinifex, add to h5parm"
    )
    parser.add_argument(
        "ms",
        type=str,
        help="Measurement set for which the tec values should be calculated.",
    )
    parser.add_argument(
        "--solset-name",
        type=str,
        help="Solset name. Default: create a new one based on first existing sol###",
    )
    parser.add_argument("--soltab-name", type=str, help="Soltab name. Default: tec###")
    parser.add_argument(
        "-o",
        "--h5parm",
        default="tec.h5",
        type=str,
        help="h5parm to which the results of the tec is added.",
    )
    parser.add_argument(
        "-a",
        "--add-to-existing-solset",
        action="store_true",
        help="Add to existing solset if it exists",
    )

    args = parser.parse_args()

    dtec = get_tecs_for_ms(args.ms)
    h5parm_tools.write_tec_to_h5parm(
        dtec,
        args.h5parm,
        solset_name=args.solset_name,
        soltab_name=args.soltab_name,
        add_to_existing_solset=args.add_to_existing_solset,
    )


if __name__ == "__main__":
    main()
