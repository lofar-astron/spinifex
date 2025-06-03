#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from spinifex.vis_tools.ms_tools import (
    cli_get_dtec_h5parm_from_ms,
    cli_get_rm_h5parm_from_ms,
)


def get_parser():
    parser = argparse.ArgumentParser(
        description="This Spinifex command-line interface can make an H5Parm with ionospheric data (TEC or RM), using the metadata from a measurement set."
    )
    subparsers = parser.add_subparsers(dest="command")

    # Subparser for get_rm_h5parm_from_ms
    parser_rm = subparsers.add_parser(
        "get_rm_h5parm_from_ms",
        description="Calculate RM values using spinifex, write to h5parm",
    )
    parser_rm.add_argument(
        "ms",
        type=Path,
        help="Measurement set for which the RM values should be calculated.",
    )
    parser_rm.add_argument(
        "--iono-model-name",
        type=str,
        default="ionex",
        help="iono mode name",
        choices=["ionex", "ionex_iri", "tomion"],
    )
    parser_rm.add_argument(
        "--solset-name",
        type=str,
        help="Solset name. Default: create a new one based on first existing sol###",
    )
    parser_rm.add_argument(
        "--soltab-name", type=str, help="Soltab name. Default: rotationmeasure"
    )
    parser_rm.add_argument(
        "-o",
        "--h5parm",
        default="rotationmeasure.h5",
        type=Path,
        help="h5parm to which the results of the RotationMeasure is added.",
    )
    parser_rm.add_argument(
        "-a",
        "--add-to-existing-solset",
        action="store_true",
        help="Add to existing solset if it exists",
    )
    parser_rm.set_defaults(func=cli_get_rm_h5parm_from_ms)

    # Subparser for get_tec_h5parm_from_ms
    parser_tec = subparsers.add_parser(
        "get_tec_h5parm_from_ms",
        description="Calculate tec values using spinifex, write to h5parm",
    )
    parser_tec.add_argument(
        "ms",
        type=Path,
        help="Measurement set for which the tec values should be calculated.",
    )
    parser_tec.add_argument(
        "--iono-model-name",
        type=str,
        default="ionex",
        help="iono mode name",
        choices=["ionex", "ionex_iri", "tomion"],
    )
    parser_tec.add_argument(
        "--solset-name",
        type=str,
        help="Solset name. Default: create a new one based on first existing sol###",
    )
    parser_tec.add_argument(
        "--soltab-name", type=str, help="Soltab name. Default: tec###"
    )
    parser_tec.add_argument(
        "-o",
        "--h5parm",
        default="tec.h5",
        type=Path,
        help="h5parm to which the results of the tec is added.",
    )
    parser_tec.add_argument(
        "-a",
        "--add-to-existing-solset",
        action="store_true",
        help="Add to existing solset if it exists",
    )
    parser_tec.set_defaults(func=cli_get_dtec_h5parm_from_ms)

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
