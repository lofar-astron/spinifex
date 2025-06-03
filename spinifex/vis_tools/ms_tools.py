from __future__ import annotations

import argparse
from concurrent.futures import Future, ThreadPoolExecutor
from functools import partial
from pathlib import Path
from typing import Any, Literal, NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex import h5parm_tools
from spinifex.get_dtec import DTEC, _get_dtec_from_skycoord
from spinifex.get_rm import (
    DEFAULT_IONO_HEIGHT,
    RM,
    _get_rm_from_skycoord,
)
from spinifex.ionospheric import ModelDensityFunction
from spinifex.ionospheric.iri_density import IRI_HEIGHTS
from spinifex.ionospheric.models import O, parse_iono_kwargs, parse_iono_model
from spinifex.ionospheric.tomion_parser import TOMION_HEIGHTS
from spinifex.logger import logger
from spinifex.magnetic import MagneticFieldFunction
from spinifex.magnetic.models import magnetic_models, parse_magnetic_model

try:
    from casacore.tables import table as _casacore_table
    from casacore.tables import taql
except ImportError as e:
    MSG = "casacore is not installed! To operate on MeasurementSets, install spinifex[casacore]."
    raise ImportError(MSG) from e

# Disable acknowledgement from opening casacore tables
table = partial(_casacore_table, ack=False)


class MsMetaData(NamedTuple):
    """Metadata from a Measurement Set"""

    times: Time
    locations: EarthLocation
    station_names: list[str]
    name: str
    source: SkyCoord


def get_metadata_from_ms(ms_path: Path) -> MsMetaData:
    """open measurement set and get metadata from it

    Parameters
    ----------
    ms_path : Path
        measurement set

    Returns
    -------
    MsMetaData
        object with metadata
    """
    timerange = list(
        taql("select gmin(TIME_CENTROID), gmax(TIME_CENTROID) from $ms_path")[
            0
        ].values()
    )
    with table(ms_path.as_posix()) as my_ms:
        timestep = my_ms.getcell("INTERVAL", 0)
        times = Time(
            np.arange(timerange[0], timerange[1] + 0.5 * timestep, timestep)
            / (24 * 3600),
            format="mjd",
        )
        pointing = table(my_ms.getkeyword("FIELD")).getcell("PHASE_DIR", 0)[0]
        stations = table(my_ms.getkeyword("ANTENNA")).getcol("NAME")
        station_pos = table(my_ms.getkeyword("ANTENNA")).getcol("POSITION")
        locations = EarthLocation.from_geocentric(*station_pos.T, unit=u.m)
        return MsMetaData(
            times=times,
            locations=locations,
            station_names=stations,
            name=ms_path.as_posix(),
            source=SkyCoord(pointing[0] * u.rad, pointing[1] * u.rad),
        )


def get_columns_from_ms(ms_path: Path) -> list[str]:
    """Get the columns from a MeasurementSet"""
    with table(ms_path.as_posix()) as my_ms:
        return list(my_ms.colnames())


def get_average_location(location: EarthLocation) -> EarthLocation:
    # TODO; implement correctly in NE plane
    """Get first location from N locations


    Parameters
    ----------
    location : EarthLocation
        N locations

    Returns
    -------
    EarthLocation
        first location
    """
    logger.warning("Using first location from list of locations - not a true average!")
    return location[0]


def get_rm_from_ms(
    ms_path: Path,
    timestep: u.Quantity | None = None,
    use_stations: list[int] | list[str] | Literal["all", "average"] = "average",
    iono_model_name: str | None = None,
    magnetic_model_name: str = "ppigrf",
    **iono_kwargs: Any,
) -> dict[str, RM]:
    """Get rotation measures for a measurement set

    Parameters
    ----------
    ms_path : Path
        measurement set
    timestep : u.Quantity | None, optional
        only calculate rotation measure every timestep, by default None
    use_stations : list[int  |  str] | None, optional
        list of stations (index or name) to use,
        if None use first of the measurement set, by default None
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    magnetic_model_name : str, optional
        geomagnetic model name, by default "ppigrf". Must be a supported geomagnetic model.
    iono_kwargs : dict, optional
        arguments for the ionospheric model, by default None

    Returns
    -------
    dict[str, RM]
        dictionary with RM object per station
    """
    if iono_model_name is None:
        iono_model_name = "ionex"

    iono_model = parse_iono_model(iono_model_name)
    height_array = DEFAULT_IONO_HEIGHT
    if iono_model_name == "tomion":
        height_array = TOMION_HEIGHTS
    elif iono_model_name == "ionex_iri":
        height_array = IRI_HEIGHTS
    iono_options = parse_iono_kwargs(iono_model=iono_model, **iono_kwargs)
    magnetic_model = parse_magnetic_model(magnetic_model_name)
    return _get_iono_from_ms(
        ms_path=ms_path,
        iono_type="rm",
        timestep=timestep,
        use_stations=use_stations,
        height_array=height_array,
        iono_model=iono_model,
        magnetic_model=magnetic_model,
        iono_options=iono_options,
    )


def get_dtec_from_ms(
    ms_path: Path,
    timestep: u.Quantity | None = None,
    use_stations: list[int] | list[str] | Literal["all", "average"] = "average",
    iono_model_name: str | None = None,
    **iono_kwargs: Any,
) -> dict[str, NDArray]:
    """Get rotation measures for a measurement set

    Parameters
    ----------
    ms_path : Path
        measurement set
    timestep : u.Quantity | None, optional
        only calculate rotation measure every timestep, by default None
    use_stations : list[int  |  str] | None, optional
        list of stations (index or name) to use,
        if None use first of the measurement set, by default None
    iono_model_name : str, optional
        ionospheric model name, by default "ionex". Must be a supported ionospheric model.
    iono_kwargs : dict, optional
        arguments for the ionospheric model, by default None

    Returns
    -------
    dict[str, NDArray]
        dictionary with electron_density_profiles per station
    """

    if iono_model_name is None:
        iono_model_name = "ionex"

    iono_model = parse_iono_model(iono_model_name)
    height_array = DEFAULT_IONO_HEIGHT
    if iono_model_name == "tomion":
        height_array = TOMION_HEIGHTS
    elif iono_model_name == "ionex_iri":
        height_array = IRI_HEIGHTS

    iono_options = parse_iono_kwargs(iono_model, **iono_kwargs)
    return _get_iono_from_ms(
        ms_path=ms_path,
        iono_type="dtec",
        timestep=timestep,
        use_stations=use_stations,
        height_array=height_array,
        iono_model=iono_model,
        iono_options=iono_options,
    )


def _get_iono_from_ms(
    ms_path: Path,
    iono_model: ModelDensityFunction[O],
    iono_type: Literal["dtec", "rm"],
    timestep: u.Quantity | None = None,
    use_stations: list[int] | list[str] | Literal["all", "average"] = "average",
    height_array: NDArray[np.float64] = DEFAULT_IONO_HEIGHT,
    iono_options: O | None = None,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
) -> dict[str, RM] | dict[str, NDArray]:
    """Get ionospheric values for a measurement set

    Parameters
    ----------
    ms_path : Path
        measurement set
    iono_type : Literal["dtec", "rm"]
        type of ionospheric value to calculate
    timestep : u.Quantity | None, optional
        only calculate rotation measure every timestep, by default None
    use_stations : list[int  |  str] | None, optional
        list of stations (index or name) to use,
        if None use first of the measurement set, by default None
    height_array : NDArray[np.float64], optional
        array of ionospheric altitudes, by default DEFAULT_IONO_HEIGHT
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    magnetic_model : MagneticFieldFunction, optional
        geomagnetic model, by default magnetic_models.ppigrf
    iono_options : OptionType | None, optional
        arguments for the ionospheric model, by default None

    Returns
    -------
    dict[str, RM] | dict[str, NDArray]
        dictionary with ionospheric values
    """

    if iono_type not in ["dtec", "rm"]:
        msg = f"iono_type should be 'dtec' or 'rm'. Got {iono_type}"
        raise ValueError(msg)

    get_func = (
        _get_dtec_from_skycoord
        if iono_type == "dtec"
        # Set the magnetic model here for get_rm
        # all other arguments are the same
        else partial(_get_rm_from_skycoord, magnetic_model=magnetic_model)
    )

    result_dict: dict[str, RM | DTEC] = {}

    ms_metadata = get_metadata_from_ms(ms_path)
    if timestep is not None:
        # dtime = ms_metadata.times[1].mjd - ms_metadata.times[0].mjd
        dtime_in_days = timestep.to(u.hr).value / 24
        times = Time(
            np.arange(
                ms_metadata.times[0].mjd - 0.5 * dtime_in_days,
                ms_metadata.times[-1].mjd + 0.5 * dtime_in_days,
                dtime_in_days,
            ),
            format="mjd",
        )
    else:
        times = ms_metadata.times

    if isinstance(use_stations, str):
        if use_stations == "average":
            location = get_average_location(ms_metadata.locations)
            result_dict["average_station_pos"] = get_func(  # type: ignore[operator]
                loc=location,
                times=times,
                source=ms_metadata.source,
                height_array=height_array,
                iono_model=iono_model,
                iono_options=iono_options,
            )
            return result_dict
        if use_stations == "all":
            station_list: list[str] | list[int] = ms_metadata.station_names
        else:
            msg = f"`use_stations` should be a list of stations, or string literal 'average' or 'all'. Got {use_stations}"
            raise ValueError(msg)

    elif isinstance(use_stations, list):
        station_list = use_stations
    else:
        msg = f"`use_stations` should be a list of stations, or string literal 'average' or 'all'. Got {use_stations}"
        raise ValueError(msg)
    # get rm per station
    logger.info("Getting DTEC per station")

    # get rm per station
    # Submit jobs to ThreadPoolExecutor
    # Gets Future objects for each station
    future_dict: dict[str, Future[RM | DTEC]] = {}
    key_names: list[str] = []

    # Execute the first station in serial
    # Avoids race conidition to download ionex files

    stat = station_list[0]
    istat = ms_metadata.station_names.index(stat) if isinstance(stat, str) else stat
    key_name = ms_metadata.station_names[istat]
    result_dict[key_name] = get_func(  # type: ignore[operator]
        loc=ms_metadata.locations[istat],
        times=times,
        source=ms_metadata.source,
        height_array=height_array,
        iono_model=iono_model,
        iono_options=iono_options,
    )

    with ThreadPoolExecutor() as executor:
        # Skip the first station
        for stat in station_list[1:]:
            istat = (
                ms_metadata.station_names.index(stat) if isinstance(stat, str) else stat
            )
            key_name = ms_metadata.station_names[istat]
            key_names.append(key_name)
            future_dict[key_name] = executor.submit(
                get_func,  # type: ignore[arg-type]
                loc=ms_metadata.locations[istat],
                times=times,
                source=ms_metadata.source,
                height_array=height_array,
                iono_model=iono_model,
                iono_options=iono_options,
            )

    # Get concrete results from futures
    # Resolve Future objects to get the results
    for key_name in key_names:
        result_dict[key_name] = future_dict[key_name].result()
    return result_dict


def cli_get_rm_h5parm_from_ms(args: argparse.Namespace) -> None:
    ms_path: Path = Path(args.ms)
    h5parm_path: str = args.h5parm
    solset_name: str | None = args.solset_name
    soltab_name: str | None = args.soltab_name
    add_to_existing_solset: bool = args.add_to_existing_solset
    iono_model_name: str | None = args.iono_model_name

    rm_dict = get_rm_from_ms(
        ms_path=ms_path, use_stations="all", iono_model_name=iono_model_name
    )
    h5parm_tools.write_rm_to_h5parm(
        rm_dict,
        h5parm_path,
        solset_name=solset_name,
        soltab_name=soltab_name,
        add_to_existing_solset=add_to_existing_solset,
    )


def cli_get_dtec_h5parm_from_ms(args: argparse.Namespace) -> None:
    ms_path: Path = Path(args.ms)
    h5parm_path: str = args.h5parm
    solset_name: str | None = args.solset_name
    soltab_name: str | None = args.soltab_name
    add_to_existing_solset: bool = args.add_to_existing_solset
    iono_model_name: str | None = args.iono_model_name

    dtec = get_dtec_from_ms(
        ms_path=ms_path, use_stations="all", iono_model_name=iono_model_name
    )
    h5parm_tools.write_tec_to_h5parm(
        dtec,
        h5parm_path,
        solset_name=solset_name,
        soltab_name=soltab_name,
        add_to_existing_solset=add_to_existing_solset,
    )
