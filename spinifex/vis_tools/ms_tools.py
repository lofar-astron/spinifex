from __future__ import annotations

from pathlib import Path
from typing import Any, NamedTuple

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from numpy.typing import NDArray

from spinifex.get_rm import DEFAULT_IONO_HEIGHT, RM, get_rm_from_skycoord
from spinifex.ionospheric import ModelDensityFunction, ionospheric_models
from spinifex.magnetic import MagneticFieldFunction, magnetic_models

try:
    import casacore.tables as tab
except ImportError as e:
    MSG = "casacore is not installed! To operate on MeasurementSets, install spinifex[casa]."
    raise ImportError(MSG) from e


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
        tab.taql("select gmin(TIME_CENTROID), gmax(TIME_CENTROID) from $ms_path")[
            0
        ].values()
    )
    with tab.table(ms_path.as_posix()) as my_ms:
        timestep = my_ms.getcell("INTERVAL", 0)
        times = Time(
            np.arange(timerange[0], timerange[1] + 0.5 * timestep, timestep)
            / (24 * 3600),
            format="mjd",
        )
        pointing = tab.table(my_ms.getkeyword("FIELD")).getcell("PHASE_DIR", 0)[0]
        stations = tab.table(my_ms.getkeyword("ANTENNA")).getcol("NAME")
        station_pos = tab.table(my_ms.getkeyword("ANTENNA")).getcol("POSITION")
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
    with tab.table(ms_path.as_posix()) as my_ms:
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
    return location[0]


def get_rm_from_ms(
    ms_path: Path,
    timestep: u.Quantity | None = None,
    use_stations: list[int | str] | None = None,
    height_array: NDArray[np.float64] = DEFAULT_IONO_HEIGHT,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
    magnetic_model: MagneticFieldFunction = magnetic_models.ppigrf,
    iono_kwargs: dict[str, Any] | None = None,
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
    height_array : NDArray[np.float64], optional
        array of ionospheric altitudes, by default DEFAULT_IONO_HEIGHT
    iono_model : ModelDensityFunction, optional
        ionospheric model, by default ionospheric_models.ionex
    magnetic_model : MagneticFieldFunction, optional
        geomagnetic model, by default magnetic_models.ppigrf
    iono_kwargs : dict[str, Any] | None, optional
        arguments for the ionospheric model, by default None

    Returns
    -------
    dict[str, RM]
        dictionary with RM object per station
    """
    iono_kwargs = iono_kwargs or {}
    rm_dict = {}
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
    # TODO: implement use_stations is all (default?)
    if use_stations is None:
        location = get_average_location(ms_metadata.locations)
        rm_dict["average_station_pos"] = get_rm_from_skycoord(
            loc=location,
            times=times,
            source=ms_metadata.source,
            height_array=height_array,
            iono_model=iono_model,
            magnetic_model=magnetic_model,
            iono_kwargs=iono_kwargs,
        )
    else:
        # get rm per station
        for stat in use_stations:
            if isinstance(stat, str):
                istat = ms_metadata.station_names.index(stat)
            else:
                istat = stat
            rm_dict[ms_metadata.station_names[istat]] = get_rm_from_skycoord(
                loc=ms_metadata.locations[istat],
                times=times,
                source=ms_metadata.source,
                height_array=height_array,
                iono_model=iono_model,
                magnetic_model=magnetic_model,
                iono_kwargs=iono_kwargs,
            )
    return rm_dict
