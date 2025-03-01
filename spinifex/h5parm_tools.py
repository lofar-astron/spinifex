from __future__ import annotations

from pathlib import Path
from typing import Any

import astropy.units as u
import h5py
import numpy as np
from numpy.typing import NDArray

from spinifex.get_dtec import DTEC
from spinifex.get_rm import RM


def get_minimal_solset_name(h5file: h5py.File) -> str:
    """Find the best name for your new solset in an h5parm file. Default name is sol###

    Parameters
    ----------
    h5file : h5py.File
        the h5parm file

    Returns
    -------
    str
        frst available name for the solset
    """
    for solset_idx in range(1000):
        if f"sol{solset_idx:03d}" not in h5file:
            return f"sol{solset_idx:03d}"
    return "sol9999"


def get_minimal_soltab_name(solset: h5py.Dataset, soltab_type: str) -> str:
    """Find the best name for your new soltab in an h5parm file. Default name is [soltab_type]###

    Parameters
    ----------
    solset : h5py.Dataset
        the solset in an h5parm
    soltab_type : str
        type of the solution table

    Returns
    -------
    str
         frst available name for the soltab
    """

    for soltab_idx in range(1000):
        if f"{soltab_type}{soltab_idx:03d}" not in solset:
            return f"{soltab_type}{soltab_idx:03d}"
    return f"{soltab_type}9999"


def create_solset(h5file: h5py.File, solset_name: str | None = None) -> h5py.Dataset:
    """Create new empty solset in your h5parm file. A solset has an antenna and source table. Default name is sol###"

    Parameters
    ----------
    h5file : h5py.File
        the h5parm file, must be writable
    solset_name : str, optional
        name of the solset, if None the first available from sol### will be used, by default None

    Returns
    -------
    h5py.Dataset
        pointer to the new solset
    """
    if solset_name is None:
        solset_name = get_minimal_solset_name(h5file)
    if solset_name in h5file:
        message = f"{solset_name} already exists in {h5file.name}."
        raise RuntimeError(message)

    solset = h5file.create_group(solset_name)
    solset.attrs.h5parm_version = "1.0"
    solset.attrs["CLASS"] = np.bytes_("GROUP")
    solset.attrs["FILTERS"] = 0
    solset.attrs["TITLE"] = np.bytes_("")
    solset.attrs["VERSION"] = np.bytes_("1.0")
    solset.create_dataset(
        "antenna",
        shape=(1,),
        maxshape=(None,),
        chunks=True,
        dtype=[("name", "S16"), ("position", "<f4", (3,))],
    )
    solset.create_dataset(
        "source",
        shape=(1,),
        maxshape=(None,),
        chunks=True,
        dtype=[("name", "S16"), ("dir", "<f4", (2,))],
    )
    return solset


def add_antenna_info(
    solset: h5py.Dataset, antenna_names: NDArray[Any], antenna_pos: NDArray[Any]
) -> None:
    """Add antenna metadata to a solset

    Parameters
    ----------
    solset : h5py.Dataset
        A solset in an h5parm
    antenna_names : NDArray[Any]
        array with antenna names
    antenna_pos : NDArray[Any]
        array ant x 3 with ITRF antenna positions
    """
    ant_meta = solset["antenna"]
    ant_meta.resize((len(antenna_names),))
    for antenna_number, (antenna_name, pos) in enumerate(
        zip(antenna_names, antenna_pos)
    ):
        ant_meta[antenna_number] = (antenna_name, pos)


def add_source_info(
    solset: h5py.Dataset, source_names: NDArray[Any], source_dirs: NDArray[Any]
) -> None:
    """Add source metadata to a solset

    Parameters
    ----------
    solset : h5py.Dataset
        A solset in an h5parm
    source_names : NDArray[Any]
        array with source names
    source_dirs : NDArray[Any]
        array source x 2 with Ra,DEC directions of the sources
    """
    source_meta = solset["source"]
    source_meta.resize((len(source_names),))
    for source_number, (source_name, source_dir) in enumerate(
        zip(source_names, source_dirs)
    ):
        source_meta[source_number] = (source_name, source_dir)


def add_soltab(
    solset: h5py.Dataset,
    soltab_type: str,
    val: NDArray[Any],
    weight: NDArray[Any],
    soltab_axes: list[Any],
    axes_values: dict[str, Any],
    soltab_name: str | None = None,
) -> h5py.Dataset:
    """Add asolution table to an existing solset

    Parameters
    ----------
    solset : h5py.Dataset
        solset
    soltab_type : str
        type of the solutions
    val : NDArray[Any]
        array with solutions values, must match the shape of the axes
    weight : NDArray[Any]
        array of weights, must match the shape of val
    soltab_axes : list
        ordered list of axes as they are defined for val/wegiht
    axes_values : dict[Any]
        dictionary with the values for each axis in soltab_axes
    soltab_name : str|None, optional
        name of the solution table, if None or existing will default to <soltab_type>###, by default None

    Returns
    -------
    h5py.Dataset
        the solution table
    """
    if soltab_name is None:
        soltab_name = get_minimal_soltab_name(solset, soltab_type)
    if soltab_name in solset:
        message = f"{soltab_name} already exists in {solset.name}."
        raise RuntimeError(message)

    soltab = solset.create_group(soltab_name)
    soltab.attrs["CLASS"] = np.bytes_("GROUP")
    soltab.attrs["FILTERS"] = 0
    soltab.attrs["TITLE"] = np.bytes_(soltab_type)
    soltab.attrs["VERSION"] = np.bytes_("1.0")
    soltab_axes_str = ",".join(soltab_axes)
    for axis in soltab_axes:
        soltab.create_dataset(axis, data=axes_values[axis])

    soltab_val = soltab.create_dataset("val", shape=val.shape, dtype="<f8")
    soltab_val.attrs["soltype"] = soltab_type
    soltab_axes_dtype = _zero_terminated_string(len(soltab_axes_str) + 1)
    soltab_val.attrs["AXES"] = np.array(soltab_axes_str, dtype=soltab_axes_dtype)
    soltab_val[:] = val
    soltab_weight = soltab.create_dataset("weight", shape=weight.shape, dtype="<f2")
    soltab_weight.attrs["AXES"] = np.array(soltab_axes_str, dtype=soltab_axes_dtype)
    soltab_weight[:] = weight  # Todo: np.isnan?
    return soltab


def check_h5parm(h5file: h5py.File) -> bool:
    """Check if given h5file is suitable as h5parm

    Parameters
    ----------
    h5file : h5py.File
        possible h5parm file

    Returns
    -------
    bool
        True if the file is conform h5parm definition
    """
    return "TITLE" in h5file.attrs


def _get_station_metadata(rms: dict[str, Any]) -> tuple[list[Any], list[Any]]:
    """helper function to generate arrays of station names and positions

    Parameters
    ----------
    rms : dict[RM|DTEC]
        dictionary with RM object per station

    Returns
    -------
    tuple
        station_names and station_positions (ITRF xyz)
    """
    station_names = sorted(rms.keys())
    station_pos = [
        [
            rms[st_name].loc.to(u.m).x.value,
            rms[st_name].loc.to(u.m).y.value,
            rms[st_name].loc.to(u.m).z.value,
        ]
        for st_name in station_names
    ]
    return station_names, station_pos


def create_empty_h5parm(h5parm_name: str) -> None:
    """create an empty h5parm file

    Parameters
    ----------
    h5parm_name : str
        name of the file

    Raises
    ------
    TypeError
       error if h5parm_name is existing and not a valid h5parm
    """
    if Path(h5parm_name).is_file():
        if not h5py.is_hdf5(h5parm_name):
            msg = f"{h5parm_name} exists and not a valid hdf5 file"
            raise ValueError(msg)

        with h5py.File(h5parm_name, "r") as h5parm:
            if not check_h5parm(h5parm):
                msg = f"{h5parm_name} exists and not a valid h5parm file"
                raise ValueError(msg)
    else:
        with h5py.File(h5parm_name, "w") as h5parm:
            h5parm.attrs["h5parm_version"] = np.bytes_("1.0")
            h5parm.attrs["CLASS"] = np.bytes_("GROUP")
            h5parm.attrs["FILTERS"] = 0
            h5parm.attrs["TITLE"] = np.bytes_("")
            h5parm.attrs["VERSION"] = np.bytes_("1.0")


def _zero_terminated_string(size: int = 10) -> h5py.Datatype:
    tid = h5py.h5t.C_S1.copy()
    tid.set_size(size)
    return h5py.Datatype(tid)


def write_rm_to_h5parm(
    rms: dict[str, RM],
    h5parm_name: str,
    solset_name: str | None = None,
    soltab_name: str | None = None,
    add_to_existing_solset: bool = False,
) -> None:
    """writes a dictionary of RM values per station to a new or existing h5parm file

    Parameters
    ----------
    rms : dict[RM]
        rm values per station
    h5parm_name : str
        name of the h5parm file
    solset_name : str | None, optional
        name of the solset if None it  will default to sol###, by default None
    soltab_name : str | None, optional
        name of the soltab if None it  will default to 'rotationmeasure###'
    add_to_existing_solset : bool = False
        whether to append to an existing solset, if it exists. If True, the user
        is responsible for having consistent antennas and sources.

    Raises
    ------
    TypeError
        error if h5parm_name is existing and not a valid h5parm
    """
    create_empty_h5parm(h5parm_name=h5parm_name)
    station_names, station_pos = _get_station_metadata(rms)
    with h5py.File(h5parm_name, "a") as h5parm:
        if solset_name is not None and solset_name in h5parm and add_to_existing_solset:
            solset = h5parm[solset_name]
        else:
            solset = create_solset(h5parm, solset_name=solset_name)
            add_antenna_info(solset, station_names, station_pos)

        soltab_axes = ["ant", "time"]
        axes_values = {}
        ant_dtype = _zero_terminated_string(max(map(len, station_names)) + 1)
        axes_values["ant"] = np.array(station_names, dtype=ant_dtype)
        axes_values["time"] = (
            rms[station_names[0]].times.mjd * 24 * 3600.0
        )  # mjd in seconds?
        rm_values = np.array([rms[stname].rm for stname in station_names])
        weights = np.ones(rm_values.shape, dtype=bool)
        weights[np.isnan(rm_values)] = 0
        add_soltab(
            solset=solset,
            soltab_type="rotationmeasure",
            val=rm_values,
            weight=weights,
            soltab_axes=soltab_axes,
            axes_values=axes_values,
            soltab_name=soltab_name,
        )


def write_tec_to_h5parm(
    dtec: dict[str, DTEC],
    h5parm_name: str,
    solset_name: str | None = None,
    soltab_name: str | None = None,
    add_to_existing_solset: bool = False,
) -> None:
    """writes a dictionary of RM values per station to a new or existing h5parm file

    Parameters
    ----------
    dtec : dict[DTEC]
        electron density profiles per station
    h5parm_name : str
        name of the h5parm file
    solset_name : str | None, optional
        name of the solset if None it  will default to sol###, by default None
    soltab_name : str | None, optional
        name of the soltab if None it  will default to 'rotationmeasure###'
    add_to_existing_solset : bool = False
        whether to append to an existing solset, if it exists. If True, the user
        is responsible for having consistent antennas and sources.

    Raises
    ------
    TypeError
        error if h5parm_name is existing and not a valid h5parm
    """
    create_empty_h5parm(h5parm_name=h5parm_name)
    station_names, station_pos = _get_station_metadata(dtec)
    with h5py.File(h5parm_name, "a") as h5parm:
        if solset_name is not None and solset_name in h5parm and add_to_existing_solset:
            solset = h5parm[solset_name]
        else:
            solset = create_solset(h5parm, solset_name=solset_name)
            add_antenna_info(solset, station_names, station_pos)

        soltab_axes = ["ant", "time"]
        axes_values = {}
        ant_dtype = _zero_terminated_string(max(map(len, station_names)) + 1)
        axes_values["ant"] = np.array(station_names, dtype=ant_dtype)
        axes_values["time"] = (
            dtec[station_names[0]].times.mjd * 24 * 3600.0
        )  # mjd in seconds?
        dtec_values = np.array(
            [
                np.sum(dtec[stname].electron_density * dtec[stname].airmass, axis=-1)
                for stname in station_names
            ]
        )
        weights = np.ones(dtec_values.shape, dtype=bool)
        weights[np.isnan(dtec_values)] = 0
        add_soltab(
            solset=solset,
            soltab_type="tec",
            val=dtec_values,
            weight=weights,
            soltab_axes=soltab_axes,
            axes_values=axes_values,
            soltab_name=soltab_name,
        )
