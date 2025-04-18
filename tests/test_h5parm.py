#  Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
#  SPDX-License-Identifier: Apache-2.0
"""Testing of the main module"""

from __future__ import annotations

from astropy.utils import iers

iers.conf.auto_download = False

import astropy.units as u
import h5py
import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from spinifex import h5parm_tools as h5pt
from spinifex.get_rm import RM
from spinifex.ionospheric.tec_data import ElectronDensity
from spinifex.magnetic.models import MagneticProfile


def test_create_h5parm(tmpdir):
    """test that h5parm is created"""
    h5parm_name = f"{tmpdir}/test.h5"
    h5pt.create_empty_h5parm(h5parm_name)
    with h5py.File(h5parm_name, "a") as h5parm:
        solset = h5pt.create_solset(h5parm, "mysolset")
        assert "mysolset" in h5parm
        station_names = ["station1"]
        station_pos = [[0, 0, 0]]
        h5pt.add_antenna_info(
            solset=solset, antenna_names=station_names, antenna_pos=station_pos
        )
        assert solset["antenna"]["name"][0].decode() == "station1"
        source_names = ["source1"]
        source_dir = [[0, 0]]
        h5pt.add_source_info(
            solset=solset, source_names=source_names, source_dirs=source_dir
        )
        assert solset["source"]["name"][0].decode() == "source1"
        axes_values = {}
        axes_values["axis1"] = [1, 2]
        axes_values["axis2"] = [3, 4, 5]
        soltab_axes = ["axis1", "axis2"]
        val = np.ones((2, 3))
        weight = np.ones(val.shape, dtype=bool)
        h5pt.add_soltab(
            solset=solset,
            soltab_type="rotationmeasure",
            val=val,
            weight=weight,
            soltab_axes=soltab_axes,
            axes_values=axes_values,
        )
        assert "rotationmeasure000" in solset
        assert solset["rotationmeasure000"]["val"].shape == (2, 3)


def test_write_rm_h5parm(tmpdir):
    """test that writing rotation measures does not crash"""
    h5parm_name = f"{tmpdir}/test.h5"
    rms = {}
    lon = 6.367 * u.deg
    lat = 52.833 * u.deg
    dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
    rm = RM(
        rm=np.arange(0, 1, 0.1),
        rm_error=np.arange(0, 1, 0.1) * 0.1,
        times=Time("2025-02-08") + np.arange(10) * u.hr,
        b_parallel=MagneticProfile(
            magnetic_field=np.zeros((10, 1)), magnetic_field_error=np.zeros((10, 1))
        ),
        electron_density=ElectronDensity(
            electron_density=np.zeros((10, 1)), electron_density_error=np.zeros((10, 1))
        ),
        height=np.array([350.0]),
        azimuth=np.zeros((10,)),
        elevation=np.zeros((10,)),
        loc=dwingeloo,
    )
    stations = ["CS001LBA", "CS002LBA"]
    for stat in stations:
        rms[stat] = rm
    h5pt.write_rm_to_h5parm(rms, h5parm_name=h5parm_name)
    # check sol000 exists, check sol000/rotationmeasure000 exists
    with h5py.File(h5parm_name) as myh5:
        assert "sol000" in myh5
        assert "sol000/rotationmeasure000" in myh5
    h5pt.write_rm_to_h5parm(rms, h5parm_name=h5parm_name)
    # check sol001 exists, check sol001/rotationmeasure000 exists
    with h5py.File(h5parm_name) as myh5:
        assert "sol001" in myh5
        assert "sol001/rotationmeasure000" in myh5
    with pytest.raises(RuntimeError, match="sol000 already exists"):
        h5pt.write_rm_to_h5parm(rms, h5parm_name=h5parm_name, solset_name="sol000")
        # check runtime error
    h5pt.write_rm_to_h5parm(
        rms, h5parm_name=h5parm_name, solset_name="sol000", add_to_existing_solset=True
    )
    # check sol000 exists, check sol000/rotationmeasure001 exists
    with h5py.File(h5parm_name) as myh5:
        assert "sol000" in myh5
        assert "sol000/rotationmeasure001" in myh5
    with pytest.raises(RuntimeError, match="rotationmeasure000 already exists"):
        h5pt.write_rm_to_h5parm(
            rms,
            h5parm_name=h5parm_name,
            solset_name="sol000",
            soltab_name="rotationmeasure000",
            add_to_existing_solset=True,
        )
        # check runtime error
    h5pt.write_rm_to_h5parm(
        rms,
        h5parm_name=h5parm_name,
        solset_name="sol000",
        soltab_name="RMextract",
        add_to_existing_solset=True,
    )
    # check sol000 exists, check sol000/RMextract exists
    with h5py.File(h5parm_name) as myh5:
        assert "sol000" in myh5
        assert "sol000/RMextract" in myh5
    h5pt.write_rm_to_h5parm(rms, h5parm_name=h5parm_name, solset_name="sol042")
    # check sol042 exists, check sol042/rotationmeasure000 exists
    with h5py.File(h5parm_name) as myh5:
        assert "sol042" in myh5
        assert "sol042/rotationmeasure000" in myh5
    h5pt.write_rm_to_h5parm(
        rms,
        h5parm_name=h5parm_name,
        soltab_name="RMextract",
        add_to_existing_solset=True,
    )
    # check sol002 exists, check sol002/RMextract exists
    with h5py.File(h5parm_name) as myh5:
        assert "sol002" in myh5
        assert "sol002/RMextract" in myh5
