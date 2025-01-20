# ruff: noqa :T201
from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spinifex import get_rm

"""Example how to use get_rm_from_skycoord"""
times = Time("2020-01-08T01:00:00") + np.arange(10) * 25 * u.min
# create source object
source = SkyCoord(ra=350.85 * u.deg, dec=58.815 * u.deg)
# create Earth Location
lon = 6.367 * u.deg
lat = 52.833 * u.deg
dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
# get rotation measures
rm = get_rm.get_rm_from_skycoord(loc=dwingeloo, times=times, source=source)
# print to screen
rotation_measures = rm.rm
rm_times = rm.times
elevations = rm.elevation
azimuths = rm.azimuth
print("time      RM (rad/lambda^2)      azimuth (deg)      elevation (deg)")
for myrm, tm, az, el in zip(rotation_measures, rm_times, azimuths, elevations):
    print(f"{tm.isot} {myrm} {az:3.2f} {el:2.2f}")
