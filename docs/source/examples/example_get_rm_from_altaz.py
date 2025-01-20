from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, EarthLocation
from astropy.time import Time
from spinifex import get_rm

"""Example how to use get_rm_from_altaz"""
times = Time("2020-01-08T01:00:00")
# create AltAz grid
azangles = np.linspace(0, 360, 30)
altangles = np.linspace(10, 90, 8)
altazangles = np.meshgrid(azangles, altangles)
az = altazangles[0] * u.deg
alt = altazangles[1] * u.deg
# create Earth Location
lon = 6.367 * u.deg
lat = 52.833 * u.deg
dwingeloo = EarthLocation(lon=lon, lat=lat, height=0 * u.km)
# make altaz object, including location and times
# get_rm expects a flattened array of coordinates, you can reshape the returned results
altaz = AltAz(az=az.flatten(), alt=alt.flatten(), location=dwingeloo, obstime=times)
# get rotation measures
rm = get_rm.get_rm_from_altaz(loc=dwingeloo, altaz=altaz)
rotation_measures = rm.rm.reshape(az.shape)
rm_times = rm.times
