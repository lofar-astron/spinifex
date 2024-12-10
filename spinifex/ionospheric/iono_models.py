"""Several implementations of Ionospheric Models. They all should have the get_density function"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time

import spinifex.ionospheric.iri_density as iri
from spinifex.ionospheric.ionex_manipulation import _read_ionex_stuff


def _read_tomion_stuff() -> np.ndarray[float]:
    raise NotImplementedError  # ionex_data


def get_density_ionex(loc: EarthLocation, times: Time, **kwargs) -> np.ndarray[float]:
    return _read_ionex_stuff(loc=loc, times=times, **kwargs)


def get_density_ionex_single_layer(
    loc: EarthLocation, times: Time, height=350 * u.km, **kwargs
) -> np.ndarray[float]:
    ntimes = loc.shape[1]  # we assume time is second axis
    index = np.argmin(np.abs(loc.height.to(u.km).value - height.to(u.km).value), axis=0)
    single_layer_loc = loc[index, np.arange(ntimes)]
    return _read_ionex_stuff(loc=single_layer_loc, times=times, **kwargs)


def get_density_ionex_iri(
    loc: EarthLocation, times: Time, **kwargs
) -> np.ndarray[float]:
    profile = iri.get_profile(loc, times)
    tec = get_density_ionex_single_layer(
        loc=loc, times=times, **kwargs
    )  # get tec at single altitude
    return tec[np.newaxis] * profile


# def get_density_tomion(loc: EarthLocation, times: Time, **kwargs) -> np.ndarray[float]:
#    return _read_tomion_stuff()
