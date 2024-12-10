"""Several implementations of Ionospheric Models. They all should have the get_density function"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.time import Time
from numpy.typing import ArrayLike

import spinifex.ionospheric.iri_density as iri
from spinifex.ionospheric.ionex_manipulation import _read_ionex_stuff


class ModelDesnsityFunction(Protocol):
    """Model density callable"""
    def __call__(self, loc: EarthLocation, times: Time) -> ArrayLike: ...


@dataclass
class IonosphericModels:
    """Names space for different ionospheric ionospheric_models. An ionospheric model should be
    a callable get_density"""

    ionex: ModelDesnsityFunction
    ionex_iri: ModelDesnsityFunction


def _read_tomion_stuff() -> ArrayLike:
    return None  # ionex_data


def get_density_ionex(loc: EarthLocation, times: Time) -> ArrayLike:
    return _read_ionex_stuff(loc=loc, times=times)


def get_density_ionex_single_layer(
    loc: EarthLocation,
    times: Time,
    height: u.Quantity = 350 * u.km,
) -> ArrayLike:
    n_times = loc.shape[1]  # we assume time is second axis
    index = np.argmin(np.abs(loc.height.to(u.km).value - height.to(u.km).value), axis=0)
    single_layer_loc = loc[index, np.arange(n_times)]
    return _read_ionex_stuff(loc=single_layer_loc, times=times)


def get_density_ionex_iri(loc: EarthLocation, times: Time) -> ArrayLike:
    profile = iri.get_profile(loc, times)
    tec = get_density_ionex_single_layer(
        loc=loc, times=times
    )  # get tec at single altitude
    return tec[np.newaxis] * profile


# def get_density_tomion(loc: EarthLocation, times: Time) -> ArrayLike:
#    return _read_tomion_stuff()

ionospheric_models = IonosphericModels(
    ionex=get_density_ionex_single_layer,
    ionex_iri=get_density_ionex_iri,
)


def get_ionosphere(
    loc: EarthLocation,
    time: Time,
    iono_model: ModelDesnsityFunction = ionospheric_models.ionex,
) -> u.Quantity[u.m**-3]:
    """Prints a nice message"""

    return iono_model(loc=loc, times=time)
