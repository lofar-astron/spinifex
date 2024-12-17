"""Several implementations of Ionospheric Models. They all should have the get_density function"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import astropy.units as u
import numpy as np
from numpy.typing import ArrayLike

import spinifex.ionospheric.iri_density as iri
from spinifex.geometry.get_ipp import IPP
from spinifex.ionospheric.ionex_manipulation import _read_ionex_stuff


class ModelDensityFunction(Protocol):
    """Model density callable"""

    def __call__(self, ipp: IPP, iono_kwargs: dict | None = None) -> ArrayLike: ...


@dataclass
class IonosphericModels:
    """Names space for different ionospheric ionospheric_models. An ionospheric model should be
    a callable get_density"""

    ionex: ModelDensityFunction
    ionex_iri: ModelDensityFunction


def _read_tomion_stuff() -> ArrayLike:
    raise NotImplementedError


def get_density_ionex(ipp: IPP) -> ArrayLike:
    return _read_ionex_stuff(ipp)


def get_density_ionex_single_layer(
    ipp: IPP, height: u.Quantity = 350 * u.km, iono_kwargs: dict | None = None
) -> ArrayLike:
    """gets the ionex files and interpolate values for a single altitude, thin screen assumption

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    height : u.Quantity, optional
        altitude of the thin scrreen, by default 350*u.km
    iono_kwargs: dict, optional

    Returns
    -------
    ArrayLike
        interpolated vTEC values at ipp, zeros everywhere apart from the altitude
        closest to the specified height
    """
    iono_kwargs = iono_kwargs or {}
    n_times = ipp.times.shape[0]  # we assume time is first axis
    index = np.argmin(
        np.abs(ipp.loc.height.to(u.km).value - height.to(u.km).value), axis=1
    )
    single_layer_loc = ipp.loc[np.arange(n_times), index]
    ipp_single_layer = IPP(
        loc=single_layer_loc,
        times=ipp.times,
        los=ipp.los,
        airmass=ipp.airmass[:, index],
        altaz=ipp.altaz,
    )
    result = np.zeros(ipp.loc.shape, dtype=float)
    result[np.arange(n_times), index] = _read_ionex_stuff(
        ipp_single_layer, iono_kwargs=iono_kwargs
    )
    return result


def get_density_ionex_iri(ipp: IPP, iono_kwargs: dict | None = None) -> ArrayLike:
    iono_kwargs = iono_kwargs or {}
    profile = iri.get_profile(ipp)
    tec = get_density_ionex_single_layer(
        ipp, height=350 * u.km, iono_kwargs=iono_kwargs
    )
    # get tec at single altitude
    return np.sum(tec, keepdims=True, axis=1) * profile


# def get_density_tomion(ipp: IPP) -> ArrayLike:
#    return _read_tomion_stuff()

ionospheric_models = IonosphericModels(
    ionex=get_density_ionex_single_layer,
    ionex_iri=get_density_ionex_iri,
)


def get_ionosphere(
    ipp: IPP,
    iono_model: ModelDensityFunction = ionospheric_models.ionex,
) -> u.Quantity[u.m**-3]:
    """Prints a nice message"""

    return iono_model(ipp)
