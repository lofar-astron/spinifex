"""Several implementations of Ionospheric Models. They all should have the get_density function"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Generic, Protocol, TypeVar

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from numpy.typing import NDArray
from pydantic import ValidationError

import spinifex.ionospheric.iri_density as iri
from spinifex.geometry.get_ipp import IPP
from spinifex.ionospheric.ionex_download import IonexOptions
from spinifex.ionospheric.ionex_manipulation import get_density_ionex
from spinifex.logger import logger

# TODO: add more ionospheric models
IonoOptions = IonexOptions  # Single model for now, but we can add more later e.g. IonexOptions | NewOptionType | AnotherOptionType
O = TypeVar("O", bound=IonoOptions)  # noqa: E741
O_contra = TypeVar("O_contra", contravariant=True, bound=IonoOptions)


class ModelDensityFunction(Protocol, Generic[O_contra]):
    """Model density callable"""

    def __call__(
        self,
        ipp: IPP,
        height: u.Quantity = 350 * u.km,
        options: O_contra | None = None,
    ) -> NDArray[np.float64]: ...


@dataclass
class IonosphericModels:
    """
    Names space for different ionospheric ionospheric_models. An ionospheric model should be
        a callable get_density
    """

    ionex: ModelDensityFunction[IonexOptions]
    ionex_iri: ModelDensityFunction[IonexOptions]


def get_density_ionex_single_layer(
    ipp: IPP, height: u.Quantity = 350 * u.km, options: IonexOptions | None = None
) -> NDArray[np.float64]:
    """gets the ionex files and interpolate values for a single altitude, thin screen assumption

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    height : u.Quantity, optional
        altitude of the thin scrreen, by default 350*u.km
    ionex_options: IonexOptions | None, optional
        options for the ionospheric model, by default None

    Returns
    -------
    NDArray
        interpolated vTEC values at ipp, zeros everywhere apart from the altitude
        closest to the specified height
    """
    n_times = ipp.times.shape[0]  # we assume time is first axis
    index = np.argmin(
        np.abs(ipp.loc.height.to(u.km).value - height.to(u.km).value), axis=1
    )
    single_layer_loc = EarthLocation(ipp.loc[np.arange(n_times), index])
    ipp_single_layer = IPP(
        loc=single_layer_loc,
        times=ipp.times,
        los=ipp.los,
        airmass=ipp.airmass[:, index],
        altaz=ipp.altaz,
        station_loc=ipp.station_loc,
    )
    result = np.zeros(ipp.loc.shape, dtype=float)
    result[np.arange(n_times), index] = get_density_ionex(
        ipp_single_layer,
        ionex_options=options,
    )
    return result


def get_density_ionex_iri(
    ipp: IPP,
    height: u.Quantity = 350 * u.km,
    options: IonexOptions | None = None,
) -> NDArray[np.float64]:
    profile = iri.get_profile(ipp)
    tec = get_density_ionex_single_layer(ipp, height=height, options=options)
    # get tec at single altitude
    return np.array(np.sum(tec, keepdims=True, axis=1) * profile)


ionospheric_models = IonosphericModels(
    ionex=get_density_ionex_single_layer,
    ionex_iri=get_density_ionex_iri,
)


# TODO: Use the TypeVar `O` when we have more than one ionospheric model
# i.e.
# >>> def parse_iono_kwargs(
# >>>     iono_model: ModelDensityFunction[O], **kwargs: Any
# >>> ) -> O | None:
def parse_iono_kwargs(
    iono_model: ModelDensityFunction[IonoOptions], **kwargs: Any
) -> IonoOptions:
    """parse ionospheric options

    Parameters
    ----------
    iono_model : ModelDensityFunction
        ionospheric model
    **kwargs : Any
        options for the ionospheric model

    Raises
    ------
    TypeError
        Incorrect arguments for ionospheric model

    Returns
    -------
    IonoOptions
        ionospheric model options

    """

    try:
        if iono_model in (ionospheric_models.ionex, ionospheric_models.ionex_iri):
            options = IonexOptions(**kwargs)
            logger.info(f"Using ionospheric model {iono_model} with options {options}")
            return options

        # TODO: add more ionospheric models and their options
        msg = f"Unknown ionospheric model {iono_model}."
        raise TypeError(msg)

    except ValidationError as e:
        msg = f"Incorrect arguments {kwargs} for ionospheric model {iono_model}"
        raise TypeError(msg) from e


def parse_iono_model(iono_model_name: str) -> ModelDensityFunction[IonoOptions]:
    """parse ionospheric model name

    Parameters
    ----------
    iono_model_name : str
        name of the ionospheric model

    Returns
    -------
    ModelDensityFunction
        ionospheric model

    Raises
    ------
    TypeError
        if the ionospheric model is not known

    """

    try:
        return getattr(ionospheric_models, iono_model_name)  # type: ignore[no-any-return]
    except AttributeError as e:
        msg = f"Unknown ionospheric model {iono_model_name}. Supported models are {list(ionospheric_models.__annotations__.keys())}"
        raise TypeError(msg) from e
