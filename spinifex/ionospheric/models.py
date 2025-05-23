"""Several implementations of Ionospheric Models. They all should have the get_density function"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Generic, Protocol, TypeVar, cast

import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from numpy.typing import NDArray
from pydantic import ValidationError

import spinifex.ionospheric.iri_density as iri
from spinifex.geometry.get_ipp import IPP
from spinifex.ionospheric.ionex_manipulation import get_density_ionex
from spinifex.ionospheric.tec_data import ElectronDensity, IonexOptions, TomionOptions
from spinifex.ionospheric.tomion_parser import get_density_dual_layer
from spinifex.logger import logger

O = TypeVar("O", IonexOptions, TomionOptions)  # noqa: E741
O_contra = TypeVar("O_contra", IonexOptions, TomionOptions, contravariant=True)


class ModelDensityFunction(Protocol, Generic[O_contra]):
    """Model density callable"""

    def __call__(
        self,
        ipp: IPP,
        options: O_contra | None = None,
    ) -> ElectronDensity: ...


@dataclass
class IonosphericModels:
    """
    Names space for different ionospheric ionospheric_models. An ionospheric model should be
        a callable get_density
    """

    ionex: ModelDensityFunction[IonexOptions]
    ionex_iri: ModelDensityFunction[IonexOptions]
    tomion: ModelDensityFunction[TomionOptions]


def get_density_ionex_single_layer(
    ipp: IPP, options: IonexOptions | None = None
) -> ElectronDensity:
    """gets the ionex files and interpolate values for a single altitude, thin screen assumption

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    ionex_options: IonexOptions | None, optional
        options for the ionospheric model, by default None

    Returns
    -------
    NDArray
        interpolated vTEC values at ipp, zeros everywhere apart from the altitude
        closest to the specified height
    """
    if options is None:
        options = IonexOptions()

    n_times = ipp.times.shape[0]  # we assume time is first axis
    index = np.argmin(
        np.abs(ipp.loc.height.to(u.km).value - options.height.to(u.km).value), axis=1
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
    tec = get_density_ionex(
        ipp_single_layer,
        ionex_options=options,
    )
    electron_density = np.zeros(ipp.loc.shape, dtype=float)
    electron_density[np.arange(n_times), index] = tec.electron_density
    return ElectronDensity(
        electron_density=electron_density,
        electron_density_error=tec.electron_density_error.reshape((-1, 1)),
    )


def get_density_ionex_iri(
    ipp: IPP,
    options: IonexOptions | None = None,
) -> ElectronDensity:
    """gets the ionex files and interpolate values for a single altitude, then multiply with a
    normalised density profile from iri

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints
    height : u.Quantity, optional
        altitude of the thin screen, by default 350*u.km
    ionex_options: IonexOptions | None, optional
        options for the ionospheric model, by default None

    Returns
    -------
    NDArray
        interpolated vTEC values at ipp
    """
    profile = iri.get_profile(ipp)
    tec = get_density_ionex_single_layer(ipp, options=options)
    # get tec at single altitude
    return ElectronDensity(
        electron_density=np.sum(tec.electron_density, keepdims=True, axis=1) * profile,
        electron_density_error=tec.electron_density_error,  # only save the rms  of the single layer
    )


# TODO: move height to IonexOptions
def get_density_tomion(
    ipp: IPP, options: TomionOptions | None = None
) -> NDArray[np.float64]:
    tec = get_density_dual_layer(ipp, tomion_options=options)
    return ElectronDensity(
        electron_density=tec.electron_density,
        electron_density_error=tec.electron_density_error,
    )


ionospheric_models = IonosphericModels(
    ionex=get_density_ionex_single_layer,
    ionex_iri=get_density_ionex_iri,
    tomion=get_density_tomion,
)


def parse_iono_kwargs(iono_model: ModelDensityFunction[O], **kwargs: Any) -> O:
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
        if iono_model in (
            ionospheric_models.ionex,
            ionospheric_models.ionex_iri,
        ):
            ionex_options = IonexOptions(**kwargs)
            logger.info(
                f"Using ionospheric model {iono_model} with options {ionex_options}"
            )
            return cast(O, ionex_options)  # type: ignore[redundant-cast]
        if iono_model == ionospheric_models.tomion:
            tomion_options = TomionOptions(**kwargs)
            logger.info(
                f"Using ionospheric model {iono_model} with options {tomion_options}"
            )
            return cast(O, tomion_options)  # type: ignore[redundant-cast]

        msg = f"Unknown ionospheric model {iono_model}."
        raise TypeError(msg)

    except ValidationError as e:
        msg = f"Incorrect arguments {kwargs} for ionospheric model {iono_model}"
        raise TypeError(msg) from e


def parse_iono_model(
    iono_model_name: str,
) -> ModelDensityFunction[O]:
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
