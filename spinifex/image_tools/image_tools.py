"""Predict / correct ionospheric RM in images"""

# Follows report by Van Eck (2021) and
# Python implementation in FRion https://github.com/CIRADA-Tools/FRion/
# FRion License: MIT Copyright (c) 2021 Cameron Van Eck

from __future__ import annotations

from typing import NamedTuple

import astropy.units as u
import numpy as np
from astropy.constants import c as speed_of_light
from astropy.time import Time
from numpy.typing import NDArray
from scipy import integrate

from spinifex.get_rm import RM, get_average_rm
from spinifex.logger import logger


class IntegratedRM(NamedTuple):
    """Time-integrated rotation  measure"""

    theta: NDArray[np.complex128]
    """Complex time-integrated effect of the ionosphere"""
    time: Time
    """Average time"""
    b_parallel: float
    """Average parallel magnetic field"""
    electron_density: float
    """Average electron content"""
    height: float
    """Average altitude (km)"""
    azimuth: float
    """Average azimuth (degrees)"""
    elevation: float
    """Average elevation (degrees)"""


def get_integrated_rm(
    time_dep_rm: RM,
    freq_arr: u.Quantity,
) -> IntegratedRM:
    """Computed the integrated RM effect following Van Eck (2021)

    Parameters
    ----------
    time_dep_rm : RM
        Time-dependent RM object
    freq_arr : u.Quantity
        Frequency array

    Returns
    -------
    IntegratedRM
        Integrated rotation measure object
    """

    # Get required values for integration
    duration_seconds = (
        ((time_dep_rm.times.max() - time_dep_rm.times.min()).jd * u.day)
        .to(u.second)
        .value
    )
    times_mjd_seconds = (time_dep_rm.times.mjd * u.day).to(u.second).value

    wave_sq = (u.Quantity(speed_of_light) / freq_arr) ** 2
    rm_arr = time_dep_rm.rm * u.rad / u.m**2

    # Following Van Eck (2021) Eqn. 3 - compute integrated RM
    integrand = np.exp(
        2j * (wave_sq[:, np.newaxis] * rm_arr[np.newaxis]).to(u.rad).value
    )
    theta = (1 / (duration_seconds)) * integrate.simpson(
        y=integrand, x=times_mjd_seconds, axis=1
    )

    assert theta.shape == freq_arr.shape, (
        f"Wrong shape for theta ({theta.shape}), expected shape of {freq_arr.shape}"
    )

    # Following FRion implementation (https://github.com/CIRADA-Tools/FRion/)
    # report where numerical issues are likely to cause problems

    longest_wave_sq = np.max(wave_sq).to(u.m**2)
    largest_delta_rm = np.max(np.abs(np.diff(rm_arr)))
    largest_delta_pa = largest_delta_rm * longest_wave_sq

    if largest_delta_pa > 10 * u.deg:
        # Report warnings if the largest change in RM is greater than 10 degrees
        # This a slightly arbitrary threshold - FRion uses 0.5 rad ~ 28.6 degrees
        logger.warning(
            "Largest change of RM in time is greater than 10 degrees, "
            "this may cause numerical issues."
            "Consider using a smaller timestep."
        )

    for bad_fraction in (
        # Report warnings if the depolarisation is greater than 2% or 10%
        0.02,
        0.1,
    ):
        if np.min(np.abs(theta)) < bad_fraction:
            percentage = (1 - bad_fraction) * 100
            msg = (
                f"Depolarsation greater than {percentage:0.1f}% predicted!\n"
                "Images corrected by these values will likely be unreliable."
            )
            logger.warning(msg)
    averaged_rm = get_average_rm(time_dep_rm)
    return IntegratedRM(
        theta=theta,
        time=averaged_rm.times,
        b_parallel=averaged_rm.b_parallel,
        electron_density=averaged_rm.electron_density,
        height=averaged_rm.height,
        azimuth=averaged_rm.azimuth,
        elevation=averaged_rm.elevation,
    )
