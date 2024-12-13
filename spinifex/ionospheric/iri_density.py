"""Module to get iri density profile"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from numpy.typing import ArrayLike
from PyIRI import coeff_dir
from PyIRI.main_library import IRI_density_1day

from spinifex.geometry.get_ipp import IPP


def get_profile(ipp: IPP) -> ArrayLike:
    """Get the normalized electron density profile for all times an altitudes in ipp

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints

    Returns
    -------
    ArrayLike
        normalied density profile per time
    """
    loc = ipp.loc
    times = ipp.times
    aalt = loc[:, 0].height.to(u.km).value  # iri input array of heights
    hidx = np.argmin(
        np.abs(aalt - 350)
    )  # use lon lat closest to altitude of 350 km for profile
    alon = loc[hidx].lon.deg  # iri input array of longitudes
    alat = loc[hidx].lat.deg  # iri input array of latitudes
    f107 = 100
    ccir_or_ursi = 0
    year, month, day, _, _, _ = times[0].ymdhms
    ahr = times.ymdhms.hour
    # IRI_density_1day treats lon/lat and hr as two separate arrays with independent lengths
    # I do not see another solution than to loop
    edp = np.zeros((aalt.shape[0], ahr.shape[0]))
    for itime in range(ahr.shape[0]):
        _, _, _, _, _, _, edpi = IRI_density_1day(
            year,
            month,
            day,
            ahr[itime : itime + 1],
            alon[itime : itime + 1],
            alat[itime : itime + 1],
            aalt,
            f107,
            coeff_dir,
            ccir_or_ursi,
        )
        edp[:, itime] = edpi.squeeze()
    edp /= np.sum(edp, axis=0, keepdims=True)
    return edp
