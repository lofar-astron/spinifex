"""Module to get iri density profile"""

from __future__ import annotations

import astropy.units as u
import numpy as np
from numpy.typing import NDArray
from PyIRI import coeff_dir
from PyIRI.main_library import IRI_density_1day

from spinifex.geometry.get_ipp import IPP
from spinifex.times import get_unique_days

IRI_HEIGHTS = np.linspace(50, 20000, 100) * u.km


def get_profile(ipp: IPP) -> NDArray[np.float32]:
    """Get the normalized electron density profile for all times an altitudes in ipp

    Parameters
    ----------
    ipp : IPP
        ionospheric piercepoints

    Returns
    -------
    NDArray
        normalied density profile per time
    """
    unique_days = get_unique_days(times=ipp.times)
    edp = np.zeros(ipp.loc.shape, dtype=float)  # electron density profile
    altitudes = ipp.loc.height.to(u.km).value
    longitudes = ipp.loc.lon.deg
    latitudes = ipp.loc.lat.deg
    f107 = 100
    ccir_or_ursi = 0
    for u_day in unique_days:
        year, month, day, _, _, _ = u_day.ymdhms
        indices = np.floor(ipp.times.mjd) == np.floor(u_day.mjd)
        index_nr = np.where(indices)[0]
        aalt = altitudes[indices]  # iri input array of heights

        hidx = np.argmin(
            np.abs(aalt - 350), axis=1
        )  # use lon lat closest to altitude of 350 km for profile
        alon = longitudes[indices, hidx]  # iri input array of longitudes
        alat = latitudes[indices, hidx]  # iri input array of latitude
        ahr = ipp.times[indices].ymdhms.hour
        # IRI_density_1day treats lon/lat and hr as two separate arrays with independent lengths
        # I do not see another solution than to loop
        for itime, tmidx in zip(range(ahr.shape[0]), index_nr):
            _, _, _, _, _, _, edpi = IRI_density_1day(
                year,
                month,
                day,
                ahr[itime : itime + 1],
                alon[itime : itime + 1],
                alat[itime : itime + 1],
                aalt[itime],
                f107,
                coeff_dir,
                ccir_or_ursi,
            )
            edp[tmidx] = edpi.squeeze()
    edp /= np.sum(edp, axis=1, keepdims=True)
    return edp
