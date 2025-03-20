from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import EarthLocation
from astropy.io import fits
from spinifex.image_tools import (
    get_integrated_rm_from_fits,
    get_metadata_from_fits,
    get_rm_from_fits,
)
from spinifex.image_tools.fits_tools import get_freq_from_fits


def _make_fits_cube(
    extra_cards: dict[str, str] | None = None,
) -> fits.HDUList:
    header = fits.Header()
    header["NAXIS"] = 4
    header["NAXIS1"] = 10
    header["NAXIS2"] = 10
    header["NAXIS3"] = 6
    header["NAXIS4"] = 1
    header["CTYPE1"] = "RA---SIN"
    header["CRVAL1"] = 0
    header["CRPIX1"] = 5
    header["CDELT1"] = -(1 / 3600)
    header["CUNIT1"] = "deg"
    header["CTYPE2"] = "DEC--SIN"
    header["CRVAL2"] = 0
    header["CRPIX2"] = 5
    header["CDELT2"] = 1 / 3600
    header["CUNIT2"] = "deg"
    header["CTYPE3"] = "FREQ"
    header["CRVAL3"] = 1.4e9
    header["CRPIX3"] = 1
    header["CDELT3"] = 1e6
    header["CUNIT3"] = "Hz"
    header["CTYPE4"] = "STOKES"
    header["CRVAL4"] = 1
    header["CRPIX4"] = 1
    header["CDELT4"] = 1
    header["CUNIT4"] = ""
    header["BUNIT"] = "Jy/beam"
    header["DATE-OBS"] = (
        "2019-04-25T12:45:52.893302"  # pulled from a random SPICE-RACS file
    )

    if extra_cards is not None:
        for key, value in extra_cards.items():
            header[key] = value

    data = np.ones((1, 6, 10, 10), dtype=np.float32)

    hdu = fits.PrimaryHDU(data, header)

    return fits.HDUList([hdu])


@pytest.fixture
def simple_fits_cube(tmpdir):
    """Create a simple FITS cube"""

    simple_hdul = _make_fits_cube(
        {
            "DURATION": "1800",
            "TELESCOP": "ASKAP",
        }
    )

    fits_path = Path(tmpdir) / "simple.fits"
    simple_hdul.writeto(fits_path, output_verify="silentfix", overwrite=True)

    yield fits_path

    fits_path.unlink()


@pytest.mark.filterwarnings("ignore:'datfix' made the change")
def test_get_freq_from_fits(simple_fits_cube):
    freqs = get_freq_from_fits(fits_path=simple_fits_cube)
    freqs_mhz = freqs.to("MHz").value
    test_vals = np.arange(1400, 1400 + 6, 1, dtype=np.float32)
    assert np.array_equal(freqs_mhz, test_vals)


@pytest.mark.filterwarnings("ignore:'datfix' made the change")
def test_get_metadata_from_fits(simple_fits_cube):
    simple_metadata = get_metadata_from_fits(fits_path=simple_fits_cube)
    assert simple_metadata.duration.sec == 1800
    assert simple_metadata.start_time.fits == "2019-04-25T12:45:52.893"
    assert simple_metadata.name == "simple.fits"
    assert simple_metadata.location == EarthLocation.of_site("ASKAP")
    assert np.isclose(
        simple_metadata.source.ra.wrap_at(180 * u.deg).deg,
        0,
        atol=1 / 3600 * 2,
    )
    assert np.isclose(
        simple_metadata.source.dec.deg,
        0,
        atol=1 / 3600 * 2,
    )


@pytest.mark.filterwarnings("ignore:'datfix' made the change")
def test_get_rm_from_fits(simple_fits_cube, tmpdir):
    rm = get_rm_from_fits(
        fits_path=simple_fits_cube,
        output_directory=Path(tmpdir),
    )
    expected_rm = np.array([-0.69612298, -0.6634822])
    assert np.allclose(rm.rm, expected_rm)


@pytest.mark.filterwarnings("ignore:'datfix' made the change")
def test_get_integrated_rm_from_fits(simple_fits_cube, tmpdir):
    integrated_rm = get_integrated_rm_from_fits(
        fits_path=simple_fits_cube,
        output_directory=Path(tmpdir),
    )
    rm = get_rm_from_fits(
        fits_path=simple_fits_cube,
        output_directory=Path(tmpdir),
    )
    delta_angles = np.rad2deg(np.angle(integrated_rm.theta))

    freqs = get_freq_from_fits(simple_fits_cube)
    lsq = (speed_of_light / freqs).to("m") ** 2

    # Compute basic angle change from sum of angle change per timestep
    basic_angles = (
        (rm.rm[..., np.newaxis] * u.rad / u.m**2 * lsq).to("deg").sum(axis=0).value
    )

    assert np.allclose(delta_angles, basic_angles, atol=1e-1)
