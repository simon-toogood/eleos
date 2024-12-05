from pathlib import Path
from typing import NamedTuple
import astropy.units as u
import numpy as np


class SpxGeometry(NamedTuple):
    nav: int
    lat: float
    lon: float
    phase: float
    emission: float
    azimuth: float
    wgeom: float
    wavelengths: np.ndarray
    spectrum: np.ndarray
    error: np.ndarray


class SpxFileData(NamedTuple):
    fwhm: float
    lat: float
    lon: float
    geometries: tuple[SpxGeometry, ...]


def read(path: str) -> SpxFileData:
    """
    Read in a .spx file from disk.

    The returned data is a SpxFileData object containing the data from the .spx file.
    For example, you can access the latitude of the observation with `data.lat` or the
    spectrum of the first geometry with `data.geometries[0].spectrum`.

    Args:
        path: The path to the .spx file.

    Returns:
        A SpxFileData object of the form `(fwhm, lat, lon, geometries)` containing the
        data from the .spx file.
    """
    with open(path, 'r') as f:
        fwhm, hdr_lat, hdr_lon, ngeom = (float(x) for x in f.readline().split())
        ngeom = int(ngeom)
        geometries: list[SpxGeometry] = []
        for _ in range(ngeom):
            nconv = int(float(f.readline()))
            nav = int(float(f.readline()))
            lat, lon, phase, emission, azimuth, wgeom = (
                float(x) for x in f.readline().split()
            )
            wavelengths = np.full(nconv, np.nan)
            spectrum = np.full(nconv, np.nan)
            error = np.full(nconv, np.nan)
            for i in range(nconv):
                w, s, e = (float(x) for x in f.readline().split())
                wavelengths[i] = w
                spectrum[i] = s
                error[i] = e
            geometries.append(
                SpxGeometry(
                    nav=nav,
                    lat=lat,
                    lon=lon,
                    phase=phase,
                    emission=emission,
                    azimuth=azimuth,
                    wgeom=wgeom,
                    wavelengths=wavelengths,
                    spectrum=spectrum,
                    error=error,
                )
            )
    return SpxFileData(
        fwhm=fwhm,
        lat=hdr_lat,
        lon=hdr_lon,
        geometries=tuple(geometries),
    )


def convert_MJysr_to_Wcm2srum(
    wavelengths: np.ndarray,
    spectrum: np.ndarray,
) -> np.ndarray:
    """
    Convert a spectrum or cube from MJy/sr to W/cm2/sr/micron.
    """
    while len(wavelengths.shape) < len(spectrum.shape):
        # work for cubes passed to sp
        wavelengths = np.expand_dims(wavelengths, -1)

    spx_MJysr = spectrum * u.MJy / u.sr
    spx_Wcm2srum = spx_MJysr.to(
        u.W / (u.cm * u.cm) / u.sr / u.micron,
        equivalencies=u.spectral_density(wavelengths * u.micron),
    )

    return spx_Wcm2srum.value


def convert_Wcm2srum_to_MJysr(
    wavelengths: np.ndarray,
    spectrum: np.ndarray,
) -> np.ndarray:
    """
    Convert a spectrum or cube from W/cm2/sr/micron to MJy/sr.
    """
    while len(wavelengths.shape) < len(spectrum.shape):
        # work for cubes passed to sp
        wavelengths = np.expand_dims(wavelengths, -1)

    spx_Wcm2srum = spectrum * u.W / (u.cm * u.cm) / u.sr / u.micron
    spx_MJysr = spx_Wcm2srum.to(
        u.MJy / u.sr,
        equivalencies=u.spectral_density(wavelengths * u.micron),
    )
    return spx_MJysr.value


    """Write data (in MJy/sr) to a .spx file (in W/cm2/sr/um)."""
    converted_spectrum = convert_MJysr_to_Wcm2srum(wavelengths, spectrum)
    converted_error = convert_MJysr_to_Wcm2srum(wavelengths, error)
    spx_data = construct_spx(wavelengths, converted_spectrum, converted_error, **kwargs)
    write(filepath, spx_data)