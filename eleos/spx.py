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


def write(
    path: str | Path,
    data: SpxFileData | None = None,
    **kwargs,
) -> None:
    """
    Write an SPX file to disk, for example:

        spx.write(
            'spectrum.spx',
            wavelengths=wavelengths,
            spectrum=spectrum,  # In MJy/sr - the unit conversion to W/cm2/sr/um is automatic
            error=error,        # In MJy/sr - the unit conversion to W/cm2/sr/um is automatic
            fwhm=0,
            lat=20.0,
            lon=30.0,
            phase=0.0,
            emission=0.0,
            azimuth=0.0,
        )

    Args:
        path: The path to write the SPX file to.
        data: An SpxFileData object to write to disk. If `None`, the data will be
            constructed from the keyword arguments.
        kwargs: Keyword arguments to construct the SpxFileData object if `data` is
            `None`. These are passed to `construct_spx`.
    """
    if data is None:
        data = construct_spx(**kwargs)

    lines: list[str] = []
    ngeom = len(data.geometries)
    lines.append(f'{data.fwhm}\t{data.lat}\t{data.lon}\t{ngeom}')
    for g in data.geometries:
        mask = ~(np.isnan(g.spectrum) | np.isnan(g.error))
        nconv = len(g.wavelengths[mask])
        lines.append(f'{nconv}')
        lines.append(f'{g.nav}')
        lines.append(
            f'{g.lat}\t{g.lon}\t{g.phase}\t{g.emission}\t{g.azimuth}\t{g.wgeom}'
        )
        for w, s, e in zip(g.wavelengths[mask], g.spectrum[mask], g.error[mask]):
            lines.append(f'{w}\t{s}\t{e}')
    lines.append('')  # end file with a newline
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def construct_spx(
    wavelengths: np.ndarray,
    spectrum: np.ndarray,
    error: np.ndarray,
    *,
    fwhm: float,
    lat: float,
    lon: float,
    phase: float,
    emission: float,
    azimuth: float,
) -> SpxFileData:
    """
    Construct an SpxFileData object from the given data.

    The units should be MJy/sr - the unit conversion is done automatically.
    """
    return SpxFileData(
        fwhm=fwhm,
        lat=lat,
        lon=lon,
        geometries=(
            SpxGeometry(
                nav=1,
                lat=lat,
                lon=lon,
                phase=phase,
                emission=emission,
                azimuth=azimuth,
                wgeom=1,
                wavelengths=wavelengths,
                spectrum=convert_MJysr_to_Wcm2srum(wavelengths, spectrum),
                error=convert_MJysr_to_Wcm2srum(wavelengths, error),
            ),
        ),
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