import numpy as np
import pandas as pd
import glob
import planetmapper
from astropy.io import fits

from . import utils
from . import spx
from . import constants


def trim_spectra(spectra, wavelengths, min_wl=-float("inf"), max_wl=float("inf")):
    """Trim a spectra between two wavelengths
    
    Args:
        spectra (np.ndarray): The spectra data
        wavelengths (np.ndarray): The corresponding wavelengths
        min_wl (float): The wavlength below which to cut out
        max_wl (float): The wavelength above which to cut out
        
    Returns:
        spectra (np.ndarray): The trimmed spectra
        wavelengths (np.ndarray): The trimmed wavelengths"""
    
    mask = (wavelengths > min_wl) & (wavelengths < max_wl)
    return spectra[mask], wavelengths[mask]


def margin_trim(cube, margin_size=3):
    """
    Set the outer spaxels of a cube to NaN.

    Args:
        cube (np.ndarray): The spectral cube to trim with shape (wavelengths, x, y).
        margin_size (int): The size of the margin to trim from each edge.

    Returns:
        np.ndarray: The trimmed spectral cube with NaN values in the margins.
    """
    
    cube[:, :margin_size] = np.nan
    cube[:, :, :margin_size] = np.nan
    cube[:, -margin_size:] = np.nan
    cube[:, :, -margin_size:] = np.nan

    return cube


def downsample_spectra(spectra, wavelengths, num_points):
    """Downsample a spectra to a given number of spectral points. This resamples
    the data on a new, uniformly spaced grid so if there are gaps in the original
    spectra then the number of requested points may be less than num_points, but never more.
    
    Args:
        spectra (np.ndarray): The spectra data
        wavelengths (np.ndarray): The corresponding wavelengths
        num_points (int): The number of points to reduce the spectrum to
        
    Returns:
        spectra (np.ndarray): The trimmed spectra
        wavelengths (np.ndarray): The trimmed wavelengths"""
    
    idx = np.linspace(0, len(wavelengths)-1, num_points, dtype=int)
    return spectra[idx], wavelengths[idx]


def remove_nonlte_emission(spectra, wavelengths, h3p_threshold=0.1, ch4_threshold=0.1):
    """Use a model to chop out any parts of the spectrum where H3+ and CH4 non-LTE
    emission is above a threshold relative to their highest peaks.
    
    Args:
        spectra (np.ndarray): The spectra data
        wavelengths (np.ndarray): The corresponding wavelengths
        h3p_threshold (float): The threshold above which to exclude data due to H3+ emission
        ch4_threshold (float): The threshold above which to exclude data due to CH4 emission
        
    Returns:
        spectra (np.ndarray): The trimmed spectra
        wavelengths (np.ndarray): The trimmed wavelengths"""
    
    if not (2.8 < wavelengths[0] < 5.3 and 2.8 < wavelengths[-1] < 5.3):
        raise NotImplementedError("Non-LTE subtraction only works for G365H data at the moment.")
    
    h3p = pd.read_csv(constants.PATH / "data/misc/H3p_G395H_model.txt", sep=" ", names=["wavelengths", "spectrum"])
    ch4 = pd.read_csv(constants.PATH / "data/misc/CH4_G395H_model.txt", sep=" ", names=["wavelengths", "spectrum"])

    mask = []
    for w in wavelengths:
        i, _ = utils.find_nearest(h3p.wavelengths, w)
        h3p_flag = h3p.spectrum[i] / np.nanmax(h3p.spectrum) > h3p_threshold
        ch4_flag = ch4.spectrum[i] / np.nanmax(ch4.spectrum) > ch4_threshold
        mask.append(h3p_flag or ch4_flag)
    
    mask = ~np.array(mask)
    return spectra[mask], wavelengths[mask]


def multiple_cube_average(cubes):
    return np.nanmean(np.array(cubes), axis=(0,2,3))


def get_observations(pattern, add_errors=True):
    """Get a list of planetmapper.Observation objects from a glob-style file pattern
    
    Args:
        pattern (str): The pattern to use to search for .fits files
        add_errors (bool): Whether to add the error cubes as an attribute (Observation.error)
        
    Returns:
        List[planetmapper.Observation]: The found observations"""
    fps = sorted(glob.glob(pattern))
    obs = [planetmapper.Observation(fp) for fp in fps]
    if add_error_cubes:
        obs = add_error_cubes(obs)
    return obs


def add_error_cubes(observations):
    """Add the error cube from a JWST observation to a list (or a single) planetmapper.Observation object, 
    accessible using the new Observation.error attribute.
    
    Args:
        observations: planetmapper.Observation or List[planetmapper.Observation]: The observation(s) to add errors for
        
    Returns:
        observations: planetmapper.Observation or List[planetmapper.Observation]: The observation(s) with added errors
    """
    for obs in observations:
        fp = obs.path 
        hdul = fits.open(fp)
        obs.error = hdul['ERR'].data
    return observations
