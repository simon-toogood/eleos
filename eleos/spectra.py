import numpy as np
import pandas as pd
import glob
import h3ppy
import planetmapper
import scipy
import copy
from astropy.io import fits
from collections import defaultdict

from . import utils
from . import spx
from . import constants
from . import ch4


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

    if margin_size == 0:
        return cube
    
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


def get_nonlte_emission_mask(wavelengths, 
                             h3p_threshold=0.1, 
                             ch4_threshold=0.1, 
                             h3p_expand=0, 
                             ch4_expand=0):
    """Use a model to mask out any parts of the spectrum where H3+ and CH4 non-LTE
    emission is above a threshold relative to their highest peaks.
    
    Args:
        spectra (np.ndarray): The spectra data
        wavelengths (np.ndarray): The corresponding wavelengths
        h3p_threshold (float): The threshold above which to exclude data due to H3+ emission
        ch4_threshold (float): The threshold above which to exclude data due to CH4 emission
        h3p_expand (int): If non-zero, then expand and excluded regions to include the nearest `expand` wavelengths
        ch4_expand (int): If non-zero, then expand and excluded regions to include the nearest `expand` wavelengths

    Returns:
        np.ndarray: The H3+ mask
        np.ndarray: The CH4 mask
    """
    
    def expand_mask(mask, by):
        arr = mask.copy()
        for i in range(1, by+1):
            rolled_fwd = np.roll(mask, (i))
            rolled_fwd[:i] = 0
            rolled_bck = np.roll(mask, -i)
            rolled_bck[-i:] = 0
            arr = arr | rolled_fwd | rolled_bck
        return arr
    
    if not (2.8 < wavelengths[0] < 5.3 and 2.8 < wavelengths[-1] < 5.3):
        raise NotImplementedError("Non-LTE subtraction only works for G395H data at the moment.")
    
    h3p = pd.read_csv(constants.PATH / "data/misc/H3p_G395H_model.txt", sep=" ", names=["wavelengths", "spectrum"])
    ch4 = pd.read_csv(constants.PATH / "data/misc/CH4_G395H_model.txt", sep=" ", names=["wavelengths", "spectrum"])

    h3p_mask = []
    ch4_mask = []

    for w in wavelengths:
        i, _ = utils.find_nearest(h3p.wavelengths, w)
        h3p_mask.append((h3p.spectrum[i] / np.nanmax(h3p.spectrum)) > h3p_threshold)
        ch4_mask.append((ch4.spectrum[i] / np.nanmax(ch4.spectrum)) > ch4_threshold)

    h3p_mask = expand_mask(np.array(h3p_mask), h3p_expand)
    ch4_mask = expand_mask(np.array(ch4_mask), ch4_expand)

    return h3p_mask, ch4_mask


def subtract_h3p(spectrum, neutral_model, wavelengths, region=(3.5, 3.6), return_model=False):
    """Subtract H3+ emission form a spectrum, using a neutral model as the background
    
    Args:  
        spectrum (np.ndarray):       Measured spectrum in units of uW/cm2/sr/um
        neutral_model (np.ndarray):  NEMESIS-fitted neutral atmosphere model in units of uW/cm2/sr/um
        wavelengths (np.ndarray):    Wavelength grid
        region (float, float):       The region to use to fit the H3+ model (by default this is a triple of bright lines)
        return_model:                If True, then return the h3p object used for the fit
    
    Returns:
        np.ndarray: Spectrum with H3+ subtracted in units of uW/cm2/sr/um
        h3ppy.h3p:  If return_model, the h3p object that was used to fit
    """
    spectrum = copy.deepcopy(spectrum) * 0.01
    neutral_model = copy.deepcopy(neutral_model) * 0.01
    h3p_mask = (wavelengths > region[0]) & (wavelengths < region[1])
    h3p_wls = wavelengths[h3p_mask]

    h3p = h3ppy.h3p()
    h3p.set(R=2700, wave=h3p_wls, data=spectrum[h3p_mask] - neutral_model[h3p_mask])
    h3p.guess_density()
    fit = h3p.fit()
    vars, errs = h3p.get_results()
    h3p.set(**vars, wave=wavelengths)
    h3p_pred = h3p.model()
    
    subtracted = spectrum - h3p_pred

    if return_model:
        return subtracted * 100, h3p
    else:
        return subtracted * 100


def subtract_ch4(spectrum, neutral_model, wavelengths, region=(3.25, 3.4), include_linear=True, return_model=False):
    """Subtract CH4 non-LTE emission form a spectrum, using a neutral model as the background.
    It is highly recommended to subtract H3+ emission first using subtract_h3p()
    
    Args:  
        spectrum (np.ndarray):       Measured spectrum in units of uW/cm2/sr/um
        neutral_model (np.ndarray):  NEMESIS-fitted neutral atmosphere model in units of uW/cm2/sr/um
        wavelengths (np.ndarray):    Wavelength grid
        region (float, float):       The region to use to fit the H3+ model (by default this is a triple of bright lines)
        include_linear (bool):       Whether to include the linear offset in the fit or not (fits better, but is it physical?)
        return_model:                If True, then return the parameters and covaraince matrix for the fit
    
    Returns:
        np.ndarray: Spectrum with CH4 subtracted units uW/cm2/sr/um
        np.ndarray: Only if return_model, the optimised parameters
        np.ndarray: Only if return_model, the covariance matrix
        function:   Only if return_model, the CH4 fit function
    """

    ch4_mask = (wavelengths > 3.25) & (wavelengths < 3.4)
    ch4_wls = wavelengths[ch4_mask]
    ch4_data = spectrum - neutral_model
    ch4_data[~ch4_mask] = 0
    ch4lines = ch4.fit_non_LTE_CH4_JWST(wavelengths)

    if include_linear:
        def ch4_model(_, fun_scale, hot_scale, m, c):
            lines = ch4lines.ch4_fun*fun_scale + ch4lines.ch4_hot*hot_scale
            x = np.arange(0, len(ch4lines.ch4_fun), 1)
            linear = x*m + c
            return lines * linear
    else:
        def ch4_model(_, fun_scale, hot_scale):
            lines = ch4lines.ch4_fun*fun_scale + ch4lines.ch4_hot*hot_scale
            return lines
        
    popt, pcov = scipy.optimize.curve_fit(ch4_model, ch4_wls, ch4_data)
    ch4_pred = ch4_model(None, *popt)
    subtracted = spectrum - ch4_pred

    if return_model:
        return subtracted, popt, pcov, ch4_model
    else:
        return subtracted


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


def get_single_spectra(file_pattern, 
                       out_filename=None, 
                       max_emission=99, 
                       margins=0, 
                       pct_error=0, 
                       num_points=None, 
                       min_wl=-999, 
                       max_wl=999, 
                       remove_nonlte=True, 
                       **nonlte_kwargs):
    """Get a single spectrum from multiple observations, with options to restrict emission angle, downsample, and restrict wavelength range. 
       Saves the output to a .spx file.

       Args:
           file_pattern (str):   File pattern to match observation files.
           out_filename (str):   Output filename template for saving the processed spectrum (if None then don't save). 
                                 This can contain any of the parameters passed into the function surrounded 
                                 by curly braces. eg. "spectra_n{num_points}.spx" will become "spectra_n80.spx" 
                                 if num_points=80 is passed in.
           max_emission (float): Maximum allowable emission angle.
           margins (int):        Number of pixels around the edge of the image to remove.
           pct_error (float):    Flat percentage error to apply to the spectrum.
           num_points (int):     Number of points to downsample the spectrum to.
           min_wl (float):       Minimum wavelength in the spectrum.
           max_wl (float):       Maximum wavelength in the spectrum.
           remove_nonlte (bool): Whether to remove non-LTE emission from the spectrum.
           nonlte_kwargs:        Keyword arguments to pass to get_nonlte_emission_mask.
           
       Returns:
           np.ndarray: The wavelengths in the new .spx file
           np.ndarray: The radiances in the new spx file (returns MJy/sr, writes uW/cm2/sr/um to the .spx file)
           np.ndarray: The radiance error in the new spx file (same units as radiance)
    """
    
    params = locals().copy()
    observations = get_observations(file_pattern)
    wavelengths = observations[0].get_wavelengths_from_header()

    cubes = []
    weights = []
    other = defaultdict(list)

    for obs in observations:

        mask = np.full(obs.data.shape, np.nan)
        mask[:, obs.get_emission_angle_img() < max_emission] = 1
        mask = margin_trim(mask, margin_size=margins)
        if np.all(np.isnan(mask)):
            continue

        cubes.append(obs.data * mask)
        weights.append(np.sum(np.where(~np.isnan(mask), 1, 0)))

        other["lat"].append(np.nanmean(obs.get_lat_img() * mask))
        other["lon"].append(np.nanmean(obs.get_lon_img() * mask))
        other["phase"].append(np.nanmean(obs.get_phase_angle_img() * mask))
        other["emission"].append(np.nanmean(obs.get_emission_angle_img() * mask))
        other["azimuth"].append(np.nanmean(obs.get_azimuth_angle_img() * mask))

    w = wavelengths
    s = multiple_cube_average(cubes)
    s, w = trim_spectra(s, w, min_wl=min_wl, max_wl=max_wl)

    if remove_nonlte:
        h3p_mask ,ch4_mask = get_nonlte_emission_mask(w, **nonlte_kwargs)
        s = s[h3p_mask | ch4_mask]
        w = w[h3p_mask | ch4_mask]

    if num_points is not None:
        s, w = downsample_spectra(s, w, num_points=num_points)

    for name, values in other.items():
        other[name] = utils.nanaverage(values, weights)

    err = s * pct_error

    if out_filename is not None:
        spx.write(out_filename.format(**params), 
                spectrum=s, 
                error=err, 
                wavelengths=w, 
                fwhm=0, 
                **other)
    
    return w, s, err
