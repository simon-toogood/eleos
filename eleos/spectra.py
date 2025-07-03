import numpy as np
import pandas as pd
import glob
import h3ppy
import planetmapper
import scipy
import copy
from astropy.io import fits
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
from scipy.signal import butter, lfilter

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


def downsample_spectrum(wavelengths, spectrum, errors=None, num_points=100, special_regions=None):
    """Downsample a spectrum to a given number of spectral points.
    
    Args:
        wavelengths (np.ndarray): The corresponding wavelengths
        spectrum (np.ndarray):    The spectral data
        error (np.ndarray):       The error on the spectral data (optional)
        num_points (int):         The number of points to reduce the spectrum to
        special_regions (dict([int, int]: int)):
                                  A dict of regions of the form 
                                  {(start, end): num_points} where the spectrum should be downsampled 
                                  to a specific number of points. Use -1 for full-resolution sampling
        
    Returns:
        wavelengths (np.ndarray): The trimmed wavelengths
        spectrum (np.ndarray):    The trimmed spectrum
        errors (np.ndarray):      The trimmed errors (if provided)"""
    
    wavelengths = np.asarray(wavelengths)
    spectrum = np.asarray(spectrum)
    if errors is not None:
        errors = np.asarray(errors)

    # Keep track of indices to include
    selected_indices = set(np.linspace(0, len(wavelengths) - 1, num_points, dtype=int))

    # Handle special regions
    if special_regions:
        for (start, end), region_points in special_regions.items():
            region_mask = (wavelengths >= start) & (wavelengths <= end)
            region_indices = np.where(region_mask)[0]

            # Keep all points in this region
            if region_points == -1:
                selected_indices.update(region_indices)

            # Or downsample appropriately
            else:
                if region_points > len(region_indices):
                    region_points = len(region_indices)
                downsampled = np.linspace(0, len(region_indices) - 1, region_points, dtype=int)
                selected_indices.update(region_indices[downsampled])

    # Final selection
    final_indices = sorted(selected_indices)
    if errors is not None:
        return wavelengths[final_indices], spectrum[final_indices], errors[final_indices]
    return wavelengths[final_indices], spectrum[final_indices]


def get_nonlte_emission_mask(wavelengths, 
                             h3p_threshold=0.1, 
                             ch4_threshold=0.1, 
                             h3p_expand=3, 
                             ch4_expand=3):
    """Use a model to mask out any parts of the spectrum where H3+ and CH4 non-LTE
    emission is above a threshold relative to their highest peaks.
    
    Args:
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


def subtract_non_lte(wavelengths, spectrum, print_info=False):
    """Utility function for combining compute_lowpass_background, subtract_h3p and subtract_ch4.
    For more fine-grained control, use those function independently.
    
    Args:
        wavelengths (np.ndarray): Wavelength grid (microns)
        spectrum (np.ndarray):    Spectral data (in units of uW/cm2/sr/um)
        print_info (bool):        Whether to print the fitted parameters (eg. H3+ temp/density, CH4 scaling...)
        
    Returns:
        np.ndarray: The spectrum with H3+ and CH4 removed"""

    bkg = compute_lowpass_background(wavelengths, spectrum)
    h3p_subtracted, h3p = subtract_h3p(wavelengths, 
                                       spectrum, 
                                       bkg, 
                                       return_model=True)
    ch4_subtracted, popt, pcov, func = subtract_ch4(wavelengths, 
                                                    h3p_subtracted, 
                                                    bkg, 
                                                    return_model=True)
    ch4_subtracted[ch4_subtracted < 1e-3] = np.nan

    if print_info:
        print(f"H3+ Temperature: {h3p.vars['temperature']}K")
        print(f"H3+ Density: {h3p.vars['density']}m-2\n")
        print(f"CH4 Fundamental scaling: {popt[0]}")
        print(f"CH4 Hot scaling: {popt[1]}")
        print(f"CH4 Linear gradient: {popt[2]}")
        print(f"CH4 Linear intercept: {popt[3]}")


    return ch4_subtracted


def compute_spline_background(wavelengths, spectrum, bin_width=10, return_spline=False, **spline_kwargs):
    """Compute a polynomial background for a spectrum  by binning the spectrum and finding the minima
    in each bin, then fitting a polynomial to those points.
    
    Args:
        wavelengths (np.ndarray): The corresponding wavelengths
        spectra (np.ndarray):     The spectral data
        bin_width (float):        The width of the bins to use for the polynomial fit (in index units, not wavelength units)
        return_spline (bool):     If True, then return the polynomial coefficients instead of the evaluated background
        **spline_kwargs:          Additional keyword arguments to pass to the UnivariateSpline constructor

    Returns:
        np.ndarray: The polynomial background evaluated at the given wavelengths
        scipy.interpolate.UnivariateSpline: Only if return_spline, the spline object
    """
    
    def bottom_nth(wavelengths, spectrum, nm=10, nx=1):
        # Find the 5th percentile value
        lower = np.nanpercentile(spectrum, 0.01)
        upper = 99#np.nanpercentile(spectrum, 100 - nm)
        mask = (spectrum > lower) & (spectrum < upper)
        low_vals = spectrum[mask]
        low_waves = wavelengths[mask]

        # Use mean of low values and corresponding wavelengths
        avg_wave = np.nanmean(low_waves)
        avg_log_val = np.nanmean(np.log(low_vals))

        return avg_wave, avg_log_val

    num_points = len(spectrum)
    bin_edges = np.arange(0, num_points, bin_width, dtype=int)

    bin_centers = []
    bin_log_minima = []

    for i in range(len(bin_edges) - 1):
        start = bin_edges[i]
        end = bin_edges[i + 1]
        bin_slice = spectrum[start:end]
        bin_waves = wavelengths[start:end]

        if len(bin_slice) == 0:
            continue

        avg_wave, avg_log_val = bottom_nth(bin_waves, bin_slice)

        bin_centers.append(avg_wave)
        bin_log_minima.append(avg_log_val)

    spl = UnivariateSpline(bin_centers, bin_log_minima, **spline_kwargs)
    bkg = np.exp(spl(wavelengths))

    if return_spline:
        return bkg, spl
    else:
        return bkg


def compute_lowpass_background(wavelengths, spectrum):
    
    def butter_lowpass_filter(data, cutoff, fs, order=5):
        b, a = butter(order, cutoff, fs=fs, btype='low', analog=False)
        y = lfilter(b, a, data)
        return y

    h3pmask, ch4mask = get_nonlte_emission_mask(wavelengths)
    mask = h3pmask | ch4mask | np.isnan(spectrum)

    spectrum_interp = np.interp(wavelengths[mask], wavelengths[~mask], spectrum[~mask])
    spec = spectrum.copy()
    spec[mask] = spectrum_interp
    n = 17
    spec = np.concatenate((spec[n:], np.full(n, spec[-1])))

    fs = 1 / (wavelengths[1] - wavelengths[0])
    cutoff = 50
    bkg = np.exp(butter_lowpass_filter(np.log(spec), cutoff, fs))

    return bkg


def interpolate_nans(wavelengths, spectrum):
    mask = ~np.isnan(spectrum)
    return np.interp(wavelengths, wavelengths[mask], spectrum[mask])


def subtract_h3p(wavelengths, spectrum, neutral_model, region=(3.5, 3.6), return_model=False, **h3p_kwargs):
    """Subtract H3+ emission form a spectrum, using a neutral model as the background
    
    Args:  
        wavelengths (np.ndarray):    Wavelength grid
        spectrum (np.ndarray):       Measured spectrum in units of uW/cm2/sr/um
        neutral_model (np.ndarray):  NEMESIS-fitted neutral atmosphere model in units of uW/cm2/sr/um
        region (float, float):       The region to use to fit the H3+ model (by default this is a triple of bright lines)
        return_model (bool):         If True, then return the h3p object used for the fit
    
    Returns:
        np.ndarray: Spectrum with H3+ subtracted in units of uW/cm2/sr/um
        h3ppy.h3p:  If return_model, the h3p object that was used to fit
    """
    spectrum = copy.deepcopy(spectrum) * 0.01
    spectrum = interpolate_nans(wavelengths, spectrum)
    neutral_model = copy.deepcopy(neutral_model) * 0.01
    h3p_mask = (wavelengths > region[0]) & (wavelengths < region[1])
    h3p_wls = wavelengths[h3p_mask]

    h3p = h3ppy.h3p()
    h3p.set(R=2700, **h3p_kwargs, wave=h3p_wls, data=spectrum[h3p_mask])# - neutral_model[h3p_mask])
    fit = h3p.fit()
    h3p.set(wave=wavelengths)
    h3p_pred = h3p.model(nbackground=0)
    
    subtracted = spectrum - h3p_pred

    if return_model:
        return subtracted * 100, h3p
    else:
        return subtracted * 100


def subtract_ch4(wavelengths, spectrum, neutral_model, region=(3.25, 3.4), include_linear=True, return_model=False):
    """Subtract CH4 non-LTE emission form a spectrum, using a neutral model as the background.
    It is highly recommended to subtract H3+ emission first using subtract_h3p()
    
    Args:  
        wavelengths (np.ndarray):    Wavelength grid
        spectrum (np.ndarray):       Measured spectrum in units of uW/cm2/sr/um
        neutral_model (np.ndarray):  NEMESIS-fitted neutral atmosphere model in units of uW/cm2/sr/um
        region (float, float):       The region to use to fit the H3+ model (by default this is a triple of bright lines)
        include_linear (bool):       Whether to include the linear offset in the fit or not (fits better, but is it physical?)
        return_model (bool):         If True, then return the parameters and covaraince matrix for the fit
    
    Returns:
        np.ndarray: Spectrum with CH4 subtracted in units uW/cm2/sr/um
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


def zonal_average_tiles(filepaths, lat_width, error_scale=1, rmsd_threshold=10, filters=None):
    """Get the zonal average from a set of JWST navigated cubes.
    
    Args:
        filepaths (List(str)): List of filepaths to tiles to use in the averaging
        lat_width (float):     The width of the latitude bins. The given latitude is in the centre of the bin.
        error_scale (float):   Multiply the average error by this amount (unused so far)
        rmsd_threshold (int):  Reject any spaxels that are in the top x percentile for root-mean-square deviation from the median
        filters (dict):        Apply any filters to remove spaxels before averaging. Format is key=planetmapper backplane name (str)
                            and value=(min, max). For example, 'EMISSION':(0, 70) will reject any spaxels with emission angle greater than 70
                            
    Returns:
        np.ndarray:       Wavelength grid, temporary - will be added to df later
        pandas.DataFrame: DataFrame containing the zonal means and information for .spx files.
                          Columns are spectrum, error (np.ndarray's), phase, emission, azimuth, lon (floats), num_spaxels (int).
                          Use df.index to get latitudes"""
    
    if filters is None:
        filters = dict()
  
    def mask_and_reshape(arr, mask, shape=(-1,)):
        return (arr * mask).reshape(*shape) 

    def compute_rmsd(group):
        # Stack the spectra for this group
        stacked = np.stack(group['spectrum'].values)
        
        # Compute the median spectrum
        median = np.nanmedian(stacked, axis=0)
        
        # Compute RMSD for each row
        # This results in an array of shape (n_rows,)
        rmsd = np.sqrt(np.nanmean((stacked - median)**2, axis=1))
        
        # Return the group with an additional column
        group = group.copy()
        group['rmsd'] = rmsd
        return group

    def nanmean_column(group, column):
        return np.nanmean(np.stack(group[column].values), axis=0)

    # Dicts to hold individual spaxels, final zonal spectra, various angles
    spaxels = pd.DataFrame(columns=["spectrum", "error", "lat", "phase", "emission", "azimuth", "lon", "nans"])
    a = 0
    # Iterate over each filepath to bin the spaxels by latitude
    for tile in filepaths:
        print("Processing ", tile)

        # Get planetmapper Observation object
        obs = planetmapper.Observation(tile)
        obs = add_error_cubes(obs)

        # Create the mask from filters
        mask = np.ones_like(obs.data[0]).astype(bool)
        for bp_name, bounds in filters.items():
            backplane = obs.get_backplane_img(bp_name)
            single_mask = (bounds[0] < backplane) & (backplane  < bounds[1])
            mask &= single_mask

        lonimg = obs.get_lon_img()
        emiimg = obs.get_emission_angle_img()
        aziimg = obs.get_azimuth_angle_img()
        phaimg = obs.get_phase_angle_img()

        # Bin each spaxel into latitude bins
        lat_bp = obs.get_lat_img()
        for l in np.arange(-90, 90, lat_width):
            lat_mask = (l - lat_width/2 < lat_bp) & (l + lat_width/2 > lat_bp)
            final_mask = lat_mask * mask
            if not (final_mask).any():
                continue
            #print(l)

            data = mask_and_reshape(obs.data,  final_mask, shape=(obs.data.shape[0], -1))
            err  = mask_and_reshape(obs.error, final_mask, shape=(obs.data.shape[0], -1))
            lons = mask_and_reshape(lonimg, final_mask)
            emis = mask_and_reshape(emiimg, final_mask)
            azim = mask_and_reshape(aziimg, final_mask)
            phas = mask_and_reshape(phaimg, final_mask)

            # Remove any all-zero or all-NaN spaxels
            for i in range(data.shape[1]):
                if not ((data[:,i] < 1e-10) | (np.isnan(data[:,i]))).all():
                    spaxels.loc[a] = [data[:,i].copy(), err[:, i].copy(), l, phas[i], emis[i], azim[i], lons[i], np.sum(np.isnan(data[:,i]))]
                    a += 1
            

    # Remove any spaxels with a really large number of NaNs (basically the very edges of the cube)
    spaxels.drop(spaxels[spaxels.nans > 200].index, inplace=True)
    # Remove any spaxels that are a certain threshold away from the root-mean-square deviation from the median
    spaxels = spaxels.groupby('lat', group_keys=False).apply(compute_rmsd)
    spaxels.drop(spaxels[(spaxels.rmsd > np.percentile(spaxels.rmsd, 100-rmsd_threshold))].index, inplace=True)
    # Calculate final zonal spectra
    z = spaxels.groupby('lat')
    zonal = z.mean()
    zonal['spectrum'] = z.apply(nanmean_column, "spectrum")
    zonal['error'] = z.apply(nanmean_column, "error")
    zonal["num_spaxels"] = z.apply(len)
    del zonal["nans"], zonal["rmsd"]

    return obs.get_wavelengths_from_header(), zonal


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
    if isinstance(observations, planetmapper.Observation):
        observations = [observations,]
        return_single = True
    else:
        return_single = False

    for obs in observations:
        fp = obs.path 
        hdul = fits.open(fp)
        obs.error = hdul['ERR'].data

    if return_single:
        return observations[0]
    else:
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
