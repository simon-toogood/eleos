import numpy as np
import pandas as pd
import glob
import h3ppy
import planetmapper
import copy
from astropy.io import fits
from collections import defaultdict
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve

from . import utils
from . import spx
from . import constants


def margin_trim(cube, margin_size=3):
    """
    Set the outer rings of spaxels in a cube to NaN.

    Args:
        cube (np.ndarray): The spectral cube to trim with shape (wavelengths, x, y).
        margin_size (int): The size of the margin to trim from each edge (eg. 2 will set every spaxel within 2 spaxels of the edge to NaN).

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


def resample_to_new(wavelength_source, spectrum_source, wavelength_target, kind="linear", fill_value="extrapolate"):
    """Resample a spectrum onto a new wavelength grid via interpolation.

    Args:
        wavelength_source (array-like): Original wavelength values of the spectrum (must be sorted).
        spectrum_source (array-like): Spectrum values corresponding to `wavelength_source`.
        wavelength_target (array-like): Target wavelength grid to resample onto.
        kind (str, optional): Interpolation type. Options include 
            "linear", "nearest", "cubic", etc. Defaults to "linear".
        fill_value (str or float, optional): Value to use outside the range of 
            `wavelength_source`. If "extrapolate", allows extrapolation. Defaults to "extrapolate".

    Returns:
        np.ndarray: Spectrum values resampled onto the `wavelength_target` grid.


    """
    wavelength_source = np.asarray(wavelength_source)
    spectrum_source = np.asarray(spectrum_source)
    wavelength_target = np.asarray(wavelength_target)

    if len(wavelength_source) != len(spectrum_source):
        raise ValueError("wavelength_source and spectrum_source must have the same length.")

    interp_func = interp1d(wavelength_source, spectrum_source, kind=kind,
                           fill_value=fill_value, bounds_error=False)
    return interp_func(wavelength_target)


def wavelength_select(wavelengths, spectrum, errors=None, min_wl=None, max_wl=None, epsilon=0):
    """Select a range of wavelengths from a spectrum.
    
    Args:
        wavelengths (np.ndarray): The wavelength grid
        spectrum (np.ndarray):    The spectral data
        errors (np.ndarray):      The error on the spectral data (optional)
        min_wl (float):           The minimum wavelength to select
        max_wl (float):           The maximum wavelength to select
        epsilon (float):          A small fudge factor to the end of the mask to fix FPE when grouping spectra
        
    Returns:
        wavelengths (np.ndarray): The trimmed wavelengths
        spectrum (np.ndarray):    The trimmed spectrum
        errors (np.ndarray):      The trimmed errors (if provided)"""
    
    
    mask = np.ones_like(wavelengths, dtype=bool)
    if min_wl is not None:
        mask &= wavelengths >= min_wl
    if max_wl is not None:
        mask &= wavelengths <= max_wl + epsilon

    if errors is not None:
        return wavelengths[mask], spectrum[mask], errors[mask]
    return wavelengths[mask], spectrum[mask]


def interpolate_nans(wavelengths, spectrum):
    mask = ~np.isnan(spectrum)
    return np.interp(wavelengths, wavelengths[mask], spectrum[mask])


def subtract_non_lte(wavelengths, spectrum):
    """Utility function for combining subtract_h3p and subtract_ch4.
    For more fine-grained control, use those function independently.
    
    Args:
        wavelengths (np.ndarray): Wavelength grid (microns)
        spectrum (np.ndarray):    Spectral data (in units of W/cm2/sr/um)
        
    Returns:
        np.ndarray: The wavelength grid (unchanged)
        np.ndarray: The spectrum with H3+ and CH4 removed"""

    h3p_subtracted = subtract_h3p(wavelengths, spectrum)
    ch4_subtracted = subtract_ch4(wavelengths, h3p_subtracted)

    return wavelengths, ch4_subtracted


def subtract_h3p(wavelengths, spectrum, latitude=-60, region=(3.525, 3.55), return_model=False, **h3p_kwargs):
    """Subtract H3+ emission from a spectrum

    Args:  
        wavelengths (np.ndarray): Wavelength grid
        spectrum (np.ndarray):    Measured spectrum in units of W/cm2/sr/um
        latitude (float):         Optionally specify a latitude to make a more informed guess at initial values
        region (float, float):    The region to use to fit the H3+ model (by default this is a triple of bright lines)
        return_model (bool):      If True, then return the h3p object used for the fit
    
    Returns:
        np.ndarray: Spectrum with H3+ subtracted in units of W/cm2/sr/um
        h3ppy.h3p:  If return_model, the h3p object that was used to fit
    """
    def remap(x, in_min, in_max, out_min, out_max):
        return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min

    temp = remap(np.abs(latitude), 40, 90, 500, 1000)
    density = remap(np.abs(latitude), 40, 90, 1e18, 1e19)

    h3p_kwargs.setdefault("temperature", temp)
    h3p_kwargs.setdefault("density", density)
    h3p_kwargs.setdefault("R", 1700)

    spectrum = copy.deepcopy(spectrum) * 10000 # Convert to W/m2/sr/um
    nanmask = np.isnan(spectrum)
    spectrum = interpolate_nans(wavelengths, spectrum)
    w, s = wavelength_select(wavelengths, spectrum, min_wl=region[0], max_wl=region[1])

    h3p = h3ppy.h3p()
    h3p.set(**h3p_kwargs, wave=w, data=s)
    # h3p.guess_density()
    fit = h3p.fit(verbose=False)
    h3p.set(wave=wavelengths)
    plt.plot(wavelengths, h3p.model() / 10000)
    h3p_pred = h3p.model(background_0=0)
    plt.plot(wavelengths, h3p.model(background_0=0) / 10000)

    subtracted = (spectrum - h3p_pred) / 10000  # Convert back to W/cm2/sr/um
    subtracted[nanmask] = np.nan  # Restore NaNs where they were in the original spectrum
    subtracted[subtracted < 0] = np.nan  # Remove any negative values
    if return_model:
        return subtracted, h3p
    else:
        return subtracted


def subtract_ch4(wavelengths, spectrum, nlines=25, linewidth=0.008, A0=1e-6, lineshape="voigt", return_model=False):
    """Subtract the brightest N individual CH4 lines from a spectrum.
    
    Args:
        wavelengths (np.ndarray): The wavelength grid
        spectrum (np.ndarray):    The spectral data in units of W/cm2/sr/um
        nlines (int):             The number of peaks to fit
        linewidth (float):        The intial guess for the line FWHM in microns
        A0 (float):               The initial guess for the amplitude of the lines
        lineshape (str):          The shape of the lines to use, either 'gaussian', 'lorentzian', or 'voigt'
        return_model (bool):      If True, then return the model spectrum as well as the final spectrum
        
    Returns:
        np.ndarray: The spectrum with CH4 subtracted in units of W/cm2/sr/um
    """

    raise NotImplementedError("CH4 subtraction is dodgy at best. Catch this error if you are sure you want to try")

    def multifunc(x, func, p0s):
        total = np.zeros_like(x)
        for p in p0s:
            total += func(x, *p)
        return total 
    
    def fit(wavelengths, spectrum, line, linewidth, c="magenta"):
        if not np.nanmin(wavelengths) < line < np.nanmax(wavelengths):
            print("fail")
            return
        
        w, s = wavelength_select(wavelengths, 
                                 spectrum, 
                                 min_wl=line - linewidth/1.5, 
                                 max_wl=line + linewidth/1.5)
        offset_0 = np.nanmedian(np.nanpercentile(s, 25))
        plt.hlines(y=offset_0/1e8, xmin=line-linewidth/2, xmax=line+linewidth/2, color=c)
        
        try:
            popt, pcov = curve_fit(eval("utils."+lineshape), 
                                   w, s, 
                                   p0=     [A0,      line,                linewidth/2,    offset_0],
                                   bounds=([0,       line - linewidth/2,  linewidth*0.2,  offset_0*0.9], 
                                           [np.inf,  line + linewidth/2,  linewidth*2,    offset_0*1.1]),)
            plt.plot(w, eval("utils."+lineshape)(w, *popt)/1e8, color=c)
            return popt

        except (RuntimeError, ValueError) as e:
            print(f"[eleos] Failed to fit CH4 line at {line:.3f} um: {e}")


    nanmask = np.isnan(spectrum)
    spectrum = interpolate_nans(wavelengths, spectrum) # linearly interpolate nans
    spectrum *= 1e8 # numerical stability
    lines = []
    i = 0
    centre_l = 3.3
    centre_r = 3.35

    # Extract the N most intense lines from the CH4 spectrum that have a minimum spacing (prevent picking out double lines)
    for _, line in constants.CH4_LINES.iterrows():
        wavelength = line["wavelength"]
        if all(abs(wavelength - existing) > linewidth for existing in lines) and not centre_l < wavelength < centre_r:
            lines.append(wavelength)
            i += 1
        if i >= nlines:
            break

    # Fit each of those lines
    fits = []
    for line in lines:
        popt = fit(wavelengths, spectrum, line, linewidth)
        if popt is None:
            continue
        fits.append(popt)

    # Now fit the extra lines that are too weak to be captured in the line list
    # probably not a good idea!
    # if extra_lines:
    #     extralines = [3.478, 3.491, 3.504, 3.517, 3.166, 3.157, 3.149, 3.142, 3.133, 3.124]
    #     for line in extralines:
    #         popt = fit(wavelengths, spectrum, line, linewidth, c="green")
    #         if popt is None:
    #             continue
    #         fits.append(popt)

    # Set the offset to 0 for subtraction
    for i in range(len(fits)):
        fits[i][3] = 0

    # Add all the lines together to make the final model
    model = multifunc(wavelengths, eval("utils."+lineshape), fits)
    subtracted = spectrum - model

    # Restore original NaNs to the spectrum
    subtracted[nanmask] = np.nan 

    # Remove the massive central peak completely
    subtracted[(wavelengths > centre_l) & (wavelengths < centre_r)] = np.nan
    model[(wavelengths > centre_l) & (wavelengths < centre_r)] = np.nan
    subtracted /= 1e8 
    model /= 1e8

    if return_model:
        return subtracted, model
    else:
        return subtracted


def zonal_average_tiles(filepaths, lat_width, error_scale=1, rmsd_threshold=10, filters=None):
    """Get the zonal averages from a set of JWST navigated cubes.
    
    Args:
        filepaths (List(str)): List of filepaths to tiles to use in the averaging
        lat_width (float):     The width of the latitude bins. The given latitude is in the centre of the bin.
        error_scale (float):   Multiply the average error by this amount (unused so far)
        rmsd_threshold (int):  Reject any spaxels that are in the top x percentile for root-mean-square deviation from the median
        filters (dict):        Apply any filters to remove spaxels before averaging. Format is key=planetmapper backplane name (str)
                               and value=(min, max). For example, filters={'EMISSION':(0, 70)} will reject any spaxels with emission angle greater than 70
                            
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
        print("[eleos] Processing ", tile)

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
                       max_wl=999):
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
           
       Returns:
           np.ndarray: The wavelengths in the new .spx file
           np.ndarray: The radiances in the new spx file (returns MJy/sr, writes W/cm2/sr/um to the .spx file)
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
    s, w = wavelength_select(w, s, min_wl=min_wl, max_wl=max_wl)

    if num_points is not None:
        s, w = downsample_spectrum(s, w, num_points=num_points)

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


def combine_multiple_spectra(*spectra_units):
    """
    Combine multiple spectra with errors. Preserves original wavelength resolution in
    non-overlapping regions and uses weighted average (inverse variance) in overlaps.
    
    Args:
        *spectra_units: Tuples of (wavelength, spectrum, error), where each element is a 1D array.
    
    Returns:
        wavelengths (np.array): Stitched wavelength grid
        spectrum (np.array): Combined spectral values
        error (np.array): Combined errors
    """
    segments = []
    spectra = list(spectra_units)
    spectra.sort(key=lambda t: t[0][0])  # Sort by wavelength start

    for i, (wl_i, f_i, e_i) in enumerate(spectra):
        wmin_i, wmax_i = wl_i.min(), wl_i.max()
        overlapping_segments = []

        for j, (wl_j, f_j, e_j) in enumerate(spectra):
            if i == j:
                continue
            wmin_j, wmax_j = wl_j.min(), wl_j.max()

            if wmax_j > wmin_i and wmin_j < wmax_i:
                w_overlap_min = max(wmin_i, wmin_j)
                w_overlap_max = min(wmax_i, wmax_j)

                # Choose higher-res grid in overlap
                res_i = np.median(np.diff(wl_i[(wl_i >= w_overlap_min) & (wl_i <= w_overlap_max)]))
                res_j = np.median(np.diff(wl_j[(wl_j >= w_overlap_min) & (wl_j <= w_overlap_max)]))
                grid = wl_i if res_i <= res_j else wl_j
                wl_common = grid[(grid >= w_overlap_min) & (grid <= w_overlap_max)]

                # Interpolate flux and error
                fi = interp1d(wl_i, f_i, bounds_error=False, fill_value=np.nan)(wl_common)
                fj = interp1d(wl_j, f_j, bounds_error=False, fill_value=np.nan)(wl_common)

                ei = interp1d(wl_i, e_i, bounds_error=False, fill_value=np.nan)(wl_common)
                ej = interp1d(wl_j, e_j, bounds_error=False, fill_value=np.nan)(wl_common)

                # Inverse-variance weighted average
                with np.errstate(divide='ignore', invalid='ignore'):
                    wi = 1 / (ei**2)
                    wj = 1 / (ej**2)
                    weights_sum = wi + wj

                    flux_comb = np.nansum([wi * fi, wj * fj], axis=0) / weights_sum
                    error_comb = np.sqrt(1 / weights_sum)

                overlapping_segments.append((wl_common, flux_comb, error_comb))

        # Remove overlap from i-th spectrum
        mask = np.full(wl_i.shape, True)
        for wl_o, _, _ in overlapping_segments:
            mask &= ~((wl_i >= wl_o.min()) & (wl_i <= wl_o.max()))

        # Add unique portion
        if np.any(mask):
            segments.append((wl_i[mask], f_i[mask], e_i[mask]))

        # Add overlapping portions
        segments.extend(overlapping_segments)

    # Stitch all segments together
    all_wl = np.concatenate([seg[0] for seg in segments])
    all_flux = np.concatenate([seg[1] for seg in segments])
    all_err = np.concatenate([seg[2] for seg in segments])

    idx = np.argsort(all_wl)
    return all_wl[idx], all_flux[idx], all_err[idx]


def convolve_gaussian(wavelength, data, R, lambda0):
    """Convolve a spectrum with a Gaussian, assuming constant R. In reality resolution varies as a function of wavelength and
    depends on the filters used (for JWST/NIRSpec see https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-dispersers-and-filters)
    
    Args:
        wavelength (np.ndarray): Wavelength array (same length as data).

        
                        data (np.ndarray):       The data to be convolved.
        R (float):               Resolving power of the instrument
        lambda0 (float):         The wavelength at which R is defined. For G395H this is ~4um. For G235H this is ~2.4um.
    
    Returns:
        np.ndarray: The convolved data.
    """
    wavelength = np.asarray(wavelength)
        
    # FWHM and sigma in wavelength units
    fwhm = lambda0 / R
    sigma = fwhm / 2.35482

    # Build Gaussian kernel in wavelength space
    # Kernel extends to 4 sigma for good truncation
    dw = np.median(np.diff(wavelength))  # wavelength step
    half_width = int(np.ceil(4 * sigma / dw))
    x = np.arange(-half_width, half_width + 1) * dw

    kernel = np.exp(-(x**2) / (2 * sigma**2))
    kernel /= kernel.sum()  # normalize so area = 1

    # Convolve using FFT for efficiency
    convolved = fftconvolve(data, kernel, mode='same')

    return convolved


def quickview(wavelength, spectrum, errors=None, log=True, block=True):
    """Quickly view a spectrum using matplotlib.
    
    Args:
        wavelength (np.ndarray): The wavelength grid
        spectrum (np.ndarray):   The spectral data
        errors (np.ndarray):     The error on the spectral data (optional)
        log (bool):              Whether to use a logarithmic scale for the y-axis
        block (bool):            Whether to block the execution until the plot is closed (default True
    """
    
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 5))
    plt.plot(wavelength, spectrum, label='Spectrum')
    
    if errors is not None:
        plt.fill_between(wavelength, spectrum - errors, spectrum + errors, alpha=0.3, label='Error')
    
    if log:
        plt.yscale('log')

    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Radiance')
    plt.legend()
    plt.grid()
    plt.show(block=block)