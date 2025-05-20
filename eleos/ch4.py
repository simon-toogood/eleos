import numpy as np
from astroquery.hitran import Hitran
from astropy import units as u
import scipy.constants as const
from scipy.optimize import curve_fit
from astropy.table import Table
import h3ppy


def get_hitran_table(wrange = [3.1, 3.9], molecule = 6, v_values = True) : 
    
    # Ask nicely for some data
    tbl = Hitran.query_lines(molecule_number = molecule, min_frequency=10000.0 / wrange[1] / u.cm, max_frequency=10000.0 / wrange[0] / u.cm)

    if (v_values == True) : 
        # Get the v values from the quanta fields
        # https://hitran.org/media/refs/HITRAN_QN_formats.pdf
        tbl['v1upper'] =  [int(quanta[3:5]) for quanta in tbl['global_upper_quanta']]
        tbl['v2upper'] =  [int(quanta[5:7]) for quanta in tbl['global_upper_quanta']]
        tbl['v3upper'] =  [int(quanta[7:9]) for quanta in tbl['global_upper_quanta']]
        tbl['v4upper'] =  [int(quanta[9:11]) for quanta in tbl['global_upper_quanta']]
        tbl['nupper']  =  [quanta[11:13] for quanta in tbl['global_upper_quanta']]
        tbl['Cupper']  =  [quanta[13:15] for quanta in tbl['global_upper_quanta']]
        tbl['v1lower'] =  [int(quanta[3:5]) for quanta in tbl['global_lower_quanta']]
        tbl['v2lower'] =  [int(quanta[5:7]) for quanta in tbl['global_lower_quanta']]
        tbl['v3lower'] =  [int(quanta[7:9]) for quanta in tbl['global_lower_quanta']]
        tbl['v4lower'] =  [int(quanta[9:11]) for quanta in tbl['global_lower_quanta']]
        tbl['nlower']  =  [quanta[11:13] for quanta in tbl['global_lower_quanta']]
        tbl['Clower']  =  [quanta[13:15] for quanta in tbl['global_lower_quanta']]

    return tbl


class hitran_spectrum: 
    def __init__(self, hitran_data, R = 2700) : 
        '''
        Calculate the thermal spectrum from HITRAN data. 

        The equations come from: 
        https://hitran.org/docs/definitions-and-units/
        https://hitran.org/static/hapi/hapi.py

        '''
        self.hitran_data = hitran_data
        
        # Some funky CGS constants
        self.k  = 1.380648813E-16 # erg/K, CGS
        self.c  = 2.99792458e10 # cm/s, CGS
        self.h  = 6.626196e-27 # erg*s, CGS
        self.N  = 6.02214076e23 # Avrogrados constant
        self.c2 = 1.4387769 # cm K

        # The instrument resolution in wavenumbers
        self.R_wavenumber = R / np.mean(self.hitran_data['nu'])
 
    def line_intensity_sw(self, T) : 
        p1 = self.hitran_data['a'] / (8.0 * np.pi * self.c * self.hitran_data['nu']**2)
        p2 = self.hitran_data['gp'] * np.exp(-self.c2 * self.hitran_data['elower'] / T) * (1.0 - np.exp(-self.c2 * self.hitran_data['nu'] / T))
        return p1 * p2 / self.partition_function_Q(T)

    def partition_function_Q(self, T) : 
        f = [-3.44643710e-26,  5.23961509e-23,  3.43918846e-19,  9.30259324e-16, 1.41549164e-12, 4.20644404e-10, -1.13426741e-06,  6.99297791e-03, -1.02182222e+00,  1.90545502e+02]
        qs = np.poly1d(f)
        return qs(T)

    def gamma_d(self, T, molar_mass = 16.04) : 
        m = molar_mass * self.N  #* 1000 # GSM units?
        const = np.sqrt((2.0 * self.N * self.k * T * np.log(2)) / m) 
        return self.hitran_data['nu'] / self.c * const

    def doppler_broadening(self, T, nu, offset = 0) : 
        gd = self.R_wavenumber # self.gamma_d(T) * 20e13
        exponent = ((nu - (self.hitran_data['nu'] + offset))**2 * np.log(2)) / gd**2
        return np.sqrt(np.log(2)/(np.pi * gd**2)) * np.exp(-exponent)

    def absorption_coefficient(self, T, w) : 
        return self.line_intensity_sw(T) * self.doppler_broadening(T, w)

    def spectrum(self, T, wave, path_length_cm = 1e2, offset = 0) : 
        spec = np.zeros(wave.shape[0])
        intensity = self.line_intensity_sw(T)
        for i, w in np.ndenumerate(wave) : 
            #spec[i] = np.sum(self.doppler_broadening(T, 1e4/w))
        #    Alw = 1 - np.exp(-absorption_coefficient(T, 1e4/w, tbl) * l) 
        #    LBBTw = 2*const.h*const.c**2*(1e4/w)**3 / ((np.exp(const.h * const.c * (1e4/w))/(const.k * T)) - 1) * 1.0E-7
        #    Xsect = Alw * LBBTw # W/sr/cm**2/cm**-1
            wavenumber = 1e4 / ( w + offset)
            aa      = intensity * self.doppler_broadening(T, wavenumber)
            acoff   = np.sum(aa)
            Alw     = 1-np.exp(-acoff * path_length_cm)
            LBBTw   = 2 * self.h * self.c**2 * wavenumber**3 / (np.exp(self.h * self.c * wavenumber/(self.k * T)) - 1) #* 1.0E-7
            Xsect   = Alw * LBBTw # W/sr/cm**2/cm**-1
            spec[i] =  np.sum(aa)

        # Phew, finally some SI units. Convert to W/cm2/sr/µm
        h=6.626e-34
        c=2.9979e8
        k=1.3806e-23
        spec = spec*1e4/100.
        v=100.*1e4/wave # wavenumber in m
        c1=2*h*c*c
        c2=h*c/k
        a=c1*v*v*v/spec
        TB=c2*v/(np.log(a+1))
        # Now convert back to W/cm2/sr/µm
        l=wave*1e-6
        a=2*h*c*c/(l**5)
        b=np.exp(h*c/(l*k*TB)) - 1
        return (a/b)/1e4/1e6


class methane_fitting: 
    def __init__(self, tbl, nbackground = 2, R = 1800) :

        # Select the v3 -> ground lines as per Sanches-Lopez et al., (2022, A&A)
        tbl2 = tbl[tbl['v3upper'] == 1]
        tbl2 = tbl2[tbl2['v1upper'] == 0]
        tbl2 = tbl2[tbl2['v2upper'] == 0]
        tbl2 = tbl2[tbl2['v1lower'] == 0]

        # Select the v3 + v4 -> v3 lines
        tbl3 = tbl[tbl['v3upper'] == 1]
        tbl3 = tbl3[tbl3['v4upper'] == 1]
        #tbl3 = tbl3[tbl3['v1upper'] == 0]
        #tbl3 = tbl3[tbl3['v2upper'] == 0]
        tbl3 = tbl3[tbl3['v4lower'] == 1]
        #tbl3 = tbl3[tbl3['v1lower'] == 0]
        #tbl3 = tbl3[tbl3['v2lower'] == 0]
        #tbl3 = tbl3[tbl3['v3lower'] == 0]

        self.fundametal_lines = tbl2
        self.hotband_lines = tbl3
        self.R = R
        self.wavenumber_offset = 0
        self.wavelength_offset = 0

        self.background  = np.array((nbackground))
        self.nbackground = nbackground

        self.slope = 0
        self.background1 = 0
        self.background2 = 0
        self.background3 = 0

        self.temperature_fundamental = 300 # K 
        self.temperature_hotband     = 300 # K 

        self.ch4_fundamental = hitran_spectrum(self.fundametal_lines, R = self.R)
        self.ch4_hotband     = hitran_spectrum(self.hotband_lines, R = self.R)

    def methane_spectrum_fn(self, eks, level_fundamental, level_hotband, background1, background2, background3, slope) : 
        '''
            The CH4 spectral function.
        '''
        self.level_fundamental  = level_fundamental
        self.level_hotband      = level_hotband
        self.background1        = background1
        self.background2        = background2
        self.background3        = background3
        self.slope              = slope

        # Calculate the CH4 fundamental, hotband and background
        return self.methane_fundamental_fn() +  self.methane_hotband_fn() + self.background_fn()
 
    def model(self, level_fundamental = 1, level_hotband = 1, background1 = 0, background2 = 0, background3 = 0) :


        #self.wave = wave
        #self.xpixels = np.arange(0, wave.shape[0])
        
        self.level_fundamental  = level_fundamental
        self.level_hotband      = level_hotband
        self.background1        = background1
        self.background2        = background2
        self.background3        = background3

        return self.methane_fundamental_fn() +  self.methane_hotband_fn() + self.background_fn()

    def methane_spectrum_Roffset_fn(self, eks, R, offset) : 

        # Set the line width-ish
        self.ch4_fundamental.R_wavenumber = R / np.mean(self.fundametal_lines['nu'])
        self.ch4_hotband.R_wavenumber = R / np.mean(self.fundametal_lines['nu'])

        # Set the wavenumber offset 
#        self.wavenumber_offset = offset 
        self.wavelength_offset = offset
        self.wave = self.wave_original + self.wavelength_offset

        # Calculate the CH4 fundamental, hotband and background
        return self.methane_fundamental_fn() +  self.methane_hotband_fn() + self.background_fn()

    def methane_fundamental_fn(self) : 
        # should be called spectrum_fundamental
        spectrum =  self.ch4_fundamental.spectrum(self.temperature_fundamental, self.wave) #, offset = self.wavenumber_offset)


        # Hack to simulate non-LTE population of the fundamental band
        xpixels = np.copy(self.xpixels) # np.arange(wl.shape[0])
        xpixels = xpixels / np.max(xpixels)
        slope    = np.sin(np.deg2rad(self.slope)) * 10.0 * xpixels + 1.0

        return self.level_fundamental * spectrum * slope

    def methane_hotband_fn(self) : 
        spectrum =  self.ch4_hotband.spectrum(self.temperature_hotband, self.wave) #, offset = self.wavenumber_offset)
        return self.level_hotband * spectrum 

    def methane_only(self) : 
        return self.methane_fundamental_fn() +  self.methane_hotband_fn()

    def background_fn(self) : 
        bkg = self.background1 
        if (self.nbackground > 1) : 
            bkg += self.background2 * self.xpixels
        if (self.background > 2) : 
            bkg += self.background3 * self.xpixels**2
        return bkg
#        return self.background1 + self.background2 * self.xpixels + self.background3 * self.xpixels**2

    def set_data(self, wave, spec, sigma) : 
        '''
        wave : microns
            The wavelength scale (x)
        spec : W/cm2/sr/micron
            The spectrum to be fitted (y)
        sigma: W/cm2/sr/micron
            The uncertainty on the spectrum (dy)
        '''
        self.wave = wave
        self.xpixels = np.arange(0, wave.shape[0])
        self.spec = spec 

    def fit(self, wave, spec, sigma) : 

        self.xpixels = np.arange(0, wave.shape[0])
        self.wave_original = wave
        self.wave          = wave
        self.spec          = spec 

        # Get the level = 1 spectrum to use as an inital guess
        l1 = self.methane_spectrum_fn(self.wave, 1, 1, 0, 0, 0, 0)
        initial_guess = [np.max(spec)/np.max(l1), np.max(spec)/np.max(l1), 0, 0, 0, 0]
        popt, pcov = curve_fit(self.methane_spectrum_fn, self.wave, self.spec, sigma = sigma, p0 = initial_guess)
        perr = np.sqrt(np.diag(pcov))
        #print(popt, perr)
        self.vars = popt
        self.errs = perr
        # print('yas')

        #return self.methane_spectrum_fn(self.xpixels, *popt)

        if False : 
            print(popt, perr)
            # Extract the fitted parameters
            self.level_fundmental, self.level_hotband, self.background1, self.background2 = popt 

            initial_guess = [1600, 0] #[8.134e-06]
            popt, pcov = curve_fit(self.methane_spectrum_Roffset_fn, self.wave, self.spec, p0 = initial_guess, sigma = sigma)
            perr = np.sqrt(np.diag(pcov))
            print(popt, perr)
            self.wavenumber_offset = popt

            initial_guess = [self.level_fundamental, self.level_hotband, self.background1, self.background2]
            popt, _ = curve_fit(self.methane_spectrum_fn, self.wave, self.spec, p0 = initial_guess, sigma = sigma)
            print(popt)


        return self.methane_spectrum_fn(self.xpixels, *popt)


class fit_non_LTE_CH4_JWST: 
    '''
        Simple class to fit the non-LTE methane component between 3.2, and 3.6 microns in JWST NIRSpec spectra.     
    '''
    def __init__(self, wave, R = 1700) : 

        self.wave = wave

        # Generate the CH4 fundamental and hotband spectra from the HITRAN data
        ch4list = Table.read('data/extras/ch4_line_list.txt', format='ascii')
        ch4fit  = methane_fitting(ch4list, R = R)
        ch4fit.set_data(wave, wave, wave)
        model   = ch4fit.methane_spectrum_fn(wave, 1.0, 100.0, 0, 0, 0, 0)
        ch4_fun = ch4fit.methane_fundamental_fn()
        ch4_hot = ch4fit.methane_hotband_fn()

        # Normalise the spectral function
        self.ch4_fun = ch4_fun / np.nanmax(ch4_fun)
        self.ch4_hot = ch4_hot / np.nanmax(ch4_hot)

        # Set up a H3+ object
        self.h3p = h3ppy.h3p()
        self.h3p.set(T = 800, N = 1e16, R = 2700)

    def fit_function(self, wave, fun_scaling, hot_scaling, bkg1, bkg2, bkg3) : 
        '''
            This is the fitting function. It has a CH4 fundamental component, and a hotband component,
            as well as a polynomial (3rd order) background.     
        '''
        # Store the variable internally
        self.fun_scaling = fun_scaling
        self.hot_scaling = hot_scaling
#        self.wave_offset = wave_offset
        self.bkg1 = bkg1
        self.bkg2 = bkg2
        self.bkg3 = bkg3

        subwave = self.fit_wave # self.wave[self.whw]

        # Construct the CH4 spectral function 
        bkg = bkg1  + bkg2 * subwave + bkg3 * subwave**2
        ch4 = self.ch4_fun[self.whw] * fun_scaling + self.ch4_hot[self.whw] * hot_scaling
        spectrum = ch4 + bkg
        return spectrum

    def fit_function_bkg(self, wave, bkg1, bkg2, bkg3) : 
        subwave = self.fit_wave
        bkg     = bkg1  + bkg2 * subwave + bkg3 * subwave**2
        return bkg

    def fit_subregion(self, wave, spec) : 
        '''
            Fit an indivudial (limited) spectral region. Fitting non-LTE CH4 using LTE spectre
            only works well over a small wavelength range.      
        '''
        spec = np.nan_to_num(spec)

        # What wavelength region are we fitting? 
        _, self.whw, _ = np.intersect1d(self.wave, wave, return_indices = True)
#        self.whw = np.argwhere((self.wave >= np.min(wave)) & (self.wave <= np.max(wave)))
        self.fit_wave = wave

        # Work out the uncertainties - we want very large sigmas where the H3+ are located.
        # First, determine where the H3+ lines are 
        self.h3p.set(wave = wave)
        model         = self.h3p.model()
        model         = model / np.max(model) # Normalise the model

        wh_h3p        = np.argwhere(model > 0.05).flatten()
        # Set the sigma to be small where there's no H3+ lines
        sigma         = model.copy()
        sigma[:]      = 0.5 * np.nanmedian(spec)
        # Set the sigma to be large where there are H3+ lines
        sigma[wh_h3p] = 20.0 * np.nanmax(spec)
        self.sigma = sigma

        if (np.max(wave) < 3.5) :
            # Provide some intital guesses, then do the fitting
            initial_guess = [0.001, 0.0001, 0, 0, 0]
            self.vars, pcov = curve_fit(self.fit_function, wave, spec, sigma = sigma, p0 = initial_guess)
            # Re-calculate the best fit spectral function
            rawvars = self.vars
            print(rawvars)
            fit = self.fit_function(wave, *self.vars)
        else : 
            initial_guess = [0, 0, 0]
            self.vars, pcov = curve_fit(self.fit_function_bkg, wave, spec, sigma = sigma, p0 = initial_guess)
            # Re-calculate the best fit spectral function
            rawvars = self.vars
            fit = self.fit_function_bkg(wave, *self.vars)

        # Calculate the errors
        self.errs = np.sqrt(np.diag(pcov))

        bkg     = rawvars[-3]  + rawvars[-2] * wave + rawvars[-1] * wave**2

        return fit, bkg

    def fit(self, wave, spec, width = 0.1, start = 3.2, end = 3.6) : 
        '''
           Fit an extended wavelength range by splitting it up into chunks.   
        '''
        ret = spec.copy()
        ret[:] = np.nan
        bkg = spec.copy()
        bkg[:] = np.nan
        regions = np.arange(start, end, width)
        for r in regions : 
        # for i in range(5) :             
        #     if (i == 0) : whl = np.argwhere((wave > 3.29) & (wave < 3.31)).flatten()
        #     if (i == 1) : whl = np.argwhere((wave > 3.37) & (wave < 3.4)).flatten()
        #     if (i == 2) : whl = np.argwhere((wave > 3.4) & (wave < 3.43)).flatten()
        #     if (i == 3) : whl = np.argwhere((wave > 3.445) & (wave < 3.465)).flatten()
        #     if (i == 4) : whl = np.argwhere((wave > 3.526) & (wave < 3.56)).flatten()

            whw = np.argwhere((wave > r) & (wave < r + width)).flatten()            
            subwave = wave[whw]
            subspec = spec[whw]
            fit, bg = self.fit_subregion(subwave, subspec)
#            ret[whw] = ret[whw] - fit
            ret[whw] = fit
            bkg[whw] = bg
        return ret, bkg



