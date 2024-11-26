#========================================================================================
# Description
"""
Extracts and plots variables from an mre file

For optimal functioning of this code, organise the APR file as follows:
	- Temperatures (if a variable)
	- Aerosols
	- Gases

Usage
		[Input parameters below]

     		python plot_multimre_v2.py

Version 2
	- Updated the spectrum plotting code so it actually works
	- Input compatibility with new NEMESIS models
	- Improved the aerosol plotting code


"""
#========================================================================================
# Inputs


main_dir = 'core_1' # Directory where nemesis.mre is located

ref_skip = 25 # Lines to skip in the ref file (usually 25)

aero_unit = 'opacity/bar' # opacity/km or opacity/bar

plot_wave = True # Set to True to plot the spectrum

g = 23.12 # aceleration due to gravity ms-2
M = 2.22 # g/mol
R_gas = 8.314 # J/mol/K


#========================================================================================
# Imports


import numpy as np
import matplotlib.pyplot as plt
import os
import glob


#========================================================================================
# Functions


# Calculates the gas VMR profile for model 20
# I don't think this funtion does anything any more, which is good because I don't want to
def define_gas(PKNEE, PTOP, vmr, fsh, pressure):

	p_po = pressure / PKNEE

	profile = np.zeros((120))
	err = np.zeros((120))

	profile[np.where(p_po > 1.0)] = vmr
	profile[np.where(p_po < 1.0)] = vmr * p_po[np.where(p_po < 1.0)]**((1.0-fsh)/fsh)

	err[np.where(p_po > 1.0)] = dvmr
	err[np.where(p_po < 1.0)] = (p_po[np.where(p_po < 1.0)]**((1-fsh)/fsh) * dvmr) +\
	((-p_po[np.where(p_po < 1.0)]**(1/fsh) * np.log(p_po[np.where(p_po < 1.0)])) / (fsh**2 * p_po[np.where(p_po < 1.0)])) * dvmr * dfsh

	arral2 = np.where(pressure <= PTOP)

	profile[arral2] = 1e-36
	err[arral2] = 1e-36

	return profile, err


#========================================================================================
# Main


mre_file = main_dir + '/nemesis.mre'


# Sanity check - make sure the mre file actually exists

if os.path.exists(mre_file):
	print('\nMRE file found\n')
else:
	raise Exception('\nMRE file not found\n')


#---------- First reading of the mre file to obtain the number of wavelengths, lon and lat positions ------

with open(mre_file, 'r') as f:
	mre_data = f.readlines()
f.close()

if len(mre_data) < 5:
	raise Exception('MRE file is incomplete\nThe retrieval may have failed before completing\n')

for idx, line in enumerate(mre_data, 1):
	if "ispec,ngeom,ny,nx,ny" in line:
		id_array_str1 = mre_data[idx-1].split()
		id_array_str2 = mre_data[idx].split()

len_wave = int(id_array_str1[2])
lat = float(id_array_str2[0])
lon = float(id_array_str2[1])


#---------- Reading the log file to obtain chisq -----------

log_file = glob.glob(mre_file.replace('nemesis.mre', 'log_*'))[0]

with open(log_file, 'r') as f:
	log_data = f.readlines()
f.close()

chisq = float([s for s in log_data if 'chisq/ny is equal to :' in s][0].replace('chisq/ny is equal to :', ''))


# Plotting the wave

if plot_wave == True:
	print('Plotting spectral fit\n')

	wave_data = np.loadtxt(mre_file, skiprows=5, max_rows=len_wave, usecols=(1, 2, 3, 5))

	wave = wave_data[:,0]
	R_meas = wave_data[:,1]
	dR_meas = wave_data[:,2]
	R_fit = wave_data[:,3]


	fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 8))

	# Spectral plot
	ax[0].plot(wave, R_meas, color='black', linestyle='-', linewidth=0.7, label="Measurement")
	ax[0].fill_between(wave, R_meas-dR_meas, R_meas+dR_meas, color='black', alpha=0.3)
	ax[0].plot(wave, R_fit, color='tab:red', linestyle='-', linewidth=0.7, label="Model, $\chi^2$={aa:.2f}".format(aa=chisq))

	# Residual plot
	ax[1].plot(wave, R_fit-R_meas, color='black', linestyle='-', linewidth=0.7, label="Residual, $\chi^2$={aa:.2f}".format(aa=chisq))
	ax[1].fill_between(wave, -dR_meas, dR_meas, color='black', alpha=0.3)

	ax[0].set_xlabel('Wavelength (µm)')
	ax[0].set_ylabel('Radiance (µWcm$^{-2}$sr$^{-1}$µm$^{-1}$)')
	ax[0].set_title(r"$\bf{" + "a" + "}$" + " Spectral Data", loc='left')

	ax[0].grid()
	ax[0].legend(prop={'size': 8})

	ax[1].set_xlabel('Wavelength (µm)')
	ax[1].set_ylabel('Residual Radiance (µWcm$^{-2}$sr$^{-1}$µm$^{-1}$)')
	ax[1].set_title(r"$\bf{" + "b" + "}$" + " Residual Data", loc='left')

	ax[1].grid()
	ax[1].legend(prop={'size': 8})

	fig.suptitle("Spectral Data ({aa:.1f}$\degree$W, {bb:.1f}$\degree$S)".format(aa=lon, bb=abs(lat)))
	plt.subplots_adjust(wspace=0.2)

	plt.show()


# ------ Plotting the retrieved variables-------


nvar = int([s for s in mre_data if 'nvar =' in s][0].replace('nvar =', ''))

print('{} retrieved variables detected\n'.format(nvar))

pressure = np.loadtxt(main_dir + '/nemesis.ref', skiprows=ref_skip, usecols=(1)) # atm
pressure *= 1013.25 # mbar
pressure_pa = pressure * 100 # pa


# In case temperature is not a variable
temp_data = np.loadtxt(main_dir + '/tempapr.dat', skiprows=1, usecols=(1, 2))
temp_ret = temp_data[:,0]
dtemp_ret = temp_data[:,1]


M /= 1000 # kg/mol
g *= 100 # cm s-2


for idx, line in enumerate(mre_data, 1):
	if "Variable" in line:
		id_array_str = mre_data[idx].split()
		id_array = [int(numeric_string) for numeric_string in id_array_str]

		print("\n{aa} {bb} {cc} detected".format(aa=id_array[0], bb=id_array[1], cc=id_array[2]))

		if id_array[0] == 0:
			print('Temperature profile detected')

			temp_data = np.loadtxt(mre_file, skiprows=idx+3, max_rows=120, usecols=(2, 3, 4, 5))
			temp_prior = temp_data[:,0]
			dtemp_prior = temp_data[:,1]
			temp_ret = temp_data[:,2]
			dtemp_ret = temp_data[:,3]

			fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

			ax.plot(temp_prior, pressure, color='black', linestyle='-', linewidth=1.0, label='Prior')
			ax.fill_betweenx(pressure, temp_prior-dtemp_prior, temp_prior+dtemp_prior, color='black', alpha=0.3)
			ax.plot(temp_ret, pressure, color='tab:red', linestyle='-', linewidth=1.0, label='Retrieved')
			ax.fill_betweenx(pressure, temp_ret-dtemp_ret, temp_ret+dtemp_ret, color='tab:red', alpha=0.3)

			ax.set_xlabel('Temperature (K)')
			ax.set_ylabel('Pressure (mbar)')

			ax.set_yscale('log')
			ax.set_ylim(max(pressure), min(pressure))

			ax.legend()
			ax.grid()

			plt.show()


		if id_array[0] == -1:
			print('Aerosol profile detected')

			with open(main_dir + '/aerosol.prf', 'r') as f:
				aero_data = f.readlines()
			f.close()

			aero_array_str = aero_data[1].split()
			naero = int(aero_array_str[1])
			print("{} layers detected\n".format(naero))

			aero_data = np.loadtxt(main_dir + '/aerosol.prf', skiprows=2)
			altitude = aero_data[:,0] # km

			fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

			colours = ['#648FFF', '#FFB001', '#DD2680', '#FE6100']

			#total_aero = np.zeros((120))

			for kk in range(naero):
				aero = aero_data[:,kk+1] # cm2/g

				if aero_unit == 'opacity/km':
					rho = (M * pressure_pa) / (R_gas * temp_ret) # Kg/m3
					rho /= 1000 # g/cm3

					aero_new = (10**5 * aero * rho) # tau/km

					xlab = r'$\tau$km$^{-1}$'

				if aero_unit == 'opacity/bar':
					aero_new = (10**5 * aero) / g # tau/bar

					# Oliver - please don't worry about this correction too much
					# I promise it works!
					aero_new /= 0.1 # uNhInGeD Correction for g units

					xlab = r'$\tau$bar$^{-1}$'

				if aero_unit != 'opacity/km' and aero_unit != 'opacity/bar':
					raise Exception('Invalid units for aerosols')

				ax.plot(aero_new, pressure, color=colours[kk], linestyle='-', label='Profile {}'.format(kk+1))
				ax.fill_between(aero_new, pressure, color=colours[kk], alpha=0.6)

				#total_aero += aero_new

			#ax.plot(total_aero, pressure, color='#775EF0', linestyle='-', label='Total'.format(kk+1))

			ax.set_xlabel(xlab)
			ax.set_ylabel('Pressure (mbar)')

			ax.set_yscale('log')
			ax.set_ylim(3000, 80)

			ax.grid()
			ax.legend()

			plt.show()

		if id_array[0] == 444:
			print('Model 444 detected for a cloud layer')
			print('Cloud layer {} has been modified\n'.format(id_array[1]))

			aero_array_str = mre_data[idx+3].split()
			aero_array = [float(numeric_string) for numeric_string in aero_array_str]
			print('Variable        Prior        Retrieved')
			print('Radius       {aa:.2f} ± {bb:.2f}    {cc:.2f} ± {dd:.2f}'.format(aa=aero_array[2], bb=aero_array[3], cc=aero_array[4], dd=aero_array[5]))

			aero_array_str = mre_data[idx+4].split()
			aero_array = [float(numeric_string) for numeric_string in aero_array_str]
			print('Variance     {aa:.2f} ± {bb:.2f}    {cc:.2f} ± {dd:.2f}'.format(aa=aero_array[2], bb=aero_array[3], cc=aero_array[4], dd=aero_array[5]))

			aero_array_str = mre_data[idx+5].split()
			aero_array = [float(numeric_string) for numeric_string in aero_array_str]
			print('Im idx      {aa:.3f} ± {bb:.3f}  {cc:.3f} ± {dd:.3f}'.format(aa=aero_array[2], bb=aero_array[3], cc=aero_array[4], dd=aero_array[5]))


		if id_array[2] == 3 or id_array[2] == 20:
			print('Gas profile detected\n')

			with open(main_dir + '/nemesis.prf', 'r') as f:
				prf_data = f.readlines()
			f.close()

			for idx, line in enumerate(prf_data, 1):
				if " {} ".format(id_array[0]) in line:
					counter = idx

					pro_prior = np.loadtxt(main_dir + '/nemesis.ref', skiprows=ref_skip, usecols=(counter)) * 1e+6 # ppm
					pro_ret = np.loadtxt(main_dir + '/nemesis.prf', skiprows=ref_skip-1, usecols=(counter)) * 1e+6 # ppm

					fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

					ax.plot(pro_prior, pressure, color='black', linestyle='-', label='Prior')
					ax.plot(pro_ret, pressure, color='tab:red', linestyle='-', label='Retrieved')

					ax.set_xlabel('VMR (ppm)')
					ax.set_ylabel('Pressure (mbar)')
					ax.set_title('ID {aa} {bb} {cc}'.format(aa=id_array[0], bb=id_array[1], cc=id_array[2]))

					ax.set_yscale('log')

					if id_array[0] == 26 or id_array[0] == 27:
						ax.set_ylim(100, 0.001)
					else:
						ax.set_ylim(1500, 80)

					ax.grid()
					ax.legend()

					plt.show()

					break


#========================================================================================

print('\nEnd of script\n')
