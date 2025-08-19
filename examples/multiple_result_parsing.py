import matplotlib.pyplot as plt
from eleos import results


# Load all the cores - use raise_errors=False to silently skip any failed retrievals
res = results.load_multiple_cores("multiple_example/", raise_errors=False)

# Loop over all the NemesisResult objects
for r in res:

    # Extract the AerosolProfile that was used to represent the haze
    haze = r.profiles["Haze"]

    # Print out the prior and retrieved opacities, and the chi squared value of the fit
    print(haze.opacity, haze.retrieved_opacity, r.chi_sq)

    # Plot the model spectrum
    plt.plot(r.retrieved_spectrum.wavelength, r.retrieved_Spectrum.model, label=f"Haze opacity={haze.opacity}")

# Plot the data for comparison
plt.plot(res[0].retrieved_spectrum.wavelength, res[0].retrieved_spectrum.measured, label="Data")
plt.legend()
plt.show()
