import matplotlib.pyplot as plt
from eleos import results


# Load the retrieved core as a NemesisResult object
res = results.NemesisResult("example/core_1")

# Get some information about the run
res.print_summary()

# Create two plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

# Plot the retrieved spectrum and aerosol optical thickness as a function of pressure 
res.plot_spectrum(ax=ax1)
res.plot_aerosol_profiles(ax=ax2)

# Save the figure
fig.tight_layout()
fig.savefig("example.png", dpi=500)
