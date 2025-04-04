from eleos import results

# Define the parent directory of the sensitivity analysis 
cd = "sensitivity/"

# Load the results
sens = results.SensitivityAnalysis(cd)

# Generate a plot for each parameter
sens.make_parameters_plot()

# Save the figure as 'parent directory/sensitivity.png'
sens.savefig("sensitivity.png", dpi=400)