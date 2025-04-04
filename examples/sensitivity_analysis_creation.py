from eleos import cores


# Define and clear the new parent directory
cd = "sensitivity/"
# cores.clear_parent_directory(cd)

# Load the best core from a previous parameter space search
core = cores.load_from_previous("example/core_1/", cd)

# Create a sensitivity analysis and run it
cores.create_sensitivity_analysis(core)

# Use type_='sensitivity' to run the cleanup codea fter all cores have finished
# and generate plots
cores.generate_alice_job(cd, 
                         python_env_name="pythonmain", 
                         username="scat2", 
                         hours=2,
                         memory=1, 
                         type_="sensitivity")
cores.run_alice_job(cd)