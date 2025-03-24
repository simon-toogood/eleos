from pathlib import Path
from eleos import cores


# Set the directory of the new core
cd = "high_res_forward"
cores.clear_parent_directory(cd)

# Load the previous core
core = cores.load_from_previous("low_res_retrieval/core_1", cd)

# Set the .spx file to a higher spectral resolution one
core.spx_file = Path("data/spectra/high_res.spx")

# Set the core to a forward model instead of a retrieval
core.forward = True

# Generate the new core with the profiles and parameters from the previous retrieval
core.generate_core()

# Run the new core on ALICE
cores.generate_alice_job(cd, python_env_name="pythonmain", username="scat2", hours=2, memory=3)
cores.run_alice_job(cd)