import sys
from pathlib import Path
import subprocess

from . import results


"""Call signatures:

python -m eleos --make_summary <core_directory>
Print the summary tables and generate summary and iteration plots for the given core

python -m eleos --make-sensitivity-summary <parent_directory>
For a directory containing a sensitivity analysis, create a plot showing the effect changing
each parameter has on the spectrum

You can append "--run-if-finished" to any command and this will run it only if
all the cores in the parent_directory specified This is done by checking
for the existance of "Done!" in the slurm output for each core.

General command format:
       python -m eleos [DIRECTORY] [COMMAND] [RUN IF FINISHED] 
argv:          0             1          2            3

"""


def exit_msg(msg="Not all cores have finished running!"):
    print(msg)
    exit()


if sys.argv[-1] == "--run-if-finished":
    to_check = Path(sys.argv[1])

    # Pick a file present in all cores whether finished, running or not run yet
    for core in to_check.glob("**/nemesis.apr"):
        try:
            slurm = tuple((core / "..").resolve().glob("*slurm*.out"))[0]
            if "Done!" not in open(slurm).read():
                exit_msg()
        except IndexError:
            exit_msg()
    
    # If we haven't exited yet, all cores must have run.
    # There could be some weird race condition here when multiple cores finish at exactly the
    # but I'm choosing to ignore that and hope it's not a problem
    print(f"All cores have run! Running command...")
    subprocess.run(f"python -m eleos {to_check.resolve()} {sys.argv[2]}", shell=True)
    exit()


if sys.argv[2] == "--make-summary":
    res = results.NemesisResult(sys.argv[1])
    res.print_summary()
    res.make_summary_plot()
    if not res.core.forward:
        res.make_iterations_plot()
    print("Done!")

elif sys.argv[2] == "--make-sensitivity-summary":
    sens = results.SensitivityAnalysis(sys.argv[1])
    sens.make_parameters_plot()
    sens.savefig("sensitivity.png", dpi=400)
    print("Done!") 

else:
    raise ValueError("Argument not recognised!")
    


