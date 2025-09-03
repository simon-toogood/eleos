import argparse
import subprocess
from pathlib import Path

from . import results


def exit_msg(msg="Not all cores have finished running!"):
    print(msg)
    exit()


def check_all_finished(directory):
    for core in directory.glob("**/nemesis.apr"):
        try:
            slurm = tuple((core / "..").resolve().glob("*slurm*.out"))[0]
            if "Done!" not in open(slurm).read():
                exit_msg()
        except IndexError:
            exit_msg()
    print("All cores have run! Running command...")


def make_summary(core_dir, log=False, silent=False):
    res = results.NemesisResult(core_dir)
    if not silent:
        res.print_summary()
    res.make_summary_plot(log=log)
    if not res.core.forward:
        res.make_iterations_plot()
    print("Done!")


def make_sensitivity_summary(parent_dir):
    sens = results.SensitivityAnalysis(parent_dir)
    sens.make_parameters_plot()
    sens.savefig("sensitivity.png", dpi=400)
    print("Done!")


parser = argparse.ArgumentParser(
    description="Run eleos analysis commands."
)
subparsers = parser.add_subparsers(dest="command", required=True)

# make-summary command
summary_parser = subparsers.add_parser("make-summary", help="Print summary tables and generate plots for a core")
summary_parser.add_argument("directory",         help="Path to a core directory",)
summary_parser.add_argument("--run-if-finished", help="Run only if all cores in directory are finished",       action="store_true",)
summary_parser.add_argument("--log-spectrum",    help="Whether to plot the spectrum/residuals on a log scale", action="store_true",)
summary_parser.add_argument("--silent",          help="Whether to print a table of the priors/fitted parameters", action="store_true",)

# make-sensitivity-summary command
sens_parser = subparsers.add_parser("make-sensitivity-summary", help="Create a plot showing effect of changing parameters on the spectrum")
sens_parser.add_argument("directory",         help="Path to a parent directory containing sensitivity analysis",)
sens_parser.add_argument("--run-if-finished", help="Run only if all cores in directory are finished", action="store_true")


args = parser.parse_args()
directory = Path(args.directory)

if args.run_if_finished:
    check_all_finished(directory)

if args.command == "make-summary":
    make_summary(directory, log=args.log_spectrum, silent=args.silent)

elif args.command == "make-sensitivity-summary":
    make_sensitivity_summary(directory)

else:
    parser.error("Command not recognised!")

