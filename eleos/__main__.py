import argparse
import re
from pathlib import Path

from . import results, utils



def make_summary(args):
    res = results.NemesisResult(args.directory)

    if not args.silent:
        res.print_summary()
    res.make_summary_plot(log=args.log_spectrum)
    if not res.core.forward:
        res.make_iterations_plot()
    print("Done!")


def make_sensitivity_summary(args):
    sens = results.SensitivityAnalysis(args.directory)
    for mode in ("absolute", "difference", "percentage"):
        sens.make_parameters_plot(mode=mode)
        sens.savefig(f"sensitivity_{mode}.png", dpi=400)
    print("Done!")


def find_best_fit(args):
    out = []
    for fp in args.directory.glob("**/nemesis.chi"):
        with open(fp) as file:
            n = utils.get_floats_from_string(file.readlines()[-1])[0]
            id_ = int(re.search(r"core_(\d+)", str(fp)).group(1))
            out.append((id_, n))
    out = list(sorted(out, key=lambda x: x[1], reverse=True))
    for x in out:
        print(*x)


parser = argparse.ArgumentParser(description="Run eleos analysis commands.")
subparsers = parser.add_subparsers(dest="command", required=True)

# make-summary command
p = subparsers.add_parser("make-summary", help="Print summary tables and generate plots for a core")
p.set_defaults(func=make_summary)
p.add_argument("directory",         help="Path to a core directory",                                 type=Path)
p.add_argument("--log-spectrum",    help="Whether to plot the spectrum/residuals on a log scale",    action="store_true",)
p.add_argument("--silent",          help="Whether to print a table of the priors/fitted parameters", action="store_true",)

# make-sensitivity-summary command
p = subparsers.add_parser("make-sensitivity-summary", help="Create a summary for a sensitivity analysis")
p.set_defaults(func=make_sensitivity_summary)
p.add_argument("directory",         help="Path to a parent directory containing sensitivity analysis", type=Path)

# find-best-fit command
p = subparsers.add_parser("find-best-fit", help="Print all cores and their final chi-squared values in order from worst to best fit")
p.set_defaults(func=find_best_fit)
p.add_argument("directory",         help="Path to a core directory",                        type=Path)


args = parser.parse_args()
args.func(args)