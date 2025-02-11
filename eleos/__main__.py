import sys

from . import results


if sys.argv[1] == "--make-summary":
    res = results.NemesisResult(sys.argv[2])
    res.print_summary()
    res.make_summary_plot()
    if not res.core.forward:
        res.make_iterations_plot()
    print("Done!")

else:
    raise ValueError("Argument not recognised!")
    