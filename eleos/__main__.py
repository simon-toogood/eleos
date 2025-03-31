import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from . import results


if sys.argv[1] == "--make-summary":
    res = results.NemesisResult(sys.argv[2])
    res.print_summary()
    res.make_summary_plot()
    if not res.core.forward:
        res.make_iterations_plot()

elif sys.argv[1] == "--make-sensitivity-summary":
    sens = results.SensitivityAnalysis(sys.argv[2])
    sens.make_parameters_plot()
    sens.savefig("sensitivity.png", dpi=400)
    
else:
    raise ValueError("Argument not recognised!")
    

print("Done!")

