import matplotlib.pyplot as plt

from nemesis_eleos import results


res = results.NemesisResult("cores/core_1/")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
res.plot_spectrum(ax=ax1)
res.plot_temperature(ax=ax2)
plt.tight_layout()

fig.savefig("nosync/temp.png", dpi=500)
