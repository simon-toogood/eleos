import results


res = results.NemesisResult("nemesis", "/home/s/scat2/NEMESIS/2022_JupSouthPole/core_1/")
fig, ax = res.plot_spectrum()
fig.savefig("spectrum.png", dpi=500)