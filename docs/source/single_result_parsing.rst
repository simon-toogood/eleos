Parsing A Single NEMESIS Core
*****************************


Once a NEMESIS job has finished running, ``eleos`` will be invoked automatically to generate some summary plots. In order to create some custom plots we can use the ``eleos.results`` submodule. The :class:`eleos.results.NemesisResult` class holds all the information about the retrieval and contains some useful plotting routines. For a full lsit of everything the :class:`eleos.results.NemesisResult` contains, see the documentation page.


.. code-block:: python
   :linenos:
   
   import matplotlib.pyplot as plt
   from eleos import results

   # Load the retrieved core as a NemesisResult object
   res = results.NemesisResult("example_cores/core_1")

   # Get some information about the run
   res.print_summary()

   # Create two plots
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

   # Plot the retrieved spectrum and aerosol optical thickness as a function of pressure 
   res.plot_spectrum(ax=ax1)
   res.plot_aerosol_profiles(ax=ax2)

   # Save the figure
   fig.tight_layout()
   fig.savefig("example.png", dpi=500)

