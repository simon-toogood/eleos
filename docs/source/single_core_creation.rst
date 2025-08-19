Generating a single NEMESIS Core
********************************

In ``eleos``, each NEMESIS core is represented by a single :class:`eleos.cores.NemesisCore` object. The :class:`eleos.cores.NemesisCore` object contains all the information rewured to generate a core - for example, paths to the .ref and .spx files to use, number of iterations, scattering toggle, among others. 

In order to add a retrieved profile to the core, we can instantiate one of the subclasses of :class:`eleos.profiles.Profile`. There is currently support for two types of profiles; gas composition profiles (:class:`eleos.profiles.GasProfile`) and aerosol density profiles (:class:`eleos.profiles.AerosolProfile`), as well as very limited support for temperature profiles (:class:`eleos.profiles.TemperatureProfile`), though this is disabled by default.

Each of these profiles requires several parameters to set up fully, but an important one is the shape of the profile (ie the value of the VMR/density/etc as a function of altitude). These are represented by :class:`eleos.shapes.Shape` objects and their names correspond to the numbers in the NEMESIS manual - profile shape number 47 in NEMESIS is a Gaussian shape in log-pressure space which corresponds to :class:`eleos.shapes.Shape47` in ``eleos``.

The core in this example will have two retrieved gas profiles; ammonia and phosphine which will be represented by shapes 20 and 2 respectively. The indentation in this example is entirely optional, but in my opinion this is the easiest style to read.

.. code-block:: python
   :linenos:

   # Create a profile to retrieve ammonia - this is parameterised as model 20 in NEMESIS
   nh3 = profiles.GasProfile(gas_name="NH3", 
                             shape=shapes.Shape20(
                                 knee_pressure=1.0, 
                                 tropopause_pressure=0.1,
                                 deep_vmr=2e-4,           deep_vmr_error=10e-6,
                                 fsh=0.2,                 fsh_error=0.1))

   # Create a profile to retrieve phosphine - this is a simple scaling of the prior distribution
   ph3 = profiles.GasProfile(gas_name="PH3", 
                             shape=shapes.Shape2(
                                 scale_factor=2, scale_factor_error=0.1))

This core will also have a single aerosol profile representing a thick cloud deck in the tropopause, parameterised by model 48. The first snippet shows a cloud with set optical properties, but this can be changed simply by setting the ``retrieve_optical`` argument and specifying errors, as shown in the second snippet. Internally, this means ``eleos`` will create a 444 profile to retrive the particle radius, variance, and constant imaginary refractive index.


.. code-block:: python
   :linenos:

   main = profiles.AerosolProfile(label="Main Cloud", 
                                  retrieve_optical=False, 
                                  radius=3,    
                                  variance=0.5,
                                  imag_n=2e-2, 
                                  real_n=2.2,       
                                  shape=shapes.Shape48( 
                                     base_pressure=2.0,     base_pressure_error=0.2,
                                     top_pressure=0.1,      top_pressure_error=0.0001,
                                     opacity=22,            opacity_error=3,
                                     fsh=0.4,               fsh_error=0.03),)
      

.. code-block:: python
   :linenos:

   main = profiles.AerosolProfile(label="Main Cloud", 
                                  retrieve_optical=True, 
                                  radius=3,              radius_error=3,
                                  variance=0.5,          variance_error=0.5,
                                  imag_n=2e-2,           imag_n_error=2e-2,
                                  real_n=2.2,       
                                  shape=shapes.Shape48( 
                                      base_pressure=2.0,     base_pressure_error=0.2,
                                      top_pressure=0.1,      top_pressure_error=0.0001,
                                      opacity=22,            opacity_error=3,
                                      fsh=0.4,               fsh_error=0.03),)
      
Now that the profiles are defined, we can start to create the actual :class:`eleos.cores.NemesisCore` object. We begin by setting a directory for the core to reside in (relative to the working directory) and clearing it in case there are previous cores in there (don't worry - it will ask for confirmation before doing this!)

.. code-block:: python
   :linenos:

   cd = "example_cores/"
   cores.clear_parent_directory(cd)

Then we can instantiate :class:`eleos.cores.NemesisCore`. By default the core has multiple scattering enabled and is pre-configured to work with JWST/NIRSpec Jupiter observations. Other planets currently aren't supported but this will be a very simple addition. There are a set of default files in the ``eleos/data/jupiter`` which it pulls from by default, including (non-exhaustive) a 120-layer .ref file, .kls files for NIRSpec and MIRI (controlled by setting the `instrument` parameter), and a parah2.ref table. The ``reference_wavelength`` parameter gives the wavelength at which to normalise the aerosol cross-sections in order to get opacities, with the other parameters being fairly self-explanatory. For a full list of available options, see :meth:`eleos.cores.NemesisCore.__init__`.

.. code-block:: python
   :linenos:

   core = cores.NemesisCore(cd,
                            spx_file=f"data/spectra/n_80/noh3p/55.0S.spx",
                            profiles=[nh3, ph3, main],
                            fmerror_factor=2,
                            reference_wavelength=4,
                            num_layers=39,
                            min_pressure=1e-3,
                            max_pressure=10,
                            num_iterations=15)
   core.generate_core()

The final line tells ``eleos`` to create the directories and files necessary, ready for NEMESIS to be run. At this point the directory tree will look something like this (assuming this example code is in a file called ``generate.py`` in the working directory):

.. code-block:: 

   .
   |____generate.py
   |____example_cores/
   | |____core_1/
   | | |____plots/
   | | |____aerosol_names.txt
   | | |____aerosol.prf
   | | |____aerosol.ref
   | | |____cloudf1.dat
   | | |____core.pkl
   | | |____eleos_generation.py
   | | |____fcloud.prf
            ...

Inside ``example_cores/core_1`` are all the files required by NEMESIS, plus a few extra required by ``eleos``. Notable examples include ``eleos_generation.py`` (an exact copy of ``generate.py`` enabling the user to see the code which generated the core even if the original file is moved/deleted), ``aerosol_names.txt``  which gives the labels assigned to each :class:`eleos.profiles.AerosolProfile` object, and ``summary.txt`` which dumps all the attributes of the :class:`eleos.cores.NemesisCore` and all the associated profiles and shapes. Finally, to run NEMESIS with ``eleos`` we can add these lines to the bottom of the script:


.. code-block:: python
   :linenos:

   cores.generate_alice_job(cd, python_env_name="pythonmain", username="scat2", hours=12, memory=3)
   cores.run_alice_job(cd)

This will create a job submission script and give it to the scheduler. Afer NEMESIS has finished running, ``eleos`` will automatically create some summary plots for the core which are stored in the ``plots/`` directory in the core. There will also be a pretty-printed ASCII table of all the profile parameters in the slurm output which should appear at the top of the core directory in a file named like ``aa_slurm-XXXXX.out`` To work with the resulting retrieval in Python, see the next example.
