Documentation
*************


Core Generation
===============

.. autoclass:: eleos.NemesisCore
   :members:
   :class-doc-from: both

.. autofunction:: eleos.clear_parent_directory
.. autofunction:: eleos.generate_alice_job
.. autofunction:: eleos.reset_core_numbering
.. autofunction:: eleos.run_alice_job


Profiles
========

.. autoclass:: eleos.TemperatureProfile
    :members:
    :class-doc-from: both

.. autoclass:: eleos.GasProfile
    :members:
    :class-doc-from: both

.. autoclass:: eleos.AerosolProfile
    :members:
    :class-doc-from: both



Profile Shapes
==============

.. autoclass:: eleos.Shape0
    :members:
    :class-doc-from: both

.. autoclass:: eleos.Shape1
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape2
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape4
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape20
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape32
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape37
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape47
    :members:
    :class-doc-from: both
   
.. autoclass:: eleos.Shape48
    :members:
    :class-doc-from: both


Results
=======

.. autoclass:: eleos.NemesisResult
    :members:
    :class-doc-from: both

.. autofunction:: eleos.load_multiple_cores


Spectra Utilites
================
 
.. autofunction:: eleos.add_error_cubes
.. autofunction:: eleos.downsample_spectra
.. autofunction:: eleos.get_observations
.. autofunction:: eleos.margin_trim
.. autofunction:: eleos.multiple_cube_average
.. autofunction:: eleos.remove_nonlte_emission
.. autofunction:: eleos.trim_spectra