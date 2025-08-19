Generating Multiple NEMESIS Cores
*********************************


The process of generating multiple NEMESIS cores with ``eleos`` is much the same as generating a single core. This example will generate a set of cores, each with a different opacity for a stratospheric haze layer and run them as an array job on ALICE. Each core gets the same summary plots and tables created as in the first example.

.. literalinclude:: ../../examples/multiple_core_creation.py
   :language: python
   :linenos: