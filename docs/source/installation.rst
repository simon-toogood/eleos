Installation
*********************************


Installing ``eleos`` should be fairly straightforward. It does require NEMESIS to be installed, for which instructions can be found `here <https://github.com/nemesiscode/radtrancode>`_ and `here <https://github.com/nemesiscode/radtrancode/blob/master/AACOMPILE.txt>`_. Once this is installed successfully, ``eleos`` can be installed using ``pip`` (`PyPi page <https://pypi.org/project/nemesis-eleos>`_):

.. code-block:: 
   
   pip install nemesis-eleos


You can check if the installation was successful by running ``python -m eleos``. The output should be: 

.. code:: 

   Call signatures:

   python -m eleos <core_directory> --make_summary
   Print the summary tables and generate summary and iteration plots for the given core

   python -m eleos <parent_directory> --make-sensitivity-summary 
   For a directory containing a sensitivity analysis, create a plot showing the effect changing
   each parameter has on the spectrum

   You can append "--run-if-finished" to any command and this will run it only if
   all the cores in the parent_directory specified This is done by checking
   for the existence of "Done!" in the slurm output for each core.

   General command format:
         python -m eleos [DIRECTORY] [COMMAND] [RUN IF FINISHED] 
   argv:          0             1          2            3
