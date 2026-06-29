.. libefp_ref:


************************
Running EFP calculations
************************

This section deals with performing pairwise EFP calculations on a unit cell of a benzene crystal.

Preparing LIBEFP input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

After obtaining the EFP parameters in the form of a .efp file from the MAKEFP job run,
we can prepare the LIBEFP input file to run an actual EFP calculation in the software library LIBEFP.
This what a template libefp input file looks like:

.. image:: ../images/libefp_inp.png
   :width: 360


Fragment name check
~~~~~~~~~~~~~~~~~~~

It is essential that the fragment name in the input file has a .efp file **of the same name**
in the same directory. The $FRAGNAME section in the header of the .efp file also needs to have the
**same** fragment name as that in the input file (**check video tutorial for detailed instructions**),
to ensure that the EFP calculation runs without any errors.




Running a calculation
~~~~~~~~~~~~~~~~~~~~~

To run an EFP calculation on the system, we use the command: efpmd <libefp_input_file_name>

Output
~~~~~~

The output of an EFP calculation typically gives the following output:
