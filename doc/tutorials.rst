.. _tutorials:

*********
Tutorials
*********

.. _non-covalent_calcs:

EFP calculations for a non-covalent system (molecular clusters, liquids, etc.)
-----------------------------------------------------------------------------------------

Step 1. Computing EFP parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EFP parameters need to be computed via `MAKEFP` calculation in the electronic structure software GAMESS.

* Setting up GAMESS MAKEFP job: see :ref:`compute parameters` and
  video tutorials `How to run MAKEFP job <https://youtu.be/orHN362tLjI?si=yGsAYYCeCTqhlyQt>`_ and
  `How to submit MAKEFP input in GAMESS <https://youtu.be/auM76y2tdzw?si=dBOa5sSojKmn4cc7>`_

* What are EFP parameters? See :ref:`efp parameters` and
  video tutorial `Analyzing MAKEFP output <https://youtu.be/R2r_IrV6NbY?si=tgePy5k_UIpmVMpl>`_

Step 2. LibEFP calculations of system energy and properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once EFP parameters of all unique fragments are computed, perform the full-system calculation in LibEFP:

* Preparing LIBEFP input file and running EFP calculation in LIBEFP. Refer to :ref:`running efp calculations` and
  video tutorial `Obtaining EFP parameters and running EFP calculations in LIBEFP <https://youtu.be/oww9uGJmKX4?si=Ih-lwMZTIWO3cK4z>`_


.. _covalent_calcs:

EFP calculations for a covalently bound systems (proteins, polymers)
--------------------------------------------------------------------

Tutorial `EFP protein-ligand calculations <https://github.com/libefp2/BioEFP-tools/tree/main/bioMAKEFP>`_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setting up QM/EFP calculations 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
See our new :ref:`applied_efp` tutorial!



Video Tutorials
---------------

*   `How to run MAKEFP job <https://youtu.be/orHN362tLjI?si=yGsAYYCeCTqhlyQt>`_

*   `How to submit MAKEFP input in GAMESS <https://youtu.be/auM76y2tdzw?si=dBOa5sSojKmn4cc7>`_

*   `Analyzing MAKEFP output <https://youtu.be/R2r_IrV6NbY?si=tgePy5k_UIpmVMpl>`_

*   `Obtaining EFP parameters and running EFP calculations in LIBEFP <https://youtu.be/oww9uGJmKX4?si=Ih-lwMZTIWO3cK4z>`_
