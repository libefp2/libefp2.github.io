.. _learn_more:

********************
Learn more about EFP
********************

EFP is a polarizable force field for computing non-covalent interactions in a molecular system.
A solvent molecule, referred to as "fragment" or "effective fragment" is represented by a model potential (or force field)
with a set of parameters determined from a preparatory *ab initio* calculation.
There is no parameter fitting in EFP method. uch that the method is free of parameter fitting. Through its force field,
the effective fragments interact with each other and with an optional *ab initio* subsystem.


EFP calculations consist of two major steps:

* preparing effective parameters for all unique fragments of the system, as described in :ref:`parameters`
* evaluating interaction energy and properties of the molecular system, described in :ref:`efp_calcs`


EFP can be used for modeling of macromolecules and polymers, including proteins, which is described in :ref:`bioefp`.

EFP is combined with a number of *ab initio* methods for fast and efficient modeling of electronic structure in the condensed phases. See :ref:`QM/EFP`
section for details.
