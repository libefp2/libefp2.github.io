.. _efp_calcs:

****************
EFP calculations
****************

EFP energy terms
----------------

LibEFP can compute four inter-fragment energy terms:

* electrostatic
* polarization
* dispersion
* exchange-repulsion

.. _elec_energy:

Electrostatic term
^^^^^^^^^^^^^^^^^^

Electrostatic energy is computed as a combination of charge-charge, charge-dipole,
charge-quadrupole, charge-octupole, dipole-dipole, dipole-quadrupole, and quadrupole-quadrupole
contributions.

The following parts of the EFP potential are used for electrostatic energy calculations:

- :ref:`COORD`
- :ref:`MON`
- :ref:`DIP`
- :ref:`QUAD`
- :ref:`OCT`
- :ref:`SCREEN`

Out of those, only :ref:`COORD` section is mandatory as it determines fragment internal coordinate frame.
All electrostatic terms (and corresponding parameters) are optional.

Charge-penetration screening
""""""""""""""""""""""""""""

Two options exist for accounting for charge-penetration contribution to electrostatic
energy:

- exponential screening ("smearing") of charges. This is achieved by invoking screening parameters
defined in ``SCREEN2`` section of the `.efp` potential (see :ref:`SCREEN`). The charge-penetration energy
is not printed separately but included in the electrostatic energy.
- overlap-based screening. This is a separate energy term derived assuming that localized orbitals
can be modeled as spherical gaussions (the same approximation is used in the exchange-repulsion term).
This calculation will utilize exchange-repulsion parameters (:ref:`FOCK`, :ref:`WF`, :ref:`BASIS`, :ref:`LMOC`).
This overlap-based charge-penetration energy is printed as a separate energy term (see examples in :ref:`libefp`).

Detailed description of damping functions and their benchmarks are published in
`Damping functions for electrostatic term <http://dx.doi.org/10.1002/jcc.20520>`_
and `Short-range damping functions <http://dx.doi.org/10.1080/00268970802712449>`_ papers.

.. _pol_energy:

Polarization term
^^^^^^^^^^^^^^^^^

Polarization energy is computed using distributed anisotropic dipole polarizabilities. Induced dipoles,
originating at the polarizability points, are converged until self-consistency. The default procedure is to
solve for induced dipoles iteratively; the direct diagonalization of the induced dipole matrix is implemented
but not parallelized, making its applicability limited to systems with a few thousands polarizability points
(see `polarization solver` section in `efpmd manual <https://github.com/libefp2/libefp/tree/master/efpmd#readme>`_).
Detailed description of the EFP polarization term can be found in
`the first EFP paper (1996) <https://doi.org/10.1063/1.472045>`_
and `gradients of polarization energy paper <https://doi.org/10.1063/1.2378767>`_.

The relevant sections of the EFP potential are:

- :ref:`POL_POINT`
- :ref:`POLAB`
- :ref:`COORD`
- :ref:`MON`
- :ref:`DIP`
- :ref:`QUAD`
- :ref:`OCT`

:ref:`POL_POINT` groups provides coordinates and values of the polarizability tensors. Other sections specify
positions and values of electrostatic multipoles that are used to compute static electric field on polarizability points.

Polarization energies are screened at short range with the Tang-Toennies (or gaussian-type) damping functions described in
the `short-range damping functions paper <http://dx.doi.org/10.1080/00268970802712449>`_. A value of the damping parameter is controlled
by an optional :ref:`POLAB` keyword; smaller values provide stronger screening of polarization energies which might be necessary for fragments
with large multiple moments (charged or strongly polar species) or large polarizabilities (e.g., large conjugated/aromatic molecules).

.. _disp_energy:

Dispersion term
^^^^^^^^^^^^^^^

Dispersion energy term captures the London interaction between the molecules. Formally, it can be expanded in
series of (1/R) operator as :math:`E_{disp} = \frac{C_6}{R^6} + \frac{C_8}{R^8} + \frac{C_{10}}{R^{10}} + ....`
In the case of distributed approach where dispersin contributions are computed as a sum of contributions due to
individual parts of a molecules, the odd terms :math:`\frac{C_7}{R^7}, \frac{C_9}{R^9}` etc are also non-zero.

The relevant sections of the EFP potential are:

- :ref:`DYN_POINT`

:ref:`DYN_POINT` group section provides coordinates and values of anisotropic dynamic polarizability tensors for computing dispersion energy.

.. _ex_rep:

Exchange Repulsion
^^^^^^^^^^^^^^^^^^

Exchange repulsion accounts for the antisymmetry of the wave function of the fragments.It is modelled using inter-fragment kinetic and
overlap integrals, and the Fock matrices of the fragment. 

The relevant sections of the EFP potential are:

- :ref:`BASIS`
- :ref:`MULTIPLICITY`
- :ref:`WF`
- :ref:`FOCK`
- :ref:`LMOC`

:ref:`BASIS` provides details of the basis set used for calculation of the exchange repulsion energy, :ref:`MULTIPLICITY` contains information
on the multiplicity of the fragment (LibEFP works only on fragments with multiplicity 1), :ref:`WF` provides the localized wave function of the 
fragment, while :ref:`FOCK` and :ref:`LMOC` contain information regarding the elements of the Fock matrix of the fragment in the localized basis, and 
the coordinates of the localized molecular orbital, respectively.
