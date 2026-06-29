.. _efpmd:

************************
LIBEFP/EFPMD user manual
************************

.. _general functionality:

General functionality
---------------------

*EfpMD* is a molecular simulation program for *LibEFP*. It supports EFP-only (no QM region) 
calculations. The following functionality is currently available:

- single point energy 
- gradient 
- geometry optimization
- semi-numerical Hessian and normal mode analysis 
- molecular dynamics simulations in microcanonical (NVE), canonical (NVT), and isobaric-isothermal (NPT) ensembles
  
Additioanally, *EfpMD* provides calculations of

- electrostatic potential on fragment atoms
- electric field on fragment atoms
- pairwise interaction energy decomposition analysis (PIEDA) 

*EfpMD* supports calculations with periodic boundary conditions in arbitrary parallelopiped unit 
cells. Switching functions can be used for treating long-range interactions. 

Since version 2.0, *EfpMD* can perform energy and geoemtry optimization calculations 
utilizing neural network (NN) potentials through interface with `LibTorch <https://pytorch.org/cppdocs/>`_. 

Additional examples of input files can be found in the `tests` directory in
source code archive.

.. _input file format:

Input file format
-----------------

The input file generally consists of two parts: a part containing parameters of the 
requested calculation (:ref:`efpmd parameters`), and the second part containing 
coordinates of a molecular system (:ref:`fragment coordinates`). 

An example *EfpMD* input file for performing geometry optimization of a 4-water cluster 
:download:`opt_1.in <../examples/opt_1.in>`.

.. literalinclude:: ../examples/opt_1.in
  :linenos:

.. note:: Lines beginning with the ``#`` symbol are ignored during input parsing.

.. _efpmd parameters:

EFPMD parameters
-----------------

.. _Basic parameters:

Basic parameters
^^^^^^^^^^^^^^^^^

**Type of the simulation**

``run_type [sp|grad|hess|opt|md|efield|elpot|gtest|etest]``  

- `sp` - single point energy calculation.
- `grad` - energy gradient calculation.
- `hess` - semi-numerical Hessian calculation and normal mode analysis.
- `opt` - geometry optimization.
- `md` - molecular dynamics simulation.
- `efield` - compute and print electric field on all atoms.
- `elpot` - compute and print electrostatic potential on all atoms.
- `frag_elpot` - compute and print electrostatic potential on all atoms of a ``special`` fragment.
- `gtest` - compute and compare numerical and analytical gradients.
- `etest` - compute and compare total energy.   

Default value: ``sp``

**Format of fragment input**

``coord [xyzabc|points|rotmat|atoms]``

- `xyzabc` - Coordinates of the center of mass and Euler angles.
- `points` - Coordinates of three atoms for each fragment.
- `rotmat` - Coordinates of the center of mass and rotation matrix.
- `atoms` - Atom names and cartesian coordinates of all atoms.

Default value: ``points``

.. note:: See :ref:`fragment input` specification for more details. 
.. warning:: Default changed to ``points`` from ``xyzabc`` in version 1.7.3

**EFP energy terms**

``terms [elec pol disp xr qq lj]``

- `elec` - Include electrostatics energy.
- `pol` - Include polarization energy.
- `disp` - Include dispersion energy.
- `xr` - Include exchange repulsion energy.
- `qq` - Include point-charge Coulomb energy due to MM charges 
  (``MM_CHARGE`` section in `.efp` file is required).
- `lj` - Include Lennard-Jones interaction term (``MM_LJ`` section in `.efp` file is required).

Default value: ``elec pol disp xr``

.. _Short-range damping parameters:

Short-range damping parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Electrostatic damping type**

``elec_damp [screen|overlap|off]``

- `screen` - Damping formula based on ``SCREEN`` group in the EFP potential.
- `overlap` - Overlap-based damping formula. This damping correction is printed
as charge penetration energy.
- `off` - No electrostatic damping.

Default value: ``screen``

**Dispersion damping type**

``disp_damp [tt|overlap|off]``

- `tt` - Damping based on the formula by Tang and Toennies, with hard-coded coefficient :math:`b = 1.5`.
- `overlap` - Overlap-based dispersion damping.
- `off` - No dispersion damping.

Default value: ``overlap``

**Polarization damping type**

``pol_damp [tt|off]``

- `tt` - Tang and Toennies like damping formula.
- `off` - No polarization damping.

Default value: ``tt``



.. _long-range:

Long-range cutoffs
^^^^^^^^^^^^^^^^^^^^

**Enable cutoff for fragment/fragment interactions**

``enable_cutoff [true|false]``

Default value: ``false``

**Cutoff distance for fragment/fragment interactions**

``swf_cutoff <value>``

Default value: ``10.0``

Unit: Angstrom

**Cutoff distance for exchange-repulsion interactions between fragments**

``xr_cutoff <value>``

Default value: ``swf_cutoff``

Unit: Angstrom

.. _polarization solver:

Polarization solver
^^^^^^^^^^^^^^^^^^^^^

**Polarization solver**

``pol_driver [iterative|direct]``

- `iterative` - Iterative solution of system of linear equations for polarization
  induced dipoles.
- `direct` - Direct solution of system of linear equations for polarization 
  induced dipoles. This solver does not have convergence issues but is unsuitable
  for large systems (more than 2000 polarizable points). The direct solver is not 
  parallelized.

Default value: ``iterative``


.. _geometry optimization:

Geometry optimization 
^^^^^^^^^^^^^^^^^^^^^^

**Maximum number of optimization and MD steps**

``max_steps <number>``

Default value: ``100``

.. note:: ``max_steps`` specifies maximum number of steps for both geometry optimization and MD.

**Gradient tolerance**

``opt_tol <value>``

Default value: ``0.0003``

Unit: Hartree/Bohr

**Energy tolerance**

``opt_energy_tol <value>``

Default value: ``0.000001``

Unit: Hartree

.. note:: Optimization will stop when maximum gradient component is less than ``opt_tol`` 
	and energy change is less than ``opt_energy_tol``.


.. _hessian:

Hessian calculation
^^^^^^^^^^^^^^^^^^^^

**Hessian accuracy**

``hess_central [true|false]``

Default value: ``false``

If ``hess_central`` is ``true``, then the more accurate central differences will be
used for numerical Hessian calculation. Otherwise, only forward differences will
be used. Note that central differences require twice as many gradient
calculations.

**Numerical differentiation step length for distances**

``num_step_dist <value>``

Default value: ``0.001``

Unit: Angstrom

**Numerical differentiation step length for angles**

``num_step_angle <value>``

Default value: `0.01`

Unit: Radian


.. _MD:

Molecular dynamics simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Ensemble**

``ensemble [nve|nvt|npt]``

- `nve` - Microcanonical ensemble.
- `nvt` - Canonical ensemble with Nose-Hoover thermostat. For the description of 
  the algorithm, see _Phys. Rev. A 31, 1695 (1985)_.
- `npt` - Isobaric-isothermal ensemble. This also sets ``enable_pbc`` to ```true``. 
  For the description of the algorithm, see _Mol. Phys. 78, 533 (1993)_.

Default value: ``nve``

**Time step**

``time_step <value>``

Unit: Femtosecond

Default value: ``1.0``


**Print step:** Number of steps between outputs of the system state.

``print_step <value>``

Default value: ``1``


**Assign initial velocities**

``velocitize [true|false]``

Default value: ``false``

If ``true``, random initial velocities will be assigned to fragments using
Gaussian distribution. Velocity magnitudes are chosen so that initial
temperature of the system is approximately equal to the target simulation
temperature.

**Simulation temperature:** target simulation temperature for NVT and NPT thermostats.

``temperature <value>``

Unit: Kelvin

Default value: ``300.0``


**Simulation pressure**

Target simulation pressure for NPT barostat. Note that whether or not pressure
coupling is used, the pressure value will oscillate significantly. Fluctuations
on the order of hundreds of bar are typical. This variation is entirely normal
due to the fact that pressure is a macroscopic property and can only be
measured properly as time average. The actual variations of instantaneous
pressure depend on the size of the system and the values of barostat
parameters.

``pressure <value>``

Unit: Bar

Default value: ``1.0``


**Thermostat parameter:** temperature relaxation time parameter of the Nose-Hoover thermostat

``thermostat_tau <value>``

Unit: Femtosecond

Default value: ``1000.0``


**Barostat parameter:** pressure relaxation time parameter of the barostat

``barostat_tau <value>``

Units: Femtosecond

Default value: ``10000.0``

.. _different


**Enable multistep MD**

``enable_multistep [true|false]``

Default value: ``false``

**Number of steps between recomputation of slow terms in multistep MD**

``multistep_steps <number>``

Default value: ``1``

Currently only exchange-repulsion EFP term is affected.


.. _PBC:

Periodic Boundary Conditions (PBC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Enable/Disable PBC**

``enable_pbc [true|false]``

Default value: ``false``

.. note:: Setting ``enable_pbc`` to ``true`` also sets ``enable_cutoff`` to ``true``.

**Dimensions of periodic cell** 

Units: Angstroms, degrees.

``periodic_box <x> <y> <z> <alpha> <beta> <gamma>``

Default value: ``30.0 30.0 30.0 90.0 90.0 90.0``

.. note:: If only three values are given, the angles are set to 90 degrees (orthogonal box). 
	Non-orthogonal PBC are implemented only for single-point energy calculations.
	The smallest box dimension must be greater than `2 * swf_cutoff`.

**Print PBC coordinates** 

Prints coordinates of the system contained in a single periodic cell around 
a fragment specified by ``ligand`` keyword.

``print_pbc [true/false]``

Default value: ``false``


.. _PIEDA:

Pairwise interaction energy decomposition analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Enable/Disable pairwise analysis**

The PIEDA analysis will be performed with respect to the fragment specified by the ``ligand`` keyword. 

``enable_pairwise [true|false]``

Default value: ``false``

**Specify ligand**

``ligand [integer]``

Default value: ``0`` (the first fragment in the system)


.. _Lattice symmetry:

Symmetric molecular crystals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Use symmetry** 

If ``true``, effectively performs calculations only on symmetry-unique fragments, which speeds up
calculations of symmetric crystal systems with PBC.
.. warning:: Implemented for single-point energy calculations. Not parallelized. 
See ``symm_frag`` keyword for specifying symmetry-identical fragments.

``symmetry [true|false]``

Default value: ``false``

**Specifying symmetry-identical fragments**

``symm_frag [frag | list]``

- `frag` - assumes that all fragments of the same type are identical.
- `list` - fragments provided in a list are symmetric (not fully implemented).

Default value: ``frag``

.. _libtorch:

Calculations with NN potentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Starting from version 2.0, calculations with NN potentials are enabled through interface 
with `LibTorch` library. Standard ANI potentials (`TorchANI <https://aiqm.github.io/torchani/index.html>`_) 
and custom models with embedded electrostatic 
potentials (see ``enable_elpot``) can be used. 
Single point energy and geometry optimization are currently implemented, see ``opt_special_frag`` for detail.

**Enable NN calculations**

``enable_torch [true|false]``

Default value: ``false``

**Specify NN fragment:** no default value, need to specify the fragment to be treated with NN.

``special_fragment <value>``

Default value: ``none``

**Enable electrostatic potential embedding in NN calculations, i.e., NN/EFP**

``enable_elpot [true|false]``

- `true` - provide custom NN using ``custom_nn`` and ``aev_nn`` 
- `false` - provide one of ANI potentials with ``torch_nn``

Default value: ``false``


**Optimization protocol in NN calculations**

``opt_special_frag [0|1]``

- 0 - optimize special (NN) fragment only
- 1 - optimize the whole system

Default value: ``0``

**Atomic gradients in NN calculations**

``atom_gradient [mm|frag]``

- `mm` - use direct atomic gradients from MM model
- `frag` - gradient on special fragment atoms computed by distributing the COM force/torque gradient

Default value: ``frag``

.. note:: `mm` algorithm can be used only with MM (``qq`` and ``lj``) potentials.

.. note:: `frag` algorithm distributes EFP center of mass force and torque into gradients on atoms. This 
	decomposition is not unique resulting in non-exact gradients. If you experience convergence issues 
	in geometry optimization, loosen optimization tolerance criteria ``opt_tol``.

**Name of the NN potential for standard ANI models**

``torch_nn <value>``

Default value: ``ani.pt``

**Name of custom NN potential for NN/EFP embedded model**

``custom_nn <value>``

Default value: ``custom_model_script.pt``

**Name of file containing atomic environment vector (AEV) for NN/EFP embedded model**

``aev_nn <value>``

Default value: ``aev_scripted.pt``

**Path to the directory with provided NN potentials**

``ml_path <value>``

Default value: `"$(prefix)/share/libefp"` (data install directory)

**Path to the directory with user NN potentials**

``userml_path <value>``

Default value: `"."` (current directory)

.. _other keywords:

Various other keywords
^^^^^^^^^^^^^^^^^^^^^^^^

**Print additional information for debugging**

``print <value>``

- 0 - standard output 
- 1 - detailed information on gradients 
- 2 - additional information on fragments; polarization convergence
- 3 - various information for debugging

Default value: ``0``

**Use single EFP parameters file**

``single_params_file [true|false]``

Default value: ``false``

**Single EFP parameters file path**

``efp_params_file <path>``

Default value: ``params.efp``

**The path to the directory with fragment library**

``fraglib_path <path>``

Default value: `"$(prefix)/share/libefp"` (data install directory)

The ``<path>`` parameter should not contain spaces or should be inside double
quotes otherwise.

**The path to the directory with user-created fragments**

``userlib_path <path>``

Default value: ``"."`` (current directory)

The ``<path>`` parameter should not contain spaces or should be in double quotes
otherwise.



.. _test jobs:

Parameters for test jobs
^^^^^^^^^^^^^^^^^^^^^^^^^

See also ``num_step_dist`` and ``num_step_angle`` in :ref:`hessian` section for 
controlling step size for numerical gradient testing.  

**Test tolerance**

``gtest_tol <value>``

Default value: ``1.0e-6``

Unit: Hartree/Bohr

**Reference energy value**

``ref_energy <value>``

Default value: ``0.0``

Unit: Hartree

.. _old:

Depricated Parameters
^^^^^^^^^^^^^^^^^^^^^^

**Geometry of the molecular-mechanics part**

``ff_geometry <path>``

Default value: ``ff.xyz``

**Molecular-mechanics force-field parameters file**

``ff_parameters <path>``

Default value: ``fraglib/params/amber99.prm``

**Enable molecular-mechanics force-field for flexible EFP links**

``enable_ff [true|false]``

Default value: ``false``


.. _fragment coordinates:

Molecular geometry
-------------------

.. _fragment name:

Fragment name
^^^^^^^^^^^^^^^

`LibEFP/EfpMD` distribution comes with a library of fragment potentials used for testing and debugging 
purposes. Those fragments are stored in the ``fraglib`` directory. 

.. warning:: `.efp` potentials from ``fraglib`` directory should not be used for any practical calculations 
	as they might produce erroneous results; 
	they are used for code testing only. 

.. note:: In order to use the `.efp` potentials from ``fraglib`` directory (defined by ``fraglib_path`` 
	parameter), the fragment ``<name>`` should contain an ``_l`` suffix, e.g., `h2o_l`, `c6h6_l`, etc.

User-defined parameters (or parameters from :ref:`efp_parameter_databases`) do not need 
suffix ``_l``, e.g., `h2o`, `c6h6`. Location of the parameter files is defined by 
``userlib_path`` keyword, with the default being the working directory ``./``. 

The name of the `.efp` file must be the same as the name of the fragment (without an ``_l`` suffix).
For example for the fragment named `H2O_L` the
parameters must be in the ``fraglib_path`` directory in the file named ``h2o.efp``.
For the fragment named `NH3` the parameters must be in the ``userlib_path``
directory in the file named ``nh3.efp``.

See :ref:`FRAGNAME` for additional information. 


.. _fragment input:

Fragment input
^^^^^^^^^^^^^^^

Fragment input should contain one or more groups starting with ``fragment <name>`` line 
(see :ref:`fragment name`) and 
followed by fragment coordinates in one of four possible formats controlled by ``coord`` 
keyword (see :ref:`Basic parameters`).

**``xyzabc`` format**

.. code-block:: default

	fragment h2o
	0.0 0.0 0.0 0.0 0.0 0.0

The numbers are coordinates of the center of mass of a fragment in Angstroms
and three Euler rotation angels in Radians.

**``points`` format**

.. code-block:: default

	fragment h2o
		0.0 0.0 0.0
		1.0 0.0 0.0
		0.0 1.0 0.0

Exactly three lines of three numbers is expected. These are x,y,z coordinates of the first three 
points (i.e., the first three atoms specified in the `.efp` potential) of a fragment in Angstroms.

**``rotmat`` format**

.. code-block:: default

	fragment h2o
		0.0 0.0 0.0
		1.0 0.0 0.0
		0.0 1.0 0.0
		0.0 0.0 1.0

The first line contain x,y,z coordinates of fragment's center of mass in Angstroms.
The next three lines contain a 3 x 3 rotation matrix determined with respect to the 
reference fragment stored in `.efp` file.

**``atoms`` format**

.. code-block:: default

	fragment nh3_l
	A01N1                0.235191     4.772273    -0.103166
	A02H2               -0.523189     4.322511     0.366654
	A03H3                0.517310     5.557741     0.446020
	A04H4                0.999952     4.130518    -0.141193

Atom name (as defined in `.efp`) and x,y,z coordinates in Angstroms for all 
fragment atoms should be given. 

.. _extra input:

Additional fragment input
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Fragment velocities**

Each fragment input can be followed by fragment velocities (for restarting MD simulations)
initiated with the ``velocity`` 
keyword, with the center of mass velocity and angular velocity in atomic units
specified on the next line. An example showing the usage: 

.. literalinclude:: ../examples/md_1.in
  :linenos:

**Fragment constraints**

Quadratic constraint on the fragment center of mass can be specified using
``constraint`` keyword with the force constant ``k`` (in a.u.) and constraint
position ``xyz`` (in Angstroms) specified on the next line. Example:

.. literalinclude:: ../examples/constraint_3.in
  :linenos:

.. _point charges:

Point charges
^^^^^^^^^^^^^^

Additionally to fragments, a system can contain a set of point charges. They
can be specified using the following format for each charge:

.. code-block:: default

	charge <q> <x> <y> <z>




