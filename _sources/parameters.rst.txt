.. _parameters:

**************
EFP parameters
**************

How to prepare EFP parameters? There are two ways, either find existing parameters in a :ref:`parameter database` or :ref:`compute parameters`.

.. _parameter database:

Parameter database
------------------

We have prepared and validated EFP parameters for the following systems:

- `monomers from S22 and S66 sets <https://github.com/libefp2/libefp/tree/master/fraglib/databases>`_
- `amino acid parameters from doi.org/10.1021/acs.jctc.0c00758 <https://github.com/libefp2/libefp/tree/master/fraglib/flex_efp_protein>`_

Learn more detail about these EFP parameters on :ref:`efp_parameter_databases`.


.. _compute parameters:

Computing parameters
--------------------

Parameters for any molecule can be computed using MAKEFP run in the GAMESS electronic structure package.

For general information on using GAMESS please refer to the
`GAMESS manual <https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/docs-input.txt>`_.

Apart of standard `.log` and `.dat` outputs, MAKEFP run creates  an additional file with an extension `.efp`
that contains all EFP parameters of your molecule.

An example of GAMESS input file for computing EFP parameters :download:`makefp.inp <../examples/makefp.inp>`.

.. literalinclude:: ../examples/makefp.inp
  :linenos:

.. _efp parameters:

`.efp` parameter file
----------------------

.. note:: `.efp` file follows general GAMESS input file conventions:

   - Group names starting with ``$`` should start at the second line character.
   - Lines should not exceed 72 characters; use ``>`` as a continuation sign for longer lines.

An example of `.efp` file :download:`nh3.efp <../examples/nh3.efp>`.

The `.efp` file has the following structure (see also `$FRAGNAME group` section in GAMESS manual):

.. _FRAGNAME:

``$FRAGNAME``
^^^^^^^^^^^^^

- ``FRAGNAME`` section provides user-specified fragment name.
- The section is *mandatory*.
- ``FRAGNAME`` tag should match the fragment name in LibEFP, GAMESS, Q-Chem inputs.
  Capitalization is not important.
- ``FRAGNAME`` tag should match the name of the `.efp` file.
  For example, ``$nh3`` fragment should be contained in `nh3.efp` file.
- ``$`` in front of the ``FRAGNAME`` tag should be positioned at the second line character.


``Commment line``
^^^^^^^^^^^^^^^^^

- The section is *mandatory*.

.. _COORD:

``COORDINATES``
^^^^^^^^^^^^^^^

- ``COORDINATES`` section provides coordinates of all atoms and bond-midpoints of the fragment.
- The section is *mandatory*.
- Each line contains name tag, x, y, z coordinates (in Bohr), atom mass, nuclear charge.
- For atoms, the name tag is produced as ``A`` + ``two-digit number`` + ``atom name`` from the MAKEFP input,
  e.g., A01N.
- For bond-midpoints, the name tag is produced as ``BO`` + ``second atom number`` + ``first atom number``,
  e.g., BO21.
- These name tags are used in ``MULTIPOLES``, ``DIPOLES``, ``QUADRUPOLES``, ``OCTUPOLES``, ``PROJECTION BASIS SET``, and
  ``SCREEN`` and ``SCREEN2`` sections.
- The name tags can be changed by user if needed (typically not needed, see below).
- GAMESS does not correctly assign name tags in fragments containing more than 99 atoms.
  Such potentials require atom and bond-midpoint renaming by user.
- The section ends with mandatory ``STOP`` line.

.. _MON:

``MONOPOLES``
^^^^^^^^^^^^^

- ``MONOPOLES`` section contains nuclear and electronic charges of atoms and bond-midpoints of the fragment.
- Each line contains name tag, electronic charge, nuclear charge.
- Name tags should match those in ``COORDINATES`` section.
- The section is mandatory in GAMESS input and optional in LibEFP and Q-Chem inputs.
- The section can contain less atoms/bond midpoints than ``COORDINATES`` section.
- The section ends with mandatory ``STOP`` line.

.. _DIP:

``DIPOLES``
^^^^^^^^^^^

- ``DIPOLES`` section contains point dipoles of atoms and bond-midpoints of the fragment.
- Each line contains name tag and x, y, z values of a point dipole in atomic units.
- Name tags should match those in ``COORDINATES`` section.
- The section is optional.
- The section can contain less atoms/bond midpoints than ``COORDINATES`` section.
- Direction of the dipoles is dependent on the fragment orientation as provided by the fragment coordinates.
- The section ends with mandatory ``STOP`` line.

.. _QUAD:

``QUADRUPOLES``
^^^^^^^^^^^^^^^

- ``QUADRUPOLES`` section contains point quadrupoles of atoms and bond-midpoints of the fragment.
- Each line contains name tag and xx, yy, zz, xy, xz, yz components of a quadrupole
  moment tensor in atomic units. Three other elements (yx, zx, zy) are identical
  to xy, xz, yz respectively and are omitted.
- Name tags should match those in ``COORDINATES`` section.
- The section is optional.
- The section can contain less atoms/bond midpoints than ``COORDINATES`` section.
- Direction of the quadrupole tensor is dependent on the fragment orientation as provided by the fragment coordinates.
- The quadrupole moments are given in the Buckingham form, i.e.,  xx, xy, xz etc components correspond to
  :math:`a_x a_x`, :math:`a_x a_y`, :math:`a_x a_z`, respectively,
  where **a** is a vector describing particle position. Conversion to multipoles in the Cartesian form
  :math:`\Theta_{\alpha \beta} = \sum_a e_a (\frac{3}{2}a_{\alpha}a_{\beta}-\frac{1}{2}a^2\delta_{\alpha \beta})`
  is conducted within LibEFP or GAMESS codes.
- The section ends with mandatory ``STOP`` line.

.. _OCT:

``OCTUPOLES``
^^^^^^^^^^^^^

- ``OCTUPOLES`` section contains point octupoles of atoms and bond-midpoints of the fragment.
- Each line contains name tag and xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz components of an octupole
  moment tensor in atomic units. Other elements are determined by symmetry and are omitted.
- Name tags should match those in ``COORDINATES`` section.
- The section is optional.
- The section can contain less atoms/bond midpoints than ``COORDINATES`` section.
- Direction of the quadrupole tensor is dependent on the fragment orientation as provided by the fragment coordinates.
- The octupole moments are given in the Buckingham form, e.g., xxx corresponds to
  :math:`a_x a_x a_x`,
  where **a** is a vector describing particle position. Conversion to the Cartesian form
  is conducted within LibEFP or GAMESS codes.
- The section ends with mandatory ``STOP`` line.

In the example below, ``$FRAGNAME`` (``$nh3`` in this case), ``COORDINATES``, ``MONOPOLES``,
``DIPOLES``, ``QUADRUPOLES``, ``OCTUPOLES`` sections are highlighted.

.. literalinclude:: ../examples/nh3.efp
  :linenos:
  :lines: 1-70
  :emphasize-lines: 3, 5, 14, 23, 32, 48

.. _SCREEN:

``SCREEN`` and ``SCREEN2``
^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``SCREEN`` and ``SCREEN2`` sections provide screening parameters
  to account for charge-penetration energy between QM and EFP regions (``SCREEN``) and
  between EFP fragments (``SCREEN2``).
- Both sections are optional.
- Each line contains name tag and :math:`\alpha` and :math:`\beta` parameters that
  effectively screen EFP electronic charges.

  For QM-EFP interactions (``SCREEN``):
  :math:`q_A \rightarrow _A(1-\alpha \exp(-\beta (r_A-r_e)^2))`,
  where :math:`r_e` is position of an electron.

  For EFP-EFP interactions (``SCREEN2``):
  :math:`q_A \rightarrow _A(1-\alpha \exp(-\beta (r_A-r_B)))`,
  where :math:`r_B` is position of another charge.

  :math:`\alpha` parameters in both cases are typically set to 1.

- ``SCREEN`` and ``SCREEN2`` sections are typically found in the very end of `.efp` file,
  even though exact placement of sections can be varied by user.
- Name tags should match those in ``COORDINATES`` section.
- Sections can contain less atoms/bond midpoints than ``COORDINATES`` section.
- Each section ends with mandatory ``STOP`` line.

.. literalinclude:: ../examples/nh3.efp
  :linenos:
  :lines: 361-
  :emphasize-lines: 1, 10

.. _POL_POINT:

``POLARIZABLE POINTS``
^^^^^^^^^^^^^^^^^^^^^^

- ``POLARIZABLE POINTS`` section provides coordinates and values of *anisotropic*
  polarizability tensors for computing polarization energy.
- The section is optional.
- The group contains a set of two-line entries each containing:

  * first line: polarizable point tag, x, y, z coordinates of the point (in Bohr)
  * second line: XX, YY, ZZ, XY, XZ, YZ, YX, ZX, ZY elements of the polarizability tensor

- The section ends with mandatory ``STOP`` line.
- Coordinates of the polarizability points and orientations of polarizability tensors
  are given with respect to the fragment coordinate frame as given in :ref:`COORD`.
- Typically coordinates of polarizability points are centroids of the localized molecular
  orbitals, but other choices (e.g., atom positions) are possible.

.. _DYN_POINT:

``DYNAMIC POLARIZABLE POINTS``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``DYNAMIC POLARIZABLE POINTS`` section provides coordinates and values of *anisotropic*
  dynamic polarizability tensors for computing dispersion energy.
- The section is optional.
- The group contains twelve entries each similar to the :ref:`POL_POINT` section, i.e.,
  each containing multiple sets of two-line entries:

  * first line: polarizable point tag, x, y, z coordinates of the point (in Bohr)
  * second line: XX, YY, ZZ, XY, XZ, YZ, YX, ZX, ZY elements of the dynamic polarizability tensor

  Each of these 12 sets contains time-dependent (dynamic) polarizabilities computed at given
  values of imaginary frequencies. Coordinates of polarizability points remain the same
  between all twelve sets.

- The section ends with mandatory ``STOP`` line.
- Coordinates of the dynamic polarizability points and orientations of tensors
  are given with respect to the fragment coordinate frame as given in :ref:`COORD`.
- Typically coordinates of dynamic polarizability points are centroids of the localized molecular
  orbitals, but other choices (e.g., atom positions) are possible.
- Only first three elements (XX, YY, ZZ) of the tensors are used in calculation of
  the :math:`\frac{C_6}{R^6}` dispersion energy (the one implemented in LibEFP).

.. literalinclude:: ../examples/nh3.efp
  :linenos:
  :lines: 71-137
  :emphasize-lines: 1, 19, 20, 36, 52

.. _BASIS:

``PROJECTION BASIS SET``
^^^^^^^^^^^^^^^^^^^^^^^^

- ``PROJECTION BASIS SET`` section contains information on the basis set used in
  calculation of the exchange-repulsion energy and overlap-based electrostatic and
  dispersion dampings.
- The section is optional.
- Basis set information for all atoms of the fragment is given, separated by empty lines.
- The first line for each atom entry contains atom name tag that should match one from
  the :ref:`COORD` section, x, y, z coordinates of the atom (in Bohr), and a number corresponding to
  the atom nuclear charge minus the number of electrons in core orbitals. For example,
  5 = 7-2 for N, 4 = 6-2 for C, 1 = 1-0 for H etc. The following information on basis functions
  is given in the GAMESS basis set format.
- The section ends with mandatory ``STOP`` line.

.. literalinclude:: ../examples/nh3.efp
  :linenos:
  :lines: 283-327

.. _MULTIPLICITY:

``MULTIPLICITY``
^^^^^^^^^^^^^^^^

- ``MULTIPLICITY`` section provides multiplicity of the fragment used in computing
  exchange-repulsion energy and overlap-based electrostatic and
  dispersion dampings.
- The section is optional.
- LibEFP works only with multiplicity 1 (closed shell fragments).
- The section ends with mandatory ``STOP`` line.

.. _WF:

``PROJECTION WAVEFUNCTION``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``PROJECTION WAVEFUNCTION`` section provides localized wave function of the fragment used in computing
  exchange-repulsion energy and overlap-based electrostatic and
  dispersion dampings.
- The section is optional.
- Format: first line contains PROJECTION WAVEFUNCTION keyword, :math:`N_{loc}`-  number of the localized occupied orbitals
  (core orbitals are excluded), :math:`N_{basis}` size of the basis set. The following lines provide orbital coefficients for
  each orbital. Each line contains orbital number, line number for the given orbital,
  five orbital coefficients not separated by comma. These entries are format sensitive!
- The number of orbital entries should match :math:`N_{loc}` (the first number in the first line).
- The number of orbital coefficients for each orbital should match :math:`N_{basis}`
  (the second number in the first line).
- Note: the section does not require ``STOP`` line.

.. _FOCK:

``FOCK MATRIX ELEMENTS``
^^^^^^^^^^^^^^^^^^^^^^^^

- ``FOCK MATRIX ELEMENTS`` section provides elements of the Fock matrix of the fragment
  in the localized basis. It is used in computing exchange-repulsion energy and overlap-based electrostatic and
  dispersion dampings.
- The section is optional.
- The Fock matrix is given in the lower triangular form, i.e., the matrix
  :math:`\begin{pmatrix}
  a_{11} & a_{12} & a_{13}\\
  a_{21} & a_{22} & a_{23}\\
  a_{31} & a_{32} & a_{33}
  \end{pmatrix}`
  is provided as :math:`a_{11}`, :math:`a_{21}`, :math:`a_{22}`,
  :math:`a_{31}`, :math:`a_{32}`, :math:`a_{33}`.
- The size of the matrix is :math:`N_{loc} \times N_{loc}`, i.e. :math:`N_{loc} * (N_{loc}-1)`
  values is expected.
- Note: the section does not require ``STOP`` line.

.. _LMOC:

``LMO CENTROIDS``
^^^^^^^^^^^^^^^^^^^^^^^^

- ``LMO CENTROIDS`` section provides coordinates of the localized molecular orbital
  centroids. Those are used in computing exchange-repulsion energy and overlap-based electrostatic and
  dispersion dampings.
- The section is optional.
- Each line contains a centroid tag name and x, y, z coordinates (in Bohr) of the centroid, with
  respect to the fragment reference frame as given in :ref:`COORD` section.
- The number of centroids should match :math:`N_{loc}` (the number of the localized occupied
  orbitals) in :ref:`WF` section.
- The section ends with mandatory ``STOP`` line.

.. literalinclude:: ../examples/nh3.efp
  :linenos:
  :lines: 328-360
  :emphasize-lines: 1, 3, 24, 28

.. _additional efp keywords

Additional `.efp` parameters
-----------------------------

.. _POLAB:

``POLAB``
^^^^^^^^^

- ``POLAB`` is a keyword specifying a short-range polarization damping parameter for the fragment.
- The section is optional.
- This section is not included in the `.efp` as prepared by the MAKEFP run and should be added manually if needed.
- The section format is:
  ``POLAB value``, followed by the mandatory ``STOP`` line,
  where ``value`` is the user-specified value of the parameter.
- Default value is 0.6.

Some wisdom on a choice of the POLAB parameter:

- polab = 0.6 (default) - works well for the majority of neutral solvents;
- polab = 0.3 might improve stability of BioEFP calculations in proteins, especially relevant for charged amino acids;
- polab = 0.1 was shown to work safely for anions like F-, Cl-, etc;
- polab = 0.1 might be useful for large aromatics like porphyrins and chlorophylls
  (their polarizabilities are huge and orbital localization schemes often struggle, resulting in very large ~100 a.u.
  polarizabilities on some sites);
- polab = 0.03 for Li+, Na+ and Mg2+.

.. _MM CHARGE:

``MM_CHARGE``
^^^^^^^^^^^^^^

- ``MM_CHARGE`` provides atomic charges, typically matching partial charges in classical MM calculations.
- Each line contains name tag and partial charge.
- Name tags should match those in ``COORDINATES`` section.
- The section is optional.
- This section is not included in the `.efp` as prepared by the MAKEFP run and should be added manually if needed.
- The section ends with mandatory ``STOP`` line.

.. _MM LJ:

``MM_LJ``
^^^^^^^^^^^^^^

- ``MM_LJ`` provides atomic Lennard-Jones parameters as in classical MM calculations.
- Each line contains name tag, :math:`\sigma`, and :math:`\epsilon`.
- Name tags should match those in ``COORDINATES`` section.
- The section is optional.
- This section is not included in the `.efp` as prepared by the MAKEFP run and should be added manually if needed.
- The section ends with mandatory ``STOP`` line.
- GROMACS units are assumed: :math:`\sigma` in nm, :math:`\epsilon` in kJ/mol.

Currently, the Lennard-Jones energy is computed as

.. math::

  E_{ij} = 4*\epsilon_{ij}*( (\frac{\sigma_{ij}}{r})^{12} - (\frac{\sigma_{ij}}{r})^6 )

using the following comboination rules

.. math::

  \sigma_{ij} = 0.5 * (\sigma_i + \sigma_j)

  \epsilon_{ij} = \sqrt{\epsilon_i * \epsilon_j}

.. literalinclude:: ../examples/tip3p_mm.efp
  :linenos:
  :lines: 73-82
  :emphasize-lines: 1, 6
