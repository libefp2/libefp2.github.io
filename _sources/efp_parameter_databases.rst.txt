.. _efp_parameter_databases:

***********************
EFP parameter databases
***********************

S22 Parameters:
^^^^^^^^^^^^^^^
Geometries from `doi.org/10.1039/B600027D <http://dx.doi.org/10.1039/b600027d>`_

- .. raw:: html

   <a href="examples/parameters/2aminopyridine.efp" target="_blank">2aminopyridine.efp</a>
- .. raw:: html

   <a href="examples/parameters/adenine-stack.efp" target="_blank">adenine-stack.efp</a>
- .. raw:: html

   <a href="examples/parameters/adenine-wc.efp" target="_blank">adenine-wc.efp</a>
- .. raw:: html

   <a href="examples/parameters/ethene.efp" target="_blank">ethene.efp</a>
- .. raw:: html

   <a href="examples/parameters/ethyne.efp" target="_blank">ethyne.efp</a>
- .. raw:: html

   <a href="examples/parameters/formicacid.efp" target="_blank">formicacid.efp</a>
- .. raw:: html

   <a href="examples/parameters/hydrogencyanide.efp" target="_blank">hydrogencyanide.efp</a>
- .. raw:: html

   <a href="examples/parameters/indole.efp" target="_blank">indole.efp</a>
- .. raw:: html

   <a href="examples/parameters/pyrazine.efp" target="_blank">pyrazine.efp</a>
- .. raw:: html

   <a href="examples/parameters/pyridone.efp" target="_blank">pyridone.efp</a>
- .. raw:: html

   <a href="examples/parameters/thymine-stack.efp" target="_blank">thymine-stack.efp</a>
- .. raw:: html

   <a href="examples/parameters/thymine-wc.efp" target="_blank">thymine-wc.efp</a>
- .. raw:: html

   <a href="examples/parameters/uracil.efp" target="_blank">uracil.efp</a>

S66 Parameters:
^^^^^^^^^^^^^^^
Geometries from `dx.doi.org/10.1021/ct200673a <http://dx.doi.org/10.1021/ct200673a>`_

- .. raw:: html

   <a href="examples/parameters/acetamide-gp.efp" target="_blank">acetamide-gp.efp</a>
- .. raw:: html

   <a href="examples/parameters/acetamide-hb.efp" target="_blank">acetamide-hb.efp</a>
- .. raw:: html

   <a href="examples/parameters/aceticacid-gp.efp" target="_blank">aceticacid-gp.efp</a>
- .. raw:: html

   <a href="examples/parameters/aceticacid-hb.efp" target="_blank">aceticacid-hb.efp</a>
- .. raw:: html

   <a href="examples/parameters/ammonia.efp" target="_blank">ammonia.efp</a>
- .. raw:: html

   <a href="examples/parameters/benzene.efp" target="_blank">benzene.efp</a>
- .. raw:: html

   <a href="examples/parameters/cyclopentane.efp" target="_blank">cyclopentane.efp</a>
- .. raw:: html

   <a href="examples/parameters/formamide.efp" target="_blank">formamide.efp</a>
- .. raw:: html

   <a href="examples/parameters/methane.efp" target="_blank">methane.efp</a>
- .. raw:: html

   <a href="examples/parameters/methanol.efp" target="_blank">methanol.efp</a>
- .. raw:: html

   <a href="examples/parameters/methylamine.efp" target="_blank">methylamine.efp</a>
- .. raw:: html

   <a href="examples/parameters/neopentane.efp" target="_blank">neopentane.efp</a>
- .. raw:: html

   <a href="examples/parameters/nmethylacetamide.efp" target="_blank">nmethylacetamide.efp</a>
- .. raw:: html

   <a href="examples/parameters/pentane.efp" target="_blank">pentane.efp</a>
- .. raw:: html

   <a href="examples/parameters/phenol.efp" target="_blank">phenol.efp</a>
- .. raw:: html

   <a href="examples/parameters/pyridine.efp" target="_blank">pyridine.efp</a>
- .. raw:: html

   <a href="examples/parameters/uracil-gp.efp" target="_blank">uracil-gp.efp</a>
- .. raw:: html

   <a href="examples/parameters/water.efp" target="_blank">water.efp</a>

Amino acid Parameters:
^^^^^^^^^^^^^^^^^^^^^^
- .. raw:: html

   Parameters in <a href="examples/parameters/1.6-31g-d" target="_blank">6-31G(d)</a> basis

- .. raw:: html

   Parameters in a hybrid <a href="examples/parameters/6.6-31g-d_6-31+xg-3df-2p" target="_blank">6-31G(d)/6-31-+G(3df,2p)</a> basis

- .. raw:: html

   Parameters in a hybrid <a href="examples/parameters/7.6-31g-d_6-311++g-3df-2p" target="_blank">6-31G(d)/6-311++G(3df,2p)</a> basis

The fragment parameters in FlexibleEFP are adjusted based on translations and rotations of local coordinate frames associated with fragment atoms to accommodate different fragment geometries. A parameter database for standard amino acids was developed to automate flexible EFP simulations in proteins using cryptochrome 1 protein (Arabidopsis thaliana's, Cry1At, PDB: 1U3D64). The applicability of flexible EFP was demonstrated in large-scale protein simulations, where binding energies, as well as vertical electron ionization and electron attachment energies of a lumiflavin chromophore in the cryptochrome 1 protein were computed. The results showed that flexible EFP closely agrees with the standard EFP procedure but with a significant reduction in computational cost. Twenty-five protein conformations from molecular dynamics trajectories of cryptochrome 1 were fragmented into individual amino acid (AA) fragments along Cα−C bonds, yielding 12,125 BioEFP AA fragments. To improve accuracy, each AA fragment was further split along the Cα−Cβ bond into a backbone group and a side-chain group, and disulfide bridges were fragmented along the S−S bond. To recombine the backbone and side-chain fragments into a complete AA moiety, parameters at a bond midpoint and the LMO centroid between Cα and Cβ were excluded to ensure the stability of the polarization self-consistent procedure. For further details please refer to:

- Database of fragments describing amino acid residues and peptide backbone groups in different geometries `Flexible EFP paper <https://doi.org/10.1021/acs.jctc.0c00758>`_.
- Strategy of rotating and shifting the parameters <https://doi.org/10.1021/acs.jctc.0c00758>`_.
- Script that will do it for you `parameter rotation script <https://github.com/libefp2/libefp/blob/master/tools/Flexible_V5.py>`_.
- Splitting of the protein into amino acid fragments and matching the parameters to specific geometries of the fragments with BioEFP and FlexEFP tutorials is described in :ref:`bioefp`.
- Note: Parameters from the database should be matched to the geometry of your system.
