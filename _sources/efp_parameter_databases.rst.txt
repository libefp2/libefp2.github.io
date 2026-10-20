.. _efp_parameter_databases:

***********************
EFP parameter databases
***********************

S22 and S66 Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The accuracy of EFP was benchmarked against SAPT  on the example of dimers 
from S22 and S66 datasets of non-covalent interactions (see `S22 EFP paper <http://dx.doi.org/10.1039/b600027d>`_).
Developed EFP parameters can be found here `EFP S22 S66 parameters <https://github.com/libefp2/EFP_parameters/tree/main/S22_S66>`_.

.. note:: S22/S66 EFP parameters are hybrid basis parameters (6-31G*/6-31+G* for electrostatics and 
   6-311++G(3df,2p) for polarization, dispersion and exchange-repulsion terms)

For convinience, these parameters are also provided below.

**S22 set** (geometries given in `S22 EFP paper <http://dx.doi.org/10.1039/b600027d>`_)

- `2aminopyridine.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/2aminopyridine.efp>`_
- `adenine-stack.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/adenine-stack.efp>`_
- `adenine-wc.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/adenine-wc.efp>`_
- `ethene.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/ethene.efp>`_
- `ethyne.efp.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/ethyne.efp.efp>`_
- `formicacid.efp.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/formicacid.efp.efp>`_
- `hydrogencyanide.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/hydrogencyanide.efp>`_
- `indole.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/indole.efp>`_
- `pyrazine.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/pyrazine.efp>`_
- `pyridone.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/pyridone.efp>`_
- `thymine-stack.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/thymine-stack.efp>`_
- `thymine-wc.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/thymine-wc.efp>`_
- `uracil.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/uracil.efp>`_

**S66 set** (original S66 geometries from `here <http://dx.doi.org/10.1021/ct200673a>`_)

- `acetamide-gp.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/acetamide-gp.efp>`_
- `acetamide-hb.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/acetamide-hb.efp>`_
- `aceticacid-gp.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/aceticacid-gp.efp>`_
- `aceticacid-hb.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/aceticacid-hb.efp>`_
- `ammonia.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/ammonia.efp>`_
- `benzene.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/benzene.efp>`_
- `cyclopentane.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/cyclopentane.efp>`_
- `formamide.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/formamide.efp>`_
- `methane.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/methane.efp>`_
- `methanol.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/methanol.efp>`_
- `methylamine.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/methylamine.efp>`_
- `neopentane.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/neopentane.efp>`_
- `nmethylacetamide.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/nmethylacetamide.efp>`_
- `pentane.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/pentane.efp>`_
- `phenol.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/phenol.efp>`_
- `pyridine.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/pyridine.efp>`_
- `uracil-gp.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/uracil-gp.efp>`_
- `water.efp <https://github.com/libefp2/EFP_parameters/blob/main/S22_S66/water.efp>`_


Amino acid parameters
^^^^^^^^^^^^^^^^^^^^^^

Parameters for amino acids were created as described in `Flexible EFP paper <https://doi.org/10.1021/acs.jctc.0c00758>`_.
The initial structures of amino acids were obtained from MD snapshots of cryptochrome 1 
protein (Arabidopsis thaliana's, Cry1At, PDB: 1U3D64), and futher geometry optimized. EFP parameters were 
computed at all unique minima. Three sets of parameters are available:

- single-basis 6-31G(d) set
- hybrid 6-31G(d)/6-31-+G(3df,2p) basis (diffuse functions are added to H atoms only)
- hybrid 6-31G(d)/6-311++G(3df,2p) basis (standard hybrid basis)

These EFP parameters can be found here: `AA flexible EFP parameters <https://github.com/libefp2/EFP_parameters/tree/main/AA_flexible_efp>`_.

This database is assumed to be used in combinatin with BioEFP and flexible EFP schemes, as described 
in :ref:`bioefp` and :ref:`flexible_efp`, respectively. 

.. note:: We have not used the quality of this database for other proteins; please use it 
   cautiosly.  
