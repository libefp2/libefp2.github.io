.. _gamess:

************
GAMESS
************

GAMESS official website: `GAMESSS <https://www.msg.chem.iastate.edu/gamess/index.html>`_

`GAMESS documentation <https://www.msg.chem.iastate.edu/gamess/documentation.html>`_

GAMESS contains original and the most complete implementation of the EFP method, as well as 
functionality to compute EFP parameter files (see :ref:`compute parameters`).

For EFP calculations in GAMESS, check the following groups in the GAMESS manual:

- generation of EFP potential
    - $MAKEFP - main group controlling generation of effective fragment potential
    - $DAMP  -  controls generation of the short-range electrosatic screening parameters  
    - $DAMPGS - provides initial guess for generation of the short-range screening parameters    

- performing EFP and QM/EFP calculations
    - $EFRAG - the main group controlling EFP calculation  
    - $FRAGNAME - specifically named fragment groups (e.g., $WATER, $CH2O, etc.), containing fragment potentials
    - $QMEFP -  QM/EFP energy decomposition       



