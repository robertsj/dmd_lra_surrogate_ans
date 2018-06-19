# Data

Analysis scripts, intermediate data (or directions to get it), and
final output scripts should be included here.  Snapshots are generated
using Detran.


## Detran

Detran as used to generate the results here is 

    commit 453153856760a31991cff35bb200d8dba252c712

Example Detran and PETSc (version 3.4.5) configuration scripts are given.

Note, these results were generated with Python 2.7 because Detran has not yet 
been extended for Python 3 support.  Further, SWIG version 2.0.12 was used, 
as SWIG 3+ variants appear not to work.

Other Python requirements (with versions used):

```
numpy                     1.14.2          py27_blas_openblas_200  [blas_openblas]  conda-forge
scipy                     1.0.0           py27_blas_openblas_201  [blas_openblas]  conda-forge
```

If the requirements are met, production of the raw data used for 
the DMD analyses can be generated via

```
  $ python 
```
