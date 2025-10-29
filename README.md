This code is used for the coupling simulation of LaMEM and FastScape.

To install this code, PETSc 3.18.6 and FastScapeLib is needed. The installation can be referred to: 

[LaMEM](https://unimainzgeo.github.io/LaMEM/dev/man/Installation/)
[FastScapeLib](https://fastscape.org/fastscapelib-fortran/)

Next you need to specify the environmental variables PETSC_OPT and FASTSCAPE_LIB
export PETSC_OPT=/path/you/install/
export FASTSCAPE_LIB=/path/you/install/

Next you can install the coupling module through:
make mode=optFS all surface=SURFACE=1

