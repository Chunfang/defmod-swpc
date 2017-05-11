## FE (defmod) installation

Download PETSc from https://www.mcs.anl.gov/petsc/download/index.html and unzip.   
Add export PETSC_DIR=/home/cmeng/petsc-[version] to ~/.bashrc (.profile)  
cd petsc-[version]

For intel-mpi cluster (change the mpi, blas/lapack and scalapack args accordingly for different clusters), Comment out all the intel/mpi modules in the .bashrc file.

export PETSC_ARCH=linux-intel-opt
./configure  --with-mpi-dir=/cm/shared/engaging/intel/intel-2013.1.046/impi/4.1.3.048/intel64 --with-blas-lapack-dir=/cm/shared/engaging/intel/intel-2013.1.046/composerxe/mkl/lib/intel64 --with-scalapack-lib="-L/cm/shared/engaging/intel/intel-2013.1.046/composerxe/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" --download-hdf5 --download-parmetis --download-mumps --download-metis --with-shared-libraries=0 --with-debugging=0 --COPTFLAGS=-O2 --FOPTFLAGS=-O2

Uncomment intel/mpi modules 

For gnu mvapich2 cluster (change the mpi and blas/lapack args accordingly for different clusters)

./configure --with-blas-lapack-dir=/usr/lib64/atlas --with-mpi-dir=/cm/shared/apps/mvapich2/gcc/64/2.0b --download-hdf5 --download-scalapack --download-parmetis --download-mumps --download-metis --with-shared-libraries=0 --with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3

For intel/openmpi Ubuntu OS (Change blas-lapack args accordingly for other Linux distros)

sudo apt-get install cmake

export PETSC_ARCH=linux-gnu-opt
./configure --with-blas-lapack-dir=/usr/lib/atlas-base --with-mpi=1 --with-cc=mpicc --with-fc=mpif90 --with-cxx=mpicxx --download-hdf5 --download-scalapack --download-parmetis --download-mumps --download-metis --download-ptscotch --download-scotch --download-hypre --download-ml --with-shared-libraries=0 --with-debugging=0 COPTFLAGS=-O3 FOPTFLAGS=-O3 CXXOPTFLAGS=-O3

To suppress gnu warnings add to the file src/Makefile FFLAGS = -Wno-maybe-uninitialized

All configurations can have --with-debugging=1 for debug. To inspect the allocatable array/pointer with gdb+gfortran:

(gdb) print *((real *)A+m)@n  
A: array/pointer name  
data type: real_8/int/long_int/  
m: 2 means print from 3rd entry  
n: number of entries to print  

With PETSC_DIR and PETSC_ARCH set correctly in ~/.bashrc
cd src/defmod && make all

## FD (swpc) sintallation

Install NetCDF and NetCDF-fortran, following http://www.unidata.ucar.edu/software/netcdf. If NetCDF is already available, only NetCDF-fortran needs to install. 

export NCDIR=[NetCDF]  
export CC=gcc(icc)  
export FC=gfortran(ifort)  
export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}  
export NFDIR=[NetCDF-fortan]  
export CPPFLAGS=-I${NCDIR}/include LDFLAGS=-L${NCDIR}/lib   
./configure --prefix=${NFDIR} --disable-fortran-type-check   
make check && make install

Add new entries to src/share/makefile.arch and src/share/makefile-tool.arch for local environment. Use existing entries, e.g. arch=mac-gfortran, as template. To set the NetCDF paths correctly.  
NCLIB   = -L[NetCDF-fortan]/lib -L[NetCDF]/lib  
NCINC   = -I[NetCDF-fortan]/include  

cd src/swpc3d  
make arch=[local arch] debug(publish)=true/false  
cd src/tool  
make arch=[local arch] debug(publish)=true/false

## Running mixed code

Note, defmod and swpc can run standalone, compatible with all the model inputs for the original codes,  
https://bitbucket.org/stali/defmod  
https://bitbucket.org/chunfangmeng/defmod-dev  
https://github.com/takuto-maeda/OpenSWPC

To run the code in FE-FD mixed mode, two additional input files are needed, [my_model]_fe.cfg [my_model]_fd.cfg

Find, SCEC[ID]_fe(fd).cfg and F3D_tet_fe(fd).cfg in example/SCEC and example/F3Dp, and follow the comments. Note, the values in [my_model]_fe.cfg should be consistent with those in [my_model].inf for swpc, and the values in [my_model]_fd.cfg, i.e. observation list, should be consistent with those in [my_model].inp for defmod. Use the two examples SCEC and F3Dp (poroelastic) as templates.

First run defmod (FE) model. -fd 1 tells the code to output the FD input.  
[MPI_cmd] [MPI-args] bin/defmod -f [my_fe_model].inp [petsc-args] -fd 1

After FE model finishes without errors, run swpc (FD) model. -r [my_fe_model] passes the code the FE model name, no extension. -e [event_ID] passes which earthquake to simulate if there are more then one event. Without -e the code will pick the first event.   
[MPI_cmd] [MPI-args] bin/swpc_3d.x   -i [my_fd_model].inf -r [my_fe_model] -e [event_ID] 

The code will display if the FE-FD mode is on via screen outputs. 

Use example/SCEC/SCEC_rup(sw).sh and example/F3Dp/F3D_rup(sw).sh as templates for command lines or job submissions. Also, use *.py files in the two examples as templates for making FE models.

To generate FE mesh, in Exodus-II (.exo) format,  
cubit(trelis) -nojournal -nographics [my_model].jou  
To generate FE inp file,  
[my_model_gen].py [my_model].exo
