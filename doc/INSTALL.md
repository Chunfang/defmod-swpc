## FE (defmod) installation

Download PETSc from https://www.mcs.anl.gov/petsc/download/index.html and unzip.   
Add "export PETSC_DIR=${HOME}/petsc-[version]" to ~/.bashrc (.profile)  
cd petsc-[version]

On CentOS cluster mkl+impi+icc, set the mpi, blas/lapack and scalapack arguments according to different Parallel Studio(R) installations. On Engaging cluster (MIT) for example:

module purge   
unset LD_LIBRARY_PATH LIBRARY_PATH  
export PATH=/bin:/usr/bin  
module add engaging/intel/2013.1.046   
unset F77  
export PETSC_DIR=$PWD  
export PETSC_ARCH=mkl-impi-icc  
./configure --with-cc=mpiicc --with-cxx=mpiicpc --with-fc=mpiifort --with-blas-lapack-dir=/cm/shared/engaging/intel/intel-2013.1.046/composerxe/mkl/lib/intel64 --with-scalapack-lib="-L/cm/shared/engaging/intel/intel-2013.1.046/composerxe/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64" --download-metis --download-parmetis --download-hdf5 --download-mumps --with-shared-libraries=0 --with-debugging=0 --COPTFLAGS="-O2 -xHost" --FOPTFLAGS="-O2 -xHost" --CXXOPTFLAGS="-O2 -xHost"  
make all  
(Note, the module/path unset is to prevent customized python and any unnecessary library linkage.)

For atlas+openmpi+gcc option: 

module purge  
unset LD_LIBRARY_PATH LIBRARY_PATH  
export PATH=/bin:/usr/bin  
module add engaging/openmpi/2.0.3   
module add gcc  
unset F77 FC CC CXX  
export PETSC_DIR=$PWD  
export PETSC_ARCH=atlas-ompi-gcc  
./configure  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpifort --with-x=0 --with-blas-lapack-dir=/usr/lib64/atlas  --download-metis --download-parmetis --download-scalapack --download-mumps --download-hdf5 --with-shared-libraries=0 --with-debugging=0 --COPTFLAGS=-O2 --FOPTFLAGS=-O2 --CXXOPTFLAGS=-O2  
make all

For Ubuntu OS mkl+gcc+openmpi, change blas-lapack arguments according to different Linux distros:

sudo apt install cmake  
export PETSC_ARCH=linux-gnu-opt  
./configure --with-blas-lapack-dir=/usr/lib/atlas-base --with-mpi=1 --with-cc=mpicc --with-fc=mpif90 --with-cxx=mpicxx --download-hdf5 --download-scalapack --download-parmetis --download-mumps --download-metis --download-ptscotch --download-scotch --download-hypre --download-ml --with-shared-libraries=0 --with-debugging=0 COPTFLAGS=-O3 FOPTFLAGS=-O3 CXXOPTFLAGS=-O3  
make all

All configurations can have --with-debugging=1 to turn on debug. To inspect the allocatable array/pointer with gdb+gfortran   
(gdb) print *((real *)A+m)@n  
A: array/pointer name  
data type: real_8/int/long_int/  
m: 2 means print from 3rd entry  
n: number of entries to print  

Add "PETSC_ARCH=[math-mpi-compiler]" to ~/.bashrc, source ~/.bashrc, cd [defmod-swpc]/src/defmod 

To suppress false compiler warnings, have "FFLAGS = -Wno-maybe-uninitialized" in Makefile. 

make all

## FD (swpc) installation

Install NetCDF and NetCDF-fortran, following http://www.unidata.ucar.edu/software/netcdf. If NetCDF is already available, only NetCDF-fortran needs to install. 

export NCDIR=[NetCDF]  
export CC=gcc(icc)  
export FC=gfortran(ifort)  
export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}  
export NFDIR=[NetCDF-fortan]  
export CPPFLAGS=-I${NCDIR}/include LDFLAGS=-L${NCDIR}/lib   
./configure --prefix=${NFDIR} --disable-fortran-type-check   
make check && make install

Add new entries to src/share/makefile.arch and src/share/makefile-tool.arch for local the environment. Use existing entries, e.g. arch=mac-gfortran, as template. To set the NetCDF paths correctly,  
NCLIB   = -L[NetCDF-fortan]/lib -L[NetCDF]/lib  
NCINC   = -I[NetCDF-fortan]/include  

cd src/swpc3d  
make arch=[local arch] debug(publish)=true/false  
cd src/tool  
make arch=[local arch] 

## Running the mixed code

defmod and swpc can run standalone, compatible with all the model inputs for the original codes:  
https://bitbucket.org/stali/defmod  
https://github.com/takuto-maeda/OpenSWPC

To run the code in FE-FD coupled mode, an additional file "[xxx]_fefd.cfg" is needed, use the files in example/SCEC as templates. Note, the MPI size in [xxx]_fefd.cfg should be consistent with those in [xxx].inf for the FD part, and the observation grid should have overlap with those in [xxx].inp for cross validation.

First run defmod (FE) model. -fd 1 tells the code to output the FD input, and -ss 0/1 is to switch off/on snapshot (VTK) output:  
[MPI_cmd] [MPI-args] [defmod-swpc]/bin/defmod -f [xxx].inp [petsc-args] -ss 0/1 -fd 1

After FE model finishes without errors, run swpc (FD) model. -r [xxx] (no extension) passes the FE model name. -e [event number] tells which event to simulate in case of multiple events. Without -e argument, swpc will pick the first event:   
[MPI_cmd] [MPI-args] [defmod-swpc]/bin/swpc_3d.x -i [fd model].inf -r [xxx] -e [event number] 

Other files, e.g. station list and velocity model, are the same for original OpenSWPC. The FE velocity will override the FE velocity, be careful about velocity consistency at the interface. Both the codes display the FE-FD binding information via screen outputs. 

To run models on cluster, use the submission files [xxx]_rup(sw).sh in example/[xxx] as template. Also, use the [xxx].py files as template for making the FE model inputs, and makeing plots/movies.

To generate FE mesh, in Exodus-II (.exo) format:  
cubit(trelis) -nojournal -nographics [xxx].jou  
To generate FE inp file:  
./[xxx prep].py [xxx].exo [options]
