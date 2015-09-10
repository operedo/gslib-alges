

MPIFC=mpif90
FC=gfortran
FFLAGS= -cpp -O3 -DUNFORMATTED 
#FFLAGS= -cpp -O3
OPENMP=-fopenmp 
#FLIBS= -L$(EXTRAE_HOME)/lib -lomptrace 
GSLIB_HOME=../../gslib
GSLIB=$(GSLIB_HOME)/gslib.a


all: seq openmp_par openmp_mpi_par

seq: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd -
	$(FC) $(FFLAGS) -c kt3d.fpp
	$(FC) $(FFLAGS) -I. kt3d.o -o kt3d_Seq $(GSLIB)
	cp kt3d_Seq ../bin

openmp_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(FC) $(FFLAGS) $(OPENMP) -c kt3d.fpp
	$(FC) $(FFLAGS) $(OPENMP) -I. kt3d.o -o kt3d_OpenMP $(GSLIB)
	cp kt3d_OpenMP ../bin

openmp_mpi_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(MPIFC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -c kt3d.fpp
	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -I. kt3d.o -o kt3d_OpenMP_MPI $(GSLIB)
	cp kt3d_OpenMP_MPI ../bin 

#ins:
#	cd $(HOME)/gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -c extrae_module.f90
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -c kt3d.fpp
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE kt3d.o -o kt3d.exe $(GSLIB) $(FLIBS) -lgomp  

clean:
	rm kt3d_Seq kt3d_OpenMP kt3d_OpenMP_MPI ../bin/* *.o *.mod
