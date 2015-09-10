MPIFC=mpif90
FC=gfortran
#FFLAGS= -cpp -O3
FFLAGS= -cpp -O3 -DDEBUG
#FFLAGS= -cpp -O3 -DUNFORMATTED 
#FFLAGS= -cpp -O3 -DUNFORMATTED -DBUFFERED
OPENMP=-fopenmp 
#FLIBS= -L$(EXTRAE_HOME)/lib -lomptrace 
GSLIB_HOME=../../gslib
GSLIB=$(GSLIB_HOME)/gslib.a


all: seq openmp_par openmp_mpi_par

seq: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd -
	$(FC) $(FFLAGS) -c sisim.fpp
	$(FC) $(FFLAGS) -I. sisim.o -o sisim_Seq $(GSLIB)
	cp sisim_Seq ../bin 

openmp_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(FC) $(FFLAGS) $(OPENMP) -c sisim.fpp
	$(FC) $(FFLAGS) $(OPENMP) -I. sisim.o -o sisim_OpenMP $(GSLIB)
	cp sisim_OpenMP ../bin 

openmp_mpi_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(MPIFC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -c sisim.fpp
	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -I. sisim.o -o sisim_OpenMP_MPI $(GSLIB)
	cp sisim_OpenMP_MPI ../bin 

#ins:
#	cd $(HOME)/gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -g -c extrae_module.f90
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -g -c sisim.fpp
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -g sisim.o -o sisim.exe $(GSLIB) -lgomp $(FLIBS)

clean:
	rm sisim_Seq sisim_OpenMP sisim_OpenMP_MPI ../bin/* *.o *.mod
