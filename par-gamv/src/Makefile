

MPIFC=mpif90
FC=gfortran
FFLAGS= -cpp -O3 
OPENMP=-fopenmp
#FLIBS= -L$(EXTRAE_HOME)/lib -lomptrace 
GSLIB_HOME=../../gslib
GSLIB=$(GSLIB_HOME)/gslib.a


all: seq openmp_par openmp_mpi_par

seq: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd -
	$(FC) $(FFLAGS) -c gamv.fpp
	$(FC) $(FFLAGS) -I. gamv.o -o gamv_Seq $(GSLIB)
	cp gamv_Seq ../bin

openmp_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(FC) $(FFLAGS) $(OPENMP) -c gamv.fpp
	$(FC) $(FFLAGS) $(OPENMP) -I. gamv.o -o gamv_OpenMP $(GSLIB)
	cp gamv_OpenMP ../bin

openmp_mpi_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(MPIFC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(MPIFC) $(FFLAGS) -DUSE_MPI=1024 $(OPENMP) -c gamv.fpp
	$(MPIFC) $(FFLAGS) -DUSE_MPI=1024 $(OPENMP) -I. gamv.o -o gamv_OpenMP_MPI $(GSLIB)
	cp gamv_OpenMP_MPI ../bin

#ins:
#	cd $(HOME)/gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -c extrae_module.f90
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -c gamv.fpp
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -c extractStatisticsFortran.fpp
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE gamv.o extractStatisticsFortran.o -o gamv.exe $(GSLIB) $(FLIBS) -lgomp  

clean:
	rm gamv_Seq gamv_OpenMP gamv_OpenMP_MPI ../bin/* *.o *.mod
