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
	$(FC) $(FFLAGS) -c sgsim.fpp
	$(FC) $(FFLAGS) -I. sgsim.o -o sgsim_Seq $(GSLIB)
	cp sgsim_Seq ../bin 

openmp_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(FC) $(FFLAGS) $(OPENMP) -c sgsim.fpp
	$(FC) $(FFLAGS) $(OPENMP) -I. sgsim.o -o sgsim_OpenMP $(GSLIB)
	cp sgsim_OpenMP ../bin 

openmp_mpi_par: 
	cd $(GSLIB_HOME); make clean; make gslib.a COMPILER="$(MPIFC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -c sgsim.fpp
	$(MPIFC) $(FFLAGS) $(OPENMP) -DUSE_MPI -I. sgsim.o -o sgsim_OpenMP_MPI $(GSLIB)
	cp sgsim_OpenMP_MPI ../bin 

#ins:
#	cd $(HOME)/gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd -
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -c extrae_module.f90
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE -c sgsim.fpp
#	$(FC) $(FFLAGS) $(OPENMP) -I. -I$(EXTRAE_HOME)/include -DTRACE sgsim.o -o sgsim.exe $(GSLIB) $(FLIBS) -lgomp  

clean:
	rm sgsim_Seq sgsim_OpenMP sgsim_OpenMP_MPI ../bin/* *.o *.mod
