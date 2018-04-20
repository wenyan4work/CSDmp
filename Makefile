# intel compiler with MKL
# lp64 interface is enough. it supports matrix of approximately 45000*45000
# lapack and bals 95 are automatically handled by intel compiler
FC = ifort
FLAGS = -O3 -xHost -parallel -std95 -ipo -qopenmp -I$(MKLROOT)/include/intel64/lp64
FLINK =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
       -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# gcc compiler with BLAS & LAPACK
#BLAS=/usr/lib
#LAPACK=/usr/lib
#FC = gfortran
#FLAGS = -O3 -march=native -fopenmp -std=f95 -I$()
#FLINK = $(FLAGS)

SD.x: main.o lub.o mob.o prutil.o fts.o in_out.o
	$(FC) $(FLAGS) main.o lub.o mob.o prutil.o fts.o in_out.o -o SD.x $(FLINK)
prutil.mod: prutil.o prutil.f90
	$(FC) $(FLAGS) prutil.f90 -c
prutil.o: prutil.f90
	$(FC) $(FLAGS) prutil.f90 -c
in_out.mod: in_out.o in_out.f90
	$(FC) $(FLAGS) in_out.f90 -c
in_out.o: in_out.f90
	$(FC) $(FLAGS) in_out.f90 -c
mob.mod: mob.o mob.f90
	$(FC) $(FLAGS) mob.f90 -c
mob.o: prutil.mod mob.f90
	$(FC) $(FLAGS) mob.f90 -c
lub.mod: lub.o lub.f90
	$(FC) $(FLAGS) lub.f90 -c
lub.o: prutil.mod in_out.mod lub.f90
	$(FC) $(FLAGS) lub.f90 -c
fts.mod: fts.o fts.f90
	$(FC) $(FLAGS) fts.f90 -c
fts.o: prutil.mod fts.f90
	$(FC) $(FLAGS) fts.f90 -c
main.o: mob.mod lub.mod prutil.mod fts.mod in_out.mod main.f90 
	$(FC) $(FLAGS) main.f90 -c
clean:  
	rm ./*.mod
	rm ./*.o
