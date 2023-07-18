FC=gfortran
FFLAGS=-g -O0 -O3 -Wall -Wextra -std=f2008
SRC=KerrEquations.f90 RK_process.f90 ray.f90 main.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o main