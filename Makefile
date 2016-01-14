### Compilers & flags
F90=gfortran

FFTWLIBS=
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = $(LAPACKLIB)


EXE = exec
F90SRC = main.f90 random.f90 declaration.f90 init.f90 modPlasma.f90 timeStep.f90 assignFunctions.f90 convergence.f90 MatrixVector.f90
F90OBJ = main.o random.o declaration.o init.o modPlasma.o timeStep.o assignFunctions.o convergence.o MatrixVector.o

### Targets
all: $(EXE)
run: $(EXE) 
	./$(EXE)

# Link object files to executables
$(EXE): $(F90OBJ)
	$(F90) -o $(EXE) $(F90OBJ) $(LIBS)

# All .o files depend on the corresponding .f90 file
%.o: %.f90
	$(F90) -c $<

# Dependencies
main.o: convergence.o
declaration.o : modPlasma.o
init.o: declaration.o random.o
timeStep.o : declaration.o assignFunctions.o MatrixVector.o
assignFunctions.o : declaration.o
convergence.o : declaration.o init.o modPlasma.o timeStep.o

clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


