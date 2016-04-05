### Compilers & flags
F90=gfortran

FFTWLIBS=~/bin/FFTW/lib/libfftw3.a
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = $(LAPACKLIB) $(FFTWLIBS)


EXE = exec
F90SRC = main.f90 constants.f90 MatrixVector.f90 modPlasma.f90 modMesh.f90 modAssign.f90 modRecord.f90 modPM1D.f90 random.f90 init.f90 modSource.f90 timeStep.f90
F90OBJ = main.o constants.o MatrixVector.o modPlasma.o modMesh.o modAssign.o modRecord.o modPM1D.o random.o init.o modSource.o timeStep.o

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
MatrixVector.o : constants.o
modPlasma.o : constants.o
modMesh.o : MatrixVector.o
modAssign.o : modPlasma.o modMesh.o
modRecord.o : modPlasma.o modMesh.o
modPM1D.o : modPlasma.o modMesh.o modAssign.o modRecord.o
init.o : modPM1D.o random.o
modSource.o : modPM1D.o
timeStep.o : modSource.o
main.o : init.o timeStep.o

clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


