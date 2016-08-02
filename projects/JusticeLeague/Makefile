COMPILER = gfortran
FLAGS = 
LINEAR_ALGEBRA = -llapack

SRC = types.f90  HF_variables.f90 HF_HartreeFock.f90 HF_HOQuadrature.f90 HF_Minnesota.f90 HF_LDA.f90 HFSolver_main.f90
OBJ = $(SRC:.f90=.o)

all: HFSolver

HFSolver: $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LINEAR_ALGEBRA)

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

HF_variables.o: HF_variables.f90 types.mod
	$(COMPILER) $(FLAGS) -c $<

HF_HartreeFock.o: HF_HartreeFock.f90 types.mod variables.mod 
	$(COMPILER) $(FLAGS) -c $<

HF_HOQuadrature.o: HF_HOQuadrature.f90 types.mod
	$(COMPILER) $(FLAGS) -c $<	

HF_Minnesota.o: HF_Minnesota.f90 types.mod variables.mod ho_quadrature.mod
	$(COMPILER) $(FLAGS) -c $<

HF_LDA.o: HF_LDA.f90 types.mod variables.mod ho_quadrature.mod minnesota.mod
	$(COMPILER) $(FLAGS) -c $<

HFSolver_main.o: HFSolver_main.f90 types.mod variables.mod hartreefock.mod minnesota.mod lda.mod
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf HFSolver
