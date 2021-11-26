#------------------------------------------------------------------------------------------------#
#                             MAKEFILE FOR RUNNING SPECIES TRANSPORT CASE                        #
#------------------------------------------------------------------------------------------------#

# Codes Names
MainCode = Main_SpeciesTransport
ExecutableName = Run_SpeciesTransport

Class_1 = Memory
Class_2 = ReadData
Class_3 = Parallel
Class_4 = Mesher
Class_5 = PostProcess
Class_6 = CFD_Solver
Class_7 = Species_Solver

Directorio_Clases = CppCodes/

default: $(MainCode)

Compilador = mpic++
Optimisation = -O

# Code Execution
$(MainCode): $(MainCode).cpp $(INC)
	$(Compilador) -o $(ExecutableName) $(MainCode).cpp $(Directorio_Clases)$(Class_1).cpp $(Directorio_Clases)$(Class_2).cpp $(Directorio_Clases)$(Class_3).cpp $(Directorio_Clases)$(Class_4).cpp $(Directorio_Clases)$(Class_5).cpp $(Directorio_Clases)$(Class_6).cpp $(Directorio_Clases)$(Class_7).cpp $(Optimisation)
	
#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN) $(ExecutableName)

run:
	
	mpirun -np 4 ./$(ExecutableName)