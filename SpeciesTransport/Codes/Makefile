#------------------------------------------------------------------------------------------------#
#                             MAKEFILE FOR RUNNING SPECIES TRANSPORT CASE                        #
#------------------------------------------------------------------------------------------------#

# Codes Names
MainCode = Main_SpeciesTransport
ExecutableName = Run_SpeciesTransport

Class_1 = Memory
Class_2 = ReadData

Directorio_Clases = CppCodes/

default: $(MainCode)

Compilador = mpic++
Optimisation = -O

# Code Execution
$(MainCode): $(MainCode).cpp $(INC)
	$(Compilador) -o $(ExecutableName) $(MainCode).cpp $(Directorio_Clases)$(Class_1).cpp $(Directorio_Clases)$(Class_2).cpp $(Optimisation)
	
#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) $(CLEAN) $(ExecutableName)

run:
	mpirun -np 4 ./$(ExecutableName)