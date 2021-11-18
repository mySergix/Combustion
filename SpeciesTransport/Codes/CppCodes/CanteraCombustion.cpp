//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR READ INPUT DATA CLASS                                 //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/CanteraCombustion.h"

using namespace std;

//Constructor del lector de datos
CanteraCombustion::CanteraCombustion(Memory M1){

    ch4 = Cantera.Species('CH4', 'C:1', 'H:4');
    
}