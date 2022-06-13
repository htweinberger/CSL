#include "itensor/all.h"
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include "correlations.h"

// sys contains this header file and this is used to check file/directory
// Use stat() and pass address of a struct stat then check if its member 
// st_mode has S_IFDIR flag set
#include <sys/stat.h>


int main(int argc, char* argv[])
{

     if (argc != 2){
        // Need to pass a parameter and sweeps input file
        std::cerr << "Usage: " << argv[0] << " NAME" << std::endl; 
        return 1; 
    }

    auto inputParameters = InputGroup(argv[1],"input");
    auto saveFileName = inputParameters.getString("SaveFileName");
    auto J1 = inputParameters.getReal("J1"); 
    auto J2 = inputParameters.getReal("J2"); 
    auto JC = inputParameters.getReal("JC");


    std::string savePath = "/home/vol00/scarf1110/simulations/CHIRAL2/SavedFiles";
    std::string fullPath = savePath + "/" + saveFileName;

    auto NuC = inputParameters.getInt("Ny"); 
    auto T = inputParameters.getInt("T");
    auto maxBondDim = inputParameters.getInt("maxBondDim");

    // Reading the MPS using a 
    SpinHalf sites;
    readFromFile(fullPath+"/sitesFinalFile",sites);
    MPS psi(sites);
    readFromFile(fullPath+"/psiFinalFile",psi);


    int xrange = 600; 


    // CORRELATION FUNCTION WHICH MULTIPLIES IN THE EDGE STATE PSI(0) AND 
    // THEN CALCULATES THE CORRELATION BETWEEN SITES 1 AND J 

    std::vector<corInd> ijCorr = pairwiseCorrelations(psi,sites,xrange,fullPath); 

    return 0; 
}
    


