#include "itensor/all.h"
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "completeEntanglementTMUpdated.h"

using namespace itensor; 

// IN THIS CODE WE CALCULATE THE ENTANGLEMENT ENTROPY
// FIRST WE REGAUGE THE MATRICES AND BUT USE THE ORIGINAL PSI(0) VALUES
// FOR THE EIGENVALUES AS THESE DO CORRESPOND TO THE HMMMMM I NEED TO THINK 
// ABOUT THIS VERY CAREFULLY AS IT IS QUITE COMPLEX

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


    std::string savePath = "../../CHIRAL2/SavedFiles";
    std::string fullPath = savePath + "/" + saveFileName;

    auto NuC = inputParameters.getInt("Ny"); 
    auto T = inputParameters.getInt("T");
    auto maxBondDim = inputParameters.getInt("maxBondDim");

    // Reading the MPS using a 
    SpinHalf sites;
    readFromFile(fullPath+"/sitesFinalFile",sites);
    MPS psi(sites);
    readFromFile(fullPath+"/psiFinalFile",psi);
    ITensor R;
    readFromFile(fullPath+"/RTensorFile",R);
    
    std::vector<ITensor> Psi = orthogonalRight(psi); 

    for (int i = NuC+2; i <=2*NuC; i++){
        Psi.insert(Psi.begin(),psi(i));
    }
    // WE REVERSE THE ORDER THAT WE HAVE PUT THE VECTORS IN AHHHHH FFS 
    std::reverse(Psi.begin(),Psi.end());

    ITensor rhoTilde  = rhoTildeConstruct(Psi); 
    std::vector<ITensor> matrix = exactDiagonalisationH(rhoTilde);
    ITensor Q = matrix[0];
    ITensor D = matrix[1];

    momentumResolvedEigenSolver(R, Q, D, fullPath);

     
    return 0; 
}





