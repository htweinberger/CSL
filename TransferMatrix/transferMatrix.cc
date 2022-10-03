#include "itensor/all.h"
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem> 
#include "arnoldiHelpers.h"

using namespace itensor; 
namespace fs = std::filesystem; 

// IN THIS CODE WE CALCULATE THE TRANSFER MATRIX BY THE ARNOLDI ITERATIVE SOLVER WHICH 
// FINDS DOMINANT AND SUB-DOMINANT EIGENVALUES. WE START BY CONTRACTING TRANSFER MATRIX 
// ONTO SOME RANDOM ITENSOR X WHICH CONSERVES SZ. 

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

    std::string fullPath = "../../SavedFiles/" + saveFileName;

    auto NuC = inputParameters.getInt("Ny"); 

    // Reading the MPS using a 
    SpinHalf sites;
    readFromFile(fullPath+"/sitesFinalFile",sites);
    MPS psi(sites);
    readFromFile(fullPath+"/psiFinalFile",psi);

    /// NOTE THAT I THINK WE ARE ALMOST CALCULATING THE LEFT EIGENVALUE AND THEN 
    // WE CAN EQUALLY CALCULATE THE OTHER EIGENVALUE

    PrintData(length(psi));

    auto s = psi(1);
    auto sdag = dag(prime(psi(1))); 

    auto index1 = inds(s)[0];
    auto index2 = inds(s)[1];
    auto index3 = inds(s)[2];

    auto index1dag = inds(sdag)[0];
    auto index2dag = inds(sdag)[1];
    auto index3dag = inds(sdag)[2];

    auto rIndex = rightLinkIndex(psi,NuC);
    auto rIndexdag = rightLinkIndex(dag(prime(psi,"Link")),NuC);


    auto [E,indexCom1] = combiner(rIndex, rIndexdag);
    auto [C,indexCom] = combiner(index3,index3dag);


    auto LM = leftMap(psi);
    auto RM = rightMap(psi);
    
    Index indexXY = sim(dag(indexCom));

    auto xL = randomITensor(QN({"Sz",0}),indexXY);
    auto yL = randomITensor(QN({"Sz",0}),indexXY);

    auto xR = randomITensor(QN({"Sz",0}),dag(indexXY));
    auto yR = randomITensor(QN({"Sz",0}),dag(indexXY));


    std::vector<ITensor> vecL = {xL,yL};
    std::vector<ITensor> vecR = {xR,yR};

    ITensor b;


    auto lambdaL = arnoldi(LM,vecL,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});
    auto lambdaR = arnoldi(RM,vecR,{"ErrGoal=",1E-14,"MaxIter=",20,"MaxRestart=",5});
    PrintData(lambdaL[0]);
    PrintData(lambdaL[1]);
    PrintData(lambdaR[0]);
    PrintData(lambdaR[1]);

    PrintData(norm(eigenCheck(vecL[0],psi) - lambdaL[0]*vecL[0]));
    PrintData(norm(eigenCheck(vecL[1],psi) - lambdaL[1]*vecL[1]));
	
    //Using std::filesystem in c++17 to create a new directory 
    fs::create_directory(fullPath + "/TransferMatrix");

    std::ofstream myfile; 
    myfile.open(fullPath + "TransferMatrix/transferMatrixEigenvalues.txt"); 
    myfile << "J1:" << J1 << std::endl; 
    myfile << "J2:" << J2 << std::endl;
    myfile << "JC:" << JC << std::endl;
    myfile << "Left Eigenvalue 1:" << lambdaL[0] << std::endl; 
    myfile << "Left Eigenvalue 2:" << lambdaL[1] << std::endl; 
    myfile << "Right Eigenvalue 1:" << lambdaR[0] << std::endl; 
    myfile << "Right Eigenvalue 2:" << lambdaR[1] << std::endl;
    myfile << "Error 0:" << norm(eigenCheck(vecL[0],psi) - lambdaL[0]*vecL[0]) << std::endl;
    myfile << "Error 1:" << norm(eigenCheck(vecL[1],psi) - lambdaL[1]*vecL[1]) << std::endl;
    myfile.close();


    return 0; 
}





