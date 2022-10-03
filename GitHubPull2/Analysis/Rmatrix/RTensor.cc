#include "itensor/all.h"
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "completeEntanglement.h"
#include <complex>
#include <filesystem>
#include <cstdlib> 



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
    std::filesystem::create_directory(fullPath + "/RTensor");

    auto NuC = inputParameters.getInt("Ny"); 
    auto T = inputParameters.getInt("T"); //Number of integer rotations {0,1,..,Ny-1}
    auto maxBondDim = inputParameters.getInt("maxBondDim");

    // Reading the MPS using a 
    SpinHalf sites;
    readFromFile(fullPath+"/sitesFinalFile",sites);
    MPS psi(sites);
    readFromFile(fullPath+"/psiFinalFile",psi);

    // WE START BY CHECKING WHETHER THE ORTHOGONALITY CONDITIONS ARE SATISFIED BY THE MPS
    //checkLeftOrthogonality(psi);
    //checkRightOrthogonality(psi);
    // WHAT WE SEE IS THAT EVERYTHING IS IN THE LEFT ORTHOGONALITY GAUGE EXCEPT THE STATE
    // WHICH IS CONTAINED IN PSI(NuC+1) WHICH IS BAFFLING 
    // WE THUS MUST NEVER MULTIPLY IN PSI(NuC+1) INTO ANY OF THE TRANSFER MATRICES 

    std::vector<ITensor> Psi = orthogonalRight(psi); 

    for (int i = NuC+2; i <=2*NuC; i++){
        Psi.insert(Psi.begin(),psi(i));
    }

    // for (auto& i : Psi){
    //     for (auto& j : inds(i)){
    //         println(j);
    //     }
    // }

    // if(dim(inds(psi(0))[0]) < 15){
    //     for (int ni=0; ni <= 2*NuC; ++ni){
    //         std::cout << "--------------- RIGHT ORTHOGONALITY OF: " << 2*NuC - ni << "--------------" << std::endl; 
    //         ITensor pa =Psi[ni];
    //         Index indc = commonIndex(pa,Psi[ni+1]);
    //         PrintData(pa*prime(dag(pa),indc));
    //     } 
    // }

    std::reverse(Psi.begin(),Psi.end()); 
    //PrintData(Psi[1]);

    // GET THE LIST OF ROTATED TENSORS 

    std::vector<std::array<Index,2>> list = translateIndices(psi,T);



    // WE NOW HAVE A NEW CONSTRUCTION OF THE MPS IN A PURELY RIGHT ORTHOGONAL GAUGE 

    Index indexLIn = inds(Psi[1])[2];
    Index indexLOut = inds(dag(prime(Psi[1],"Link")))[2]; 
    Index indexROut = uniqueIndex(Psi[2*NuC],Psi[2*NuC-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    std::vector<Index> edgeIndices;
    edgeIndices.push_back(indexLIn);
    edgeIndices.push_back(indexLOut);
    edgeIndices.push_back(indexROut);
    edgeIndices.push_back(indexRIn);

    auto [L,indexComL] = combiner(indexLOut,indexLIn);
    auto [R,indexComR] = combiner(indexROut, indexRIn);

    Index indexXY = sim(dag(indexComR));

    auto xR1 = randomITensorC(
                            QN({"Sz",0}),
                            sim(dag(indexROut)),sim(dag(indexRIn))
                            );
    
    auto xR2 = randomITensorC(
                            QN({"Sz",0}),
                            sim(dag(indexROut)),sim(dag(indexRIn))
                            );

    auto xR3 = randomITensorC(
                            QN({"Sz",0}),
                            sim(dag(indexROut)),sim(dag(indexRIn))
                            );

    // check that we obey the right canonical representation such that \sum_{s} A^{s}_{R} L A^{s}_{R} = L right canonical representations
    // Sadly we cannont explicitely construct the matrix 

    auto RM = rightMap(Psi,T,list);
    auto T_R = rightTransferMatrix(Psi);
    auto L_T = leftTransferMatrix(Psi);

    auto lambdaR = arnoldi(RM,xR1,{"ErrGoal=",1E-20,"MaxIter=",40,"MaxRestart=",10});
    auto lambdaTR = arnoldi(T_R,xR2,{"ErrGoal=",1E-20,"MaxIter=",40,"MaxRestart=",10});
    auto lambdaLT = arnoldi(L_T,xR3,{"ErrGoal=",1E-20,"MaxIter=",40,"MaxRestart=",10});

    PrintData(lambdaR);
    PrintData(lambdaTR);
    PrintData(lambdaLT);
    

    std::ofstream myfile; 
    myfile.open(fullPath + "/RTensor/entanglementEntropy.txt"); 
    myfile << "J1:" << J1 << std::endl; 
    myfile << "J2:" << J2 << std::endl;
    myfile << "JC:" << JC << std::endl;
    myfile << "Mixed Transfer Matrix Eigenvalue:" << lambdaR << std::endl; 
    myfile << "Right Eigenvalue:" << lambdaTR << std::endl; 
    myfile << "Left Eigenvalue:" << lambdaLT << std::endl;
    myfile.close();

    ///////////////////////////////////////////////////////////////////////////////////////////
    //////////////////// NOW WE WANT DENSE REPRESENATIONS OF THE QNMPS TO /////////////////
    //////////////////// EXTRACT AS MATRICES AND ANALYSE USING PYTHON     /////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    
    Index n = inds(xR1)[0];
    Index m = inds(xR1)[1];

    // SAVE THE ITENSOR TO FILE 
    writeToFile(fullPath + "/RTensorFile",xR1);

    
    std::ofstream RMatrix; 
    RMatrix.open(fullPath + "/RTensor/RMatrix.txt");
    
    for(int i=1; i <= dim(n); i++){
    	for(int j=1; j <= dim(m); j++){
		
		IndexVal iIndex = IndexVal(n,i);
                IndexVal jIndex = IndexVal(m,j);
		QN fluxi = qn(iIndex); 
		QN fluxj = qn(jIndex);
		
		if(j < dim(n)){
			RMatrix << "[" ;
			RMatrix << eltC(xR1,iIndex,jIndex) ;
			RMatrix << "," ; 
			RMatrix << fluxi.val("Sz");
			RMatrix <<"," ;
			RMatrix << fluxj.val("Sz"); 
			RMatrix << "]";
			RMatrix << ";";
		}

		else{
			RMatrix << "[" ;
			RMatrix << eltC(xR1,iIndex,jIndex) ;
			RMatrix << "," ; 
			RMatrix << fluxi.mod("Sz");
			RMatrix <<"," ;
			RMatrix << fluxj.mod("Sz"); 
			RMatrix << "]";
			RMatrix << std::endl; 
		}
 
	}
    }

    RMatrix.close();

    //////////////////////////////////////////////////////////////////////////////
    ///////////// SEE IF THIS ROTATION MATRIX ACTUALLY COMPUTES A ROTATION ///////
    //////////////////////////////////////////////////////////////////////////////

    std::ofstream checkRotationFile;
    checkRotationFile.open(fullPath+"/RTensor/rotationQN.txt");

    checkRotationFile << "---------- Right orthogonal gauge psi ----------" << std::endl;
    for (int i = 1; i < NuC ; i++){
	  checkRotationFile << overlapR(psi,translateIndices(psi,i)) << std::endl;
	
    }
    checkRotationFile << "---------- Right orthogonal after regauging ----------" << std::endl;
    for (int i = 1; i < NuC ; i++){
	  checkRotationFile << overlapR(Psi,translateIndices(psi,i)) << std::endl;
	
    }

    checkRotationFile << normalisation(psi) << std::endl; 
    checkRotationFile << normalisation(Psi) << std::endl;


    /////////////////////////////////////////////////////////////////////////////
    ///////////// CHECK R REPRESENTATION AS A ROTATION OF THE ///////////////////
    //////////// ITENSORS IN REAL SPACE                       ///////////////////
    /////////////////////////////////////////////////////////////////////////////

    checkRotationFile << "-------------- BOND REPRESENTATION OF THE MATRICES ------------" << std::endl; 
    checkRotationFile << checkRRep(psi,xR1);

    checkRotationFile.close();
    return 0;

}

