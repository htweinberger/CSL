#include "itensor/all.h"
#include "idmrg.h"
#include "MPOconstruction.h"
#include "dmrgCustomObserver.h"
#include "randomGen.h"
#include <array>
#include <string>
#include <fstream>
#include <cstdlib>
#include <filesystem>

using namespace itensor;

int main(int argc, char* argv[])
    {
    
    //Needs a parameter file as an argv argument pointer
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " NAME" << std::endl; 
        return 1; 
    }

    //Parameter declaration and using input.h file 
    //Each of these must be declared in the parameter file which is passed as an argument pointer
    InputGroup inputParameters = InputGroup(argv[1],"input"); 
    int NuC = inputParameters.getInt("Ny"); //Number of sites in one ring of the cylinder
    int multiplicity = inputParameters.getInt("Multiplicity"); //Number of rings in one unit cell 
    Real J1 = inputParameters.getReal("J1"); //Real are defined in iTensor namespace
    Real J2 = inputParameters.getReal("J2"); 
    Real JC = inputParameters.getReal("JC");
    std::string saveFileName = inputParameters.getString("SaveFileName"); 
    std::string intermediateSaveFileName = inputParameters.getString("IntermediateSaveFileName");
    bool initialSetup = inputParameters.getYesNo("RandomSetup"); //understands yes/no and returns a boolean
    bool infinite = inputParameters.getYesNo("Infinite");

    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////// Customisation by looking at file input.h /////////////////////////
    ///////////////////// as can read into a specific variable     /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////


    //Construction of MPO    
    int N = multiplicity*NuC; //Number of sites in the unit cell for iDMRG 
    SiteSet sites = SpinHalf(N,{"ConserveQNs=",true});
    MPO H = Heisenberg(sites,{"Infinite=",infinite,"J2=",J2,"J1=",J1,"Multiplicity=",multiplicity,"JC=",JC});
    auto state = InitState(sites);

    //Construction of sweep objects
    int nSweeps = inputParameters.getInt("NSweeps");
    auto sw_table = InputGroup(inputParameters,"table_name");
    auto sweeps = Sweeps(nSweeps,sw_table);

    //Save the initial state we need to have these folders already made
    std::string path = "../SavedFiles/" + saveFileName + "/";
    std::ofstream initialStateFile;
    std::filesystem::create_directory(path + "iDMRG");
    initialStateFile.open (path + "iDMRG/initialStateFile.txt");
    

    // InitialState setup using an overengineered Fisher-Yates algorithm 
    // Use cases, we need inputs to be std::strings however the randomGen.h 
    // file can work with any data format provided that N%len == 0 where N 
    // in this case is the number of lattice points in the unit cell.  

    if(initialSetup){

        std::string statesSz[] = {"Up","Dn"}; 
        constexpr int len = sizeof(statesSz)/sizeof(statesSz[0]);
        std::vector<std::string> initialStates = randomStatesSz<std::string,len>(statesSz,N).getRandomNumbers();

        for(int i = 1; i <= N; ++i){
            state.set(i,initialStates[i-1]);
            initialStateFile << "Site ";
            initialStateFile << i;
            initialStateFile << ": ";
            initialStateFile << initialStates[i-1] << std::endl; 
        }
    }
    else{

        for(int i = 1; i <= N; ++i)
        {
            if((i%2) == 1){
                state.set(i,"Up");
                initialStateFile << "Site ";
                initialStateFile << i;
                initialStateFile << ": ";
                initialStateFile << "Up" << std::endl; 
            }
            else{
                state.set(i,"Dn");
                initialStateFile << "Site ";
                initialStateFile << i;
                initialStateFile << ": ";
                initialStateFile << "Dn" << std::endl; 
            }

        }
    }

    initialStateFile.close();
    
    MPS psi = MPS(state);
    writeToFile(path + "sitesInitialFile",sites);
    writeToFile(path + "psiInitialFile",psi);

    // dmrgObserver object allows you to save intermediate data more easily 
    CustomObserver observer(psi); 

    // idmrg returns a struct holding various useful
    // things such as the energy and the "edge tensors"
    // representing the Hamiltonian projected into the infinite MPS
    // also takes in the CustomObserver object which has the intermediate
    // savefilename

    auto res = idmrg(psi,H,sweeps,observer,{"OutputLevel",
                                            1,
                                            "IntermediateSaveFileName",
                                            intermediateSaveFileName,
                                            "Path",
                                            path+"iDMRG/"});

    
    // average energy per lattice site 
    printfln("\nGround state energy / site = %.20f",res.energy/N);

    //Measure correlation function by repeating infinite MPS unit cell 
    //

    //Multiply in psi(0) which holds singular values betweent he bipartition of the infinite cylinder
    //Not the singular values of the bipartition are topologically protected since they also correspond 
    //To the left eigenvector of unity of the right transfer matrix and therefore do not decay --> global 
    //Property

    auto wf1 = psi(0)*psi(1); 
    //oi is the outer IQIndex "sticking out" of the left edge of psi(0)
    auto oi = uniqueIndex(psi(0),psi(1),"Link");
    //lcorr is the left side of the correlation function tensor
    //which grows site by site below
    auto lcorr = prime(wf1,oi)*op(sites,"Sz",1)*dag(prime(wf1));

    std::ofstream myfile2; 
    myfile2.open(path + "/iDMRG/correlationsSz.txt"); 

    println("\nj <psi|Sz_1 Sz_j|psi> = ");
    //xrange can be anything >= 2 and represents the connected correlation function
    // <Sz_i Sz_j> where i=1 and j=600 in the 1D mapped system. 
    int xrange = 600; 
    for(int j = 2; j <= xrange; ++j)
        {
        int n = (j-1)%N+1; //translate from j to unit cell site number
        //ui is the IQIndex "sticking out" of the right edge of psi(n)
        auto ui = uniqueIndex(psi(n),lcorr,"Link");
        //prime ui so it contracts with the "bra" tensor on top = dag(prime(psi(n)))
        auto val = eltC(dag(prime(psi(n)))*lcorr*prime(psi(n),ui)*op(sites,"Sz",n));
        
	myfile2 << val << std::endl;
        printfln("%d %.20f",j,val);
        lcorr *= psi(n);
        lcorr *= dag(prime(psi(n),"Link"));
        }

    myfile2.close();

    //Currently not using HDF5 however this can be implemented
    //Save files locally as file type
    writeToFile(path + "sitesFinalFile",sites);
    writeToFile(path + "psiFinalFile", psi);

    return 0;
}