#include "itensor/all.h"
#include "idmrg.h"
#include "MPOconstruction.h"
#include "dmrgCustomObserver.h"
#include "randomGen.h"
#include <array>
#include <string>
#include <fstream>

using namespace itensor;

const bool initialSetup = true;
const bool infinite = true;


int main(int argc, char* argv[])
    {
    
    if (argc != 2)
    {
        // Need to pass a parameter and sweeps input file
        std::cerr << "Usage: " << argv[0] << " NAME" << std::endl; 
        return 1; 
    }

    // Parameter declaration
    auto inputParameters = InputGroup(argv[1],"input"); 
    auto NuC = inputParameters.getInt("Ny"); 
    auto multiplicity = inputParameters.getInt("Multiplicity"); 
    auto J1 = inputParameters.getReal("J1"); 
    auto J2 = inputParameters.getReal("J2"); 
    auto JC = inputParameters.getReal("JC");
    auto saveFileName = inputParameters.getString("SaveFileName");
    auto intermediateSaveFileName = inputParameters.getString("IntermediateSaveFileName");


    // MPO Construction     
    int N = multiplicity*NuC;
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    MPO H = Heisenberg(sites,{"Infinite=",infinite,"J2=",J2,"J1=",J1,"Multiplicity=",multiplicity,"JC=",JC});
    auto state = InitState(sites);

    // Sweep Object Construction 
    auto nSweeps = inputParameters.getInt("NSweeps");
    auto sw_table = InputGroup(inputParameters,"table_name");
    auto sweeps = Sweeps(nSweeps,sw_table);

    // InitialState setup 
    if(initialSetup){
        const std::array<std::string,2> statesSz = {"Up","Dn"}; 
        auto initialStates = randomStateSz(statesSz,N); 

        for(int i = 1; i <= N; ++i){
        state.set(i,initialStates[i-1]);
        }
    }
    else{

        for(int i = 1; i <= N; ++i)
        {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        }

    }
    
    auto psi = MPS(state);

    //Save the initial state we need to have these folders already made
    std::string path = "SavedFiles/" + saveFileName + "/";
    writeToFile(path + "sitesInitialFile",sites);
    writeToFile(path + "psiInitialFile",psi);

    //Create dmrgObserver object 
    CustomObserver observer(psi); 

    // idmrg returns a struct holding various useful
    // things such as the energy and the "edge tensors"
    // representing the Hamiltonian projected into the infinite MPS
    auto res = idmrg(psi,H,sweeps,observer,{"OutputLevel",1,"IntermediateSaveFileName", intermediateSaveFileName,"Path",path});

    printfln("\nGround state energy / site = %.20f",res.energy/N);

    //Interesting to compare ground state energy to White, Huse, PRB 48, 3844 (1993).

    //
    //Measure correlation function by repeating infinite MPS unit cell 
    //

    //Multiply in psi(0) which holds singular values
    auto wf1 = psi(0)*psi(1); 
    //oi is the outer IQIndex "sticking out" of the left edge of psi(0)
    auto oi = uniqueIndex(psi(0),psi(1),"Link");
    //lcorr is the left side of the correlation function tensor
    //which grows site by site below
    auto lcorr = prime(wf1,oi)*op(sites,"Sz",1)*dag(prime(wf1));

    std::ofstream myfile2; 
    myfile2.open(path + "/correlationsSz.txt"); 

    println("\nj <psi|Sz_1 Sz_j|psi> = ");
    //xrange is how far to go in measuring <Sz_1 Sz_j>, 
    //ok to adjust xrange to any size >= 2
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

    
    //auto f = h5_open(path+"hdf5MPS.h5",'w'); //open HDF5 file in write 'w' mode
    //h5_write(f,"MPS",psi); 
    //close(f);

    writeToFile(path + "sitesFinalFile",sites);
    writeToFile(path + "psiFinalFile", psi);



    return 0;
}