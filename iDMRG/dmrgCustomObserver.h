// IT IS THE MEASURE METHOD IN THE CUSTOMOBSERVER CLASS WHICH IS A DERVIED CLASS OF THE PUBLIC DMRGOBSERVER 
// METHODS AND ATTRIBUTES. IT IS THE DERIVED METHOD CUSTOMOBSERVER::MEASURE (WHICH BASICALLY MEANS ACCESS THE
// MEASURE METHOD IN CUSTOMOBSERVER WHICH WE HAVE CREATED) WHICH DETERMINES WHAT IS PRINTED OUT AT EACH STEP. 
// IDEA HERE WILL BE TO SAVE THE ENTANGLEMENT ENTROPY, ENERGY AND OVERLAPS TO A FILENAMEINTERMEDIATE VALUES FILE 
// THIS SHOULD GIVE AN IDEA OF CONVERGENCE AND CAN BE USED TO SHOW THAT WE HAVE ACHIEVED THE CORRECT CONVERGENCE
// WE CAN ALSO DEFINE CHECK WHICH RETURNS A BOOLEAN AND TELLS US WHETHER THE SOLUTION HAS PROPERLY CONVERGED

#ifndef __CUSTOM_DMRG_OBSERVER_H
#define __CUSTOM_DMRG_OBSERVER_H

#include <iomanip>
#include <fstream> 

using namespace itensor;

class CustomObserver : public DMRGObserver{
    public:

    CustomObserver(MPS const& psi,
                   Args const& args = Args::global())
        : DMRGObserver(psi)
        {}

    void virtual
    measure(Args const& args);

}; // class CustomObserver



void CustomObserver::
measure(Args const& args = Args::global()){

    auto& psi = DMRGObserver::psi();
    DMRGObserver::measure(args);

    auto N = length(psi);
    auto b = args.getInt("AtBond",1);
    auto sw = args.getInt("Sweep",0);
    auto nsweep = args.getInt("NSweep",0);
    auto ha = args.getInt("HalfSweep",0);
    auto energy = args.getReal("Energy",0);
    auto silent = args.getBool("Silent",false);
    auto iFileName = args.getString("IntermediateSaveFileName");
    auto path = args.getString("Path");

    std::ofstream intermediateMyFile; 
    std::string saveFileString = path + iFileName;
    intermediateMyFile.open(saveFileString, std::ios_base::app);


    if(!silent && printeigs)
        {
        if(b == N/2 && ha == 2)
            {
            println();
            auto center_eigs = last_spec_.eigsKept();
            // Normalize eigs
            Real norm_eigs = 0;
            for(auto& p : center_eigs)
            norm_eigs += p;
            center_eigs /= norm_eigs;
            // Calculate entropy
            //Real S = 0;
            S = 0.0;
            for(auto& p : center_eigs)
                {
                if(p > 1E-13) S += p*log(p);
                }
            S *= -1;
            intermediateMyFile << std::setprecision(12) <<  S << ";";
            auto ten = decltype(center_eigs.size())(10);
            for(auto j : range(std::min(center_eigs.size(),ten)))
                {
                auto eig = center_eigs(j);
                if(eig < 1E-3) break;
                intermediateMyFile << std::setprecision(6) << eig << ",";
                }
            intermediateMyFile << ";";
           
            }
        }

    max_eigs = std::max(max_eigs,last_spec_.numEigsKept());
    max_te = std::max(max_te,last_spec_.truncerr());
    if(!silent)
        {
        if(b == 1 && ha == 2) 
            {
            if(!printeigs) println();
            auto swstr = (nsweep>0) ? format("%d/%d",sw,nsweep) 
                                    : format("%d",sw);
            intermediateMyFile << max_te << ";"; 
            intermediateMyFile << max_eigs <<";";
            max_eigs = -1;
            max_te = -1;
            intermediateMyFile << std::setprecision(12) << energy << ";"; 
            intermediateMyFile << std::setprecision(12) << energy/N << "; \n"; 
            }
        }

    intermediateMyFile.close();

}


#endif