// HEADER FILE TO CALCULATE I,J CORRELATIONS FOR ALL 
// LATTICE POINTS IN THE UNIT CELL 

#ifndef __CORRELATIONS
#define __CORRELATIONS

#include <vector> 
#include <array>
#include <iostream> 
#include <string> 
#include <complex>


inline std::vector<std::vector<std::string>> pauliOp = {{"Sz","Sz"}}; //{{"Sz","Sz"},{"Sp","Sm"},{"Sm","Sp"}};

class corInd{

    public: 

        int latticeSite_i; 
        int latticeSite_j; 

        // WE WILL ADD THE VALUES OF THE CORRELATIONS BETWEEN I AND J AS 
        // {SZ,SZ}, {SP,SM}, {SM,SP} WITH NO FACTORS OF HALF
        std::vector<std::complex<double>> correlators;   

        // INITIALISE KNOWN VALUES USING BASE CONSTRUCTOR
        corInd() : latticeSite_i(0), 
                latticeSite_j(0){}; 
        corInd(int,int); 

    
        void addCorrelation(std::complex <double> correlator){
            correlators.push_back(correlator);
        }

};

inline corInd::
corInd( int x, int y) : 
    latticeSite_i(x),      
    latticeSite_j(y){};   

using namespace itensor; 

// RETURN VECTOR CONTAINING CORRELATION INDEX OBJECTS
std::vector<corInd> pairwiseCorrelations(MPS psi, SiteSet const& sites, int xrange,std::string path){

    int N = length(psi); 
    // USE A VECTOR TO STORE THE CORIND OBJECTS BETWEEN SITES I AND J  
    std::vector<corInd> pairCorrelations;

    //INITIALISE THE ARRAY AND CORRELATION INDICES WITH THEIR VALUES
    for (int i = 1; i <= 1; i++){
        for(int j = 1; j <= xrange; j++){

            if(j>i){
                pairCorrelations.push_back(corInd(i,j)); 
            }
        
        }
    }

    std::ofstream outfile;
    outfile.open(path + "/correlations.txt", std::ios_base::app);//std::ios_base::app
    outfile << "\n ----------------------- NEW RUN --------------------------- \n ";

    // CLAIM: WE HAVE CYLINDRICAL SYMMETRY AND SO WE DO NOT 
    // WE ALSO HAVE ROTATIONAL SYMMETRY AND SO WE SHOULD ONLY NEED 
    // TO CONSIDER 1-N AND SO WE DO NOT NEED TO WORRY ABOUT 
    

    for (std::vector<std::string> pauli : pauliOp)
    {
        int counter = 0;
        int position = 0;
        auto wf1 = psi(0)*psi(1); 
        auto oi = uniqueIndex(psi(0),psi(1),"Link");
        auto SzOp = op(sites,pauli[0],1);

        auto lcorr = prime(wf1,oi)*SzOp*dag(prime(wf1));

        for (auto& i : inds(lcorr)){
            println(i);
        }

        
        for(int j = 2; j <= xrange; ++j){
            int n = (j-1)%N+1; 
            auto ui = uniqueIndex(psi(n),lcorr,"Link");

            auto SzOp2 = op(sites,pauli[1],n);

            pairCorrelations[position].addCorrelation(eltC(dag(prime(psi(n)))*lcorr*prime(psi(n),ui)*SzOp2));   
            lcorr *= psi(n);
            lcorr *= dag(prime(psi(n),"Link"));
            std::cout << "---" << std::endl;
            outfile << pairCorrelations[position].correlators[0]<< std::endl;
            position+=1;
        }
        counter += 1;
        
    }

    outfile.close();

    return pairCorrelations;   
}
    
    
#endif