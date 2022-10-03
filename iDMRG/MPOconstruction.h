//Next-Nearest-Neighbour MPO generator

#ifndef __MPO_CONSTRUCTION_H
#define __MPO_CONSTRUCTION_H

#include "itensor/mps/mpo.h"
#include "dressedNNObject.h"
#include "dressedChiralFinder.h"
#include <iostream>
#include <fstream>

namespace itensor {

class Heisenberg
    {
    public:

    Heisenberg(SiteSet const& sites, 
               Args const& args = Args::global());

    operator MPO() { init_(); return H; }

    private:

    //////////////////
    //
    // Data Members

    SiteSet const& sites_;
    int N_;
    int multiplicity_;
    Real J1_, 
         J2_,
         JC_;
    bool initted_,
         infinite_;
    MPO H;

    //
    //////////////////

    void 
    init_();

    }; //class Heisenberg

inline Heisenberg::
Heisenberg(SiteSet const& sites, 
           Args const& args)
  : sites_(sites), // class inheritance
    initted_(false)
    { 
    N_ = sites_.length();
    J1_ = args.getReal("J1",1.);
    J2_ = args.getReal("J2",0.);
    JC_ = args.getReal("JC",0.);
    infinite_ = args.getBool("Infinite",false);
    multiplicity_ = args.getInt("Multiplicity",2);
    }

void inline Heisenberg::
init_()
    {
    if(initted_) return;

    int Ny_ = N_/multiplicity_;

    H = MPO(sites_);

    std::vector<Index> links(N_+1);

    for(int l = 0; l <= N_; ++l) 
        {
        auto ts = format("Link,l=%d",l);
        links.at(l) = Index(QN({"Sz", 0}),3*Ny_-1+2,
                            QN({"Sz",+2}),3*Ny_-1,
                            QN({"Sz",-2}),3*Ny_-1,
                            Out,
                            ts);
        }
    
    std::ofstream MyFile("mpoTerms.txt");

    Index const& last = (infinite_ ? links.at(0) : links.at(N_));

    for(int n = 1; n <= N_; ++n)
    {
        
        #if 0
        std::cout << n << std::endl; 
        #endif 

        auto& W = H.ref(n);
        auto row = dag(links.at(n-1));
        auto col = (n==N_ ? last : links.at(n));

        W = ITensor(dag(sites_(n)),prime(sites_(n)),row,col);

        

        MyFile << "---------MPO TERM: " << std::to_string(n) << " -----------" << std::endl;

        auto chiralTerms = chiralHelperFunction(n,JC_,Ny_);

        MyFile << "--------ALL CHIRAL TERMS --------"<<std::endl;

        for(chiralCont wIndex : chiralTerms){
            W += sites_.op(wIndex.opType,n) * setElt(row(wIndex.row)) * setElt(col(wIndex.col)) * wIndex.cConst;
            MyFile << wIndex.opType << "," << wIndex.row << "," << wIndex.col << "," << wIndex.cConst << std::endl; 
        }

    

        MyFile << "--------ALL NN TERMS --------"<<std::endl;

        auto nnTerms = nnHelperFunction(n,J1_,Ny_);
        for(nnCont wIndex : nnTerms){
            W += sites_.op(wIndex.opType,n) * setElt(row(wIndex.row)) * setElt(col(wIndex.col)) * wIndex.cConst;
            MyFile << wIndex.opType << "," << wIndex.row << "," << wIndex.col << "," << wIndex.cConst << std::endl; 
        }



        MyFile << "--------ALL NNN TERMS --------"<<std::endl;

        auto nnnTerms = nnnHelperFunction(n,J2_,Ny_);
        for(nnCont wIndex : nnnTerms){
            W += sites_.op(wIndex.opType,n) * setElt(row(wIndex.row)) * setElt(col(wIndex.col)) * wIndex.cConst;
            MyFile << wIndex.opType << "," << wIndex.row << "," << wIndex.col << "," << wIndex.cConst << std::endl; 
        }

        MyFile << "--------ALL Identity TERMS --------"<<std::endl;

        auto additionalTerms = idHelperFunction(Ny_);
        for(nnCont wIndex : additionalTerms){

            W += sites_.op(wIndex.opType,n) * setElt(row(wIndex.row)) * setElt(col(wIndex.col)) * wIndex.cConst;
            MyFile << wIndex.opType << "," << wIndex.row << "," << wIndex.col << "," << wIndex.cConst << std::endl; 
        }

    }

    MyFile.close();

    auto LH = setElt(links.at(0)(2));
    auto RH = setElt(dag(last)(1));

    if(not infinite_)
        {
        //Multiply first and last
        //MPO tensor by edge vectors
        H.ref(1) *= LH;
        H.ref(N_) *= RH;
        }
    else
        {
        //Store edge vectors just before
        //and after first and last sites
        H.ref(0) = LH;
        H.ref(N_+1) = RH;
        }

    initted_ = true;
    }

} //namespace itensor

#endif