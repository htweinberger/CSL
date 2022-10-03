#ifndef __NNN_NEIGHBOURS_H
#define __NNN_NEIGHBOURS_H 

// AFTER GOING THROUGH THE OTHER NNN CODE I NOTICED THAT I HAD A COUPLE BUGS IN TERMS OF THE PLACEMENT OF THE OPERATORS IN THE 
// SCHUR CONSTRUCTION OF THE MPO WHICH WAS LEADING TO UNEXPECTED RESULTS. WITH THIS IN MIND I HAVE DECIDED TO TAKE IT UPON MYSELF 
// AND WRITE A NEW HELPER FUNCTION WHICH CONSTRUCTS THE INDIVIDUAL CONTRIBUTION OF THE MPO IN A MORE ROBUST WAY FOLLOWING THE STRATEGY 
// OF THE MPO PROCEDURE FOR THE CHIRAL TERMS 

#include <string> 
#include <vector> 
#include <array> 
#include <complex>

const extern std::string Id = "Id"; 
const extern std::string Sz = "Sz";
const extern std::string Sp = "Sp"; 
const extern std::string Sm = "Sm";

const extern std::array<std::string,3> pauliNot = {Sz,Sp,Sm};

// CHIRAL OBJECTS HOLD THE SITES, THE CHIRAL TERMS FOR EACH SITE AND SOME METADATA SUCH AS THE TYPE OF TRIANGLE FOR DEBUGGING PURPOSES
class nnObjects{
    
    public: 

        std::array<int,2> sites;
        std::array<std::string,2> nnTerm;
        std::complex<double> cConst; 

        nnObjects() : sites({0,0}),  
                    cConst(0)
                    {};

        nnObjects(std::array<int,2>, 
                std::array<std::string,2>, 
                std::complex<double>); 
};

inline nnObjects:: 
nnObjects(std::array<int,2> sites_, 
        std::array<std::string,2> nnTerm_, 
        std::complex<double> cConst_)
        : sites(sites_), 
          nnTerm(nnTerm_),  
          cConst(cConst_) {};




// I NOTE THAT THERE COULD BE POTENTIAL COMPATABILITY ISSUES BETWEEN THE COUPLINGS AS WE HAVE LISTED THIS AS A COMPLEX 
// AND CONVERSION FROM REAL TO COMPLEX IS NOT ALLOWED IN C++ UNLESS ITENSOR IS ABLE TO HANDLE THIS EFFECTIVELY WHICH I 
// DO NOT KNOW 

std::vector<nnObjects> nnTermsFinder(int i, int j, const std::complex<double>& J1_, const int& N_){

    //VECTOR TO HOLD THE NNTERM OBJECTS IMPORTANTLY THESE NEED TO HOLD THE DIFFERENCE BETWEEN THE OPERATORS AND 
    // THE FINAL OPERATOR TELLS US THE POTENTIAL ENDSTATES

    std::vector<nnObjects> nnTerms; 

    // LIST OF ALL POTENTIAL PAIRINGS 
    //LIST OF POTENTIAL ORDERINGS STORED IN THIS ARRAY 
    std::array<std::array<std::string,2>,3> allPairings;
    // NOTE THAT EVEN IN C++17 WE NEED THE EXTRA BRACKETS WHICH SEEMS SILLY BUT ANYWAYS THAT IS JUST HOW IT IS 
    allPairings = {{{Sz,Sz},{Sp,Sm},{Sm,Sp}}};

    // NEED TO MULTIPLY BY FACTOR OF 1/2 
    double div2 = 0.5;

    //THE PAIRINGS SHOULD BE GIVEN BY THE DIFFERENCE BETWEEN TWO SUCCESSIVE TERMS

    int ijDif = j - i; 

    // CONDITIONS FOR WHEN WE I AND J ARE NN COUPLINGS - WE WILL MAKE AMENDMENTS TO THIS AND INCLUDE THE NNN TERMS BECAUSE 
    // TO BE HONEST I CANNOT BE BOTHERED TO RECODE SOMETHING WHICH IS IN ESSENCE THE ENTIRELY SAME ARCHETECTURE

    //IF N AWAY AND SO NEAREST NEIGHBOUR TO THE RIGHT
    if (ijDif == N_){

        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J1_;

            if(a[0] != "Sz"){
                ncConst = J1_*div2;
            }
            nnTerms.push_back(nnObjects({i,j},a,ncConst));                        
        }

    }
    // NN COUPLING FROM SITE 1 TO SITE N FROM WRAPPING AROUND CYLINDER
    else if (ijDif == N_-1 && (j%N_==0)){

        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J1_;

            if(a[0] != "Sz"){
                ncConst = J1_*div2;
            }
            nnTerms.push_back(nnObjects({i,j},a,ncConst));                           
        }
    } 

    // NN COUPLING FROM SITE SITE I TO I+1 UNLESS WE ARE IN THE J%N_ == 1 FIRST SITES WHICH DO NOT COUPLING WITH THE SITES IMMEDIATELY BEFORE 
    else if (ijDif == 1 && (j%N_!=1)){

        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J1_;

            if(a[0] != "Sz"){
                ncConst = J1_*div2;
            }
            nnTerms.push_back(nnObjects({i,j},a,ncConst));                          
        }
    }

    else{
        // std::cout << "Check code again lol" << std::endl; 
        // DEFAULTED TO ZERO AND SO CANNOT CONTRIBUTE IN ANY CIRCUMSTANCE
        std::complex<double> theBadZero = 0;
        nnTerms.push_back(nnObjects({-1,-1},{"no","no"},theBadZero)); 
    }


    return nnTerms;
}


// NOW WE DEFINE THE EQUIVALENT HELPER FUNCTION FOR THE NNN OBJECTS --- TO BE DONE 

std::vector<nnObjects> nnnTermsFinder(int i, int j, const std::complex<double>& J2_, const int& N_){

    //VECTOR TO HOLD THE NNTERM OBJECTS IMPORTANTLY THESE NEED TO HOLD THE DIFFERENCE BETWEEN THE OPERATORS AND 
    // THE FINAL OPERATOR TELLS US THE POTENTIAL ENDSTATES

    std::vector<nnObjects> nnnTerms; 

    //LIST OF POTENTIAL ORDERINGS STORED IN THIS ARRAY
    std::array<std::array<std::string,2>,3> allPairings;
    // NOTE THAT EVEN IN C++17 WE NEED THE EXTRA BRACKETS WHICH SEEMS SILLY BUT ANYWAYS THAT IS JUST HOW IT IS 
    allPairings = {{{Sz,Sz},{Sp,Sm},{Sm,Sp}}};


    // NEED TO MULTIPLY BY FACTOR OF 1/2 
    double div2 = 0.5;

    //THE PAIRINGS SHOULD BE GIVEN BY THE DIFFERENCE BETWEEN TWO SUCCESSIVE TERMS

    int ijDif = j - i; 

    // CONDITIONS FOR WHEN WE I AND J ARE NN COUPLINGS - WE WILL MAKE AMENDMENTS TO THIS AND INCLUDE THE NNN TERMS BECAUSE 
    // TO BE HONEST I CANNOT BE BOTHERED TO RECODE SOMETHING WHICH IS IN ESSENCE THE ENTIRELY SAME ARCHETECTURE

    //IF N AWAY AND SO NEAREST NEXT NEIGHBOUR TO THE RIGHT
    if (ijDif == N_+1){

        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J2_;

            if(a[0] != "Sz"){
                ncConst = J2_*div2;
            }
            nnnTerms.push_back(nnObjects({i,j},a,ncConst));                        
        }

    }
    // NN COUPLING FROM SITE 1 TO SITE N FROM WRAPPING AROUND CYLINDER
    else if (ijDif == N_-1 && (j%N_ != 0)){

        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J2_;

            if(a[0] != "Sz"){
                ncConst = J2_*div2;
            }
            nnnTerms.push_back(nnObjects({i,j},a,ncConst));                           
        }
    } 

    // NN COUPLING FROM SITE SITE I TO I+1 UNLESS WE ARE IN THE J%N_ == 1 FIRST SITES WHICH DO NOT COUPLING WITH THE SITES IMMEDIATELY BEFORE 
    else if (ijDif == 1 && (j%N_ == 1)){

        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J2_;

            if(a[0] != "Sz"){
                ncConst = J2_*div2;
            }
            nnnTerms.push_back(nnObjects({i,j},a,ncConst));                          
        }
    }

    else if (ijDif == 2*N_-1){
        // WRAPPED AROUND AND CONNECTING SITE 1 WITH SITE ABOVE TO THE RIGHT WHICH IS THE 2N SITE AND SO THE DIFFERENCE SHOULD BE 2N-1 
        
        for (std::array<std::string,2> a : allPairings){

            std::complex<double> ncConst = J2_;

            if(a[0] != "Sz"){
                ncConst = J2_*div2;
            }
            nnnTerms.push_back(nnObjects({i,j},a,ncConst));                          
        }

    }
    else{
        // std::cout << "Check code again lol" << std::endl; 
        // DEFAULTED TO ZERO AND SO CANNOT CONTRIBUTE IN ANY CIRCUMSTANCE
        std::complex<double> theBadZero = 0;
        nnnTerms.push_back(nnObjects({-1,-1},{"no","no"},theBadZero)); 
    }


    return nnnTerms; 
}



class nnCont{
    
    public: 

        int row; 
        int col;
        std::string opType;
        std::complex<double> cConst;  
        nnCont() : row(0),
                    col(0),  
                    cConst(0)
                    {};
        nnCont(int, int, std::string, std::complex<double>); 
};

inline nnCont:: 
nnCont(int row_, int col_, std::string opType_,std::complex<double> cConst_)
        : row(row_), 
          col(col_), 
          opType(opType_), 
          cConst(cConst_) {};



// THIS MATCHES THE FIRST OPERATOR WITH A STANDARD SET OF OPERATORS AND GIVES THE STARTING SECTOR
// SECTOR 1 MEANS THAT WE START WITH AN SZ OPERATOR, SECTOR 2 MEANS THAT WE STARTED WITH AN SM, 
// SECTOR 3 MEANS THAT WE START WITH AN SP OPERATOR 
int checkOperator(std::string T, std::array<std::string,3> S){

     int counter=1;

     while(T != S[counter-1]){
         counter +=1;
     }

     return counter; 

}


// n REPRESENTS THE POSITION IN THE CIRCUMFERENCE AS THIS WILL RETURN A VECTOR CONTAINING ALL THE RELEVANT 
// CHIRAL TERMS FOR THE MPO FOR SITE n 
std::vector<nnCont> nnHelperFunction(int const& n, const std::complex<double>& J1_, const int& N_)
{
    std::vector<nnCont> siteNNTerms;
    std::vector<nnObjects> bigOlVector; 
    const std::array<std::string,3> horizontalSectorOp{Sz,Sm,Sp};
    const std::array<std::string,3> verticalSectorOp{Sz,Sp,Sm};
    int sGap = 3*N_ -1;

    int nmod = n%N_;
    if(nmod==0){
        nmod = 2*N_;
    }
    else if(nmod==1){
        nmod = N_+1;
    }
    else{
        nmod += N_;
    }
    
    for (int i = 1; i <= 2*N_ ; i++){

        std::vector<nnObjects> nnTerms = nnTermsFinder(i,nmod,J1_,N_);
        
        for(nnObjects i : nnTerms){

            if(i.sites[0] != -1){
                bigOlVector.push_back(i);
            }

        }

    }

    // NOW THE BIGOLVECTOR IS POPULATED WITH THE CORRECT NN TERMS FOR A GIVEN ORDER - WE HAD SOMETHING DIFFERENT IN THE OTHER FILE 
    // WHICH COULD POTENTIALLY BE THE BUG WITH THE OTHER CODE. WE CHECK WHAT THE VALUE OF THE FIRST OPERATOR IS IN THE CHIRALOBJECT
    // IMPORTANTLY WE NOTE THAT ALL THE SITES ARE TRANSLATED BY N_

    for (nnObjects T : bigOlVector){

        int startingSector = (sGap*checkOperator(T.nnTerm[0], horizontalSectorOp)) + 2 ; 
        std::string middleOperator = T.nnTerm[1];

        int ijDif = T.sites[1] - T.sites[0]; 

        int row = startingSector - (ijDif -1); 
        int col = 1;

        siteNNTerms.push_back(nnCont(row, col, middleOperator, T.cConst));

        //FINISHING STATE FOR THE CHIRAL TERMS IS DEFINED IN THE DRESSEDNEARESTNEIGHBOUR.H FILE 


    }

    return siteNNTerms;
}



// n REPRESENTS THE POSITION IN THE CIRCUMFERENCE AS THIS WILL RETURN A VECTOR CONTAINING ALL THE RELEVANT 
// CHIRAL TERMS FOR THE MPO FOR SITE n 
std::vector<nnCont> nnnHelperFunction(int const& n, const std::complex<double>& J2_, const int& N_)
{
    std::vector<nnCont> siteNNNTerms;
    std::vector<nnObjects> bigOlVector; 
    const std::array<std::string,3> horizontalSectorOp{Sz,Sm,Sp};
    const std::array<std::string,3> verticalSectorOp{Sz,Sp,Sm};
    int sGap = 3*N_-1;

    int nmod = n%N_;
    if(nmod==0){
        nmod = 2*N_;
    }
    else if(nmod==1){
        nmod = N_+1;
    }
    else{
        nmod += N_;
    }
    
    for (int i = 1; i <= 2*N_ ; i++){

        std::vector<nnObjects> nnnTerms = nnnTermsFinder(i,nmod,J2_,N_);
        
        for(nnObjects i : nnnTerms){

            if(i.sites[0] != -1){
                bigOlVector.push_back(i);
            }

        }

    }

    // NOW THE BIGOLVECTOR IS POPULATED WITH THE CORRECT NN TERMS FOR A GIVEN ORDER
    // IMPORTANTLY WE NOTE THAT ALL THE SITES ARE TRANSLATED BY N_

    for (nnObjects T : bigOlVector){


        int startingSector = (sGap*checkOperator(T.nnTerm[0], horizontalSectorOp)) + 2 ; 
        std::string middleOperator = T.nnTerm[1];

        // BECAUSE WE HAVE TAKEN THE MODULUS WE KNOW THAT THE DIFFERENCE IS ALWAYS GOING TO BE 
        // WITHIN 2N-1 AND CORRECT - THE DIFFERENCE DEFINES WHERE THE OPERATOR NEEDS TO BE PLACED 
        // IN THE W MATRIX 
        int ijDif = T.sites[1] - T.sites[0]; 

        int row = startingSector - (ijDif -1); 
        int col = 1;

        siteNNNTerms.push_back(nnCont(row, col, middleOperator, T.cConst));


    }

    return siteNNNTerms;
}


/// WE WILL HAVE A FUNCTION WHICH RETURNS ALL THE DIFFERENT CONTRIBUTIONS FROM IDENTITY OPERATORS AND STARTING OPERATORS WHICH HAVE NOT ALREADY BEEN DONE 
// THIS SHOULD BE A FUNCTION OF THE N_, AND THE GEOMETRY OF THE PROBLEM THAT WE ARE TACKLING TYPICALLY FOR THE CHIRAL PERTURBATION WE HAVE SOMETHING WHICH IS 
// 9N_ +2 DIMENSION SQUARE MATRIX WHICH MEANS THAT THE DIMSEP= 3N AND THIS FUNCTION SHOULD TAKE THE SEPARATION AS 3N, IF CHIRAL WE INCLUDE THE END STATES FOR THE CHIRAL 
// TERMS IN THE CHIRAL SECTOR 

std::vector<nnCont> idHelperFunction(const int& N_)
{
    std::vector<nnCont> s; 
    int sGap = 3*N_-1; 
    int mpoDim = 9*N_-1;
    int szStart = sGap + 2; 
    int smStart = sGap*2 +2; 
    int spStart = sGap*3 + 2;

    //DECLARING THE INITIAL AND FINAL STATES IN THE MPO MATRIX 
    s.push_back(nnCont(1,1,Id,1));
    s.push_back(nnCont(2,2,Id,1));
    
    //DECLARE THE INITIAL STATES FOR THE MPO 
    s.push_back(nnCont(2,szStart,Sz,1));
    s.push_back(nnCont(2,smStart,Sm,1));
    s.push_back(nnCont(2,spStart,Sp,1));

    //CHIRAL STATES FINISHING TERMS
    s.push_back(nnCont(3,1,Sz,1));
    s.push_back(nnCont(3+sGap,1,Sp,1));
    s.push_back(nnCont(3+2*sGap,1,Sm,1));

    //DECLARE ALL THE IDENTITIES ON THE OFF-DIAGONAL
    //WE NEED TO CHECK THROUGH THIS AGAIN TO ENSURE IN CORRECT PLACE - THIS SHOULD BE THE FINAL 
    int counter = 1;
    int bounds = 0;
    while (bounds <= mpoDim){

        int col = 2 + counter;
        int row = col+1; 
        
        // WE ALSO HAVE SOME MORE CONDITIONS WHICH ARE WHEN WE START A NEW SECTOR 
        // THIS MEANS THAT WHEN WE HAVE OUR FINAL CHIRAL TERM WE ALSO NEED TO NOT HAVE AN 
        // IDENTITY 
        // WE ALSO NEED A SPLIT BETWEEN THE CHIRAL AND NON CHIRAL TERMS
        if( row == mpoDim - (2*N_-2) || 
            row == mpoDim - sGap - (2*N_-2) || 
            row == mpoDim - 2*sGap - (2*N_-2) ||
            row == 3 + sGap ||
            row == 3 + 2*sGap
            )

        {
            counter +=1;
        }
        else{
            s.push_back(nnCont(row,col,Id,1));
            counter +=1;
        }

        bounds = row+1;
    }

    return s;

}




#endif 


