//
// CHIRAL HELPER FUNCTION WHICH TAKES IN THE POSITIONS OF THE OPERATORS
// AS THREE INTERNAL INDICES AND THEN ASSIGNS A SIGN DEPENDING ON WHETHER 
// THESE AGREE WITH SOME FIDUCIAL ORDERING WHICH SHOULD BE COMPLETELY ANTISYMMETRIC 
//

#ifndef __DRESSEDCHIRALFINDER_H 
#define __DRESSEDCHIRALFINDER_H 

#include <string> 
#include <vector> 
#include <array> 
#include <complex>

// WE DEFINE A FUNCTION WHICH TAKES IN AN ARRAY OF STRINGS IN THE ORDERING <"Sz","Sp","Sm"> AND COMPARES 
// THIS TO THE FIDUCIAL ORDERING AND RETURNS AN OBJECT WITH THE STARTING SECTOR,


using arrayArray = std::array<std::array<std::string,3>,6>;


// CHIRAL OBJECTS HOLD THE SITES, THE CHIRAL TERMS FOR EACH SITE AND SOME METADATA SUCH AS THE TYPE OF TRIANGLE FOR DEBUGGING PURPOSES
class chiralObject{
    
    public: 

        std::array<int,3> sites;
        std::array<std::string,3> chiralTerm;
        std::string triangle;
        std::complex<double> cConst; 

        chiralObject() : sites({0,0,0}),  
                    cConst(0)
                    {};

        chiralObject(std::array<int,3>, 
                std::array<std::string,3>, 
                std::string, 
                std::complex<double>); 
};

inline chiralObject:: 
chiralObject(std::array<int,3> sites_, 
        std::array<std::string,3> chiralTerm_, 
        std::string triangle_, 
        std::complex<double> cConst_)
        : sites(sites_), 
          chiralTerm(chiralTerm_), 
          triangle(triangle_), 
          cConst(cConst_) {};


// THIS IS A FUNCTION WHICH RETURNS AN ARRAY OF ARRAYS FOR ALL POTENTIAL COMBINATIONS OF THE OPERATORS FOR EACH TRIANGLE
arrayArray allTriplets()
{
    arrayArray pauliCombinatorics; 
    int counter = 0; 
    for( std::string l : pauliNot){
        for (std::string m : pauliNot){
            for(std::string n : pauliNot){
                
                if (m!=l && n!=l && n!=m){
                    pauliCombinatorics[counter][0] = l; 
                    pauliCombinatorics[counter][1] = m; 
                    pauliCombinatorics[counter][2] = n; 
                    counter +=1;
                }
            }
        }

    }

    return pauliCombinatorics;
}


// THIS TAKES IN A CONST ARRAY WHICH WE DEFINE TO BE THE POSITIVE ASSIGNMENT AND THEN WE GO THROUGH ALL POTENTIAL COMBINATIONS 
// OF THE PAULI MATRICES IN PAULI COMBINATORICS AND NOTE THE NUMBER OF SWAPS WE HAVE TO GET FROM THE ORIGINAL ARRAY TO THE NEW ARRAY 
// AND IF THE NUMBER OF SWAPS IS ODD THEN WE RETURN -1 AND IF THE NUMBER OF SWAPS IS EVEN THEN WE RETURN +1. 
// THIS SHOULD TAKE IN PAULICOMBINATORICS OBJECT AND THEN HAVE

std::complex<double> signPermutation(std::array<std::string,3> input, std::array<std::string,3> const& standard){

    int swaps = 0;

    for(int i = 0; i < 2; i++){
        
        if (standard[i] == input[i]){}

        else{

            //FIND THE POSITION OF THE ITH MEMBER OF THE GOAL ARRAY IN THE INPUT ARRAY AND SWAP POSITIONS 
            int j = i+1;
            while(standard[i] != input[j]){
                j+=1; 
            }
            std::string a = input[i];
            std::string b = input[j]; 

            input[i] = b; 
            input[j] = a;
            
            swaps +=1;

        }
    }

    if (swaps%2 ==0){
        std::complex<double> perm = +1;
        return perm; 
    }
        
    else{
        std::complex<double> perm = -1;
        return perm;
    }
}



// FUNCTION WHICH TAKES IN THE 3 SITE INDICES AND CAN TELL WHETHER WE HAVE A TRIANGLE BETWEEN THESE 
// THREE POINTS AND RETURNS A CHIRALOBJECT WHICH CONTAINS THE SITE INDICES, THE ORDER OF OPERATORS AND 
// WILL RETURN A VECTOR OF LENGTH 6 IF THERE IS A SITE AND 1 IF THERE ISNT A SITE
// WE HAVE SEPARATED FUNCTIONS WHICH COULD HAVE BEEN INTEGRATED PURELY FOR ROBUSTNESS 

// WE DEFINITELY HAVE A BUG IN THIS CODE SINCE WE HAVE DOUBLY DEFINED CHIRAL STATE CONDITIONS
// IMPORTANTLY THIS CAN BE SOLVED BY INCLUDING FLAGS FOR WHICH MIDDLE SITE WE ARE IN. THIS HAS NOW BEEN 
// DONE AND WE HAVE CONDITIONS WHICH SHOULD MAKE SURE THAT THIS IS MORE ROBUST

std::vector<chiralObject> triangular(int i, int j, int k, const std::complex<double>& JC_, const int& N_){

    std::vector<chiralObject> chiralTerms; 

    //LIST OF POTENTIAL CHIRAL ORDERINGS STORED IN THIS VECTOR
    arrayArray chiralCombinatorics = allTriplets();

    //IMPORTANTLY GIVEN THE DIFFERENCES BETWEEN I, J AND K WE KNOW EXACTLY WHICH TRIANGLES THESE REPRESENT 
    //THEREFORE TO CATEGORISE THE TRIANGLES WE SHOULD FIRST LOOK AT THE DIFFERENCE BETWEEN THE TRIANGLES AND
    //THEN RETURN SOME INTERNAL INFORMATION AS TO WHICH TRIANGLE WE ARE LOOKING AT

    //WE WANT TO ENFORCE SO(2) SPIN SYMMETRY IN THE CHIRAL TERM WHICH MEANS THAT WE WANT TO HAVE SOMETHING THAT LOOKS LIKE 
    // |_  WITH SpSzSm DEFINED TO BE POSITIVE AND SO CONVENTION WE ARE USING {Sm,Sp,Sz} --> {1,2,3} AND CYCLIC ROTATIONS 


    int ijDif = j - i; 
    int jkDif = k - j; 

    //THESE CAN ONLY TAKE A SET OF VALUES DEPENDING ON WHICH TRIANGLES WE ARE LOOKING AT 
    //WE WILL MAKE THESE DIAGRAMS ON LATEX TONIGHT AND MAKE SURE THAT WE HAVE A GOOD IDEA OF WHAT IS HAPPENING
    //WE FIRST LOOK AT THE SPECIAL CASES OF WHERE WE HAVE WRAPPING AND ENSURE THAT THESE ARE ONLY ACCESSIBLE IF 
    //THE MIDDLE STATE IS EITHER FROM THE 1ST SITE IN A SLICE OR THE NTH SITE IN A SLICE

    if (ijDif == N_-1 && jkDif == 1 && (j%N_==0)){
    // std::cout << "|_ : this is the case for the 1st site where we have wrapping" << std::endl; 
    // (i,j,k) : (Sz,Sp,Sm) +VE AND NOTE THIS HAS ACQUIRED A MINUS SIGN RELATIVE TO ABOVE 
    std::array<std::string,3> basePositive = {Sz,Sp,Sm};
    for (std::array<std::string,3> a : chiralCombinatorics)
    {
        // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
        std::complex<double> sign =  signPermutation(a,basePositive);
        std::string triangleT = "|_";
        chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                        
    }
    }
    else if (ijDif == N_-1 && jkDif == N_ && (j%N_==0)){
        // std::cout << "|- : this is the case for the Nth site where we have wrapping" << std::endl; 
        // (i,j,k) : (Sm,Sz,Sp) +VE AND NOTE THAT THIS HAS ACQUIRED A MINUS SIGN RELATIVE TO ABOVE
        std::array<std::string,3> basePositive = {Sm,Sz,Sp};
        for (std::array<std::string,3> a : chiralCombinatorics)
        {
            // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
            std::complex<double> sign =  signPermutation(a,basePositive);
            std::string triangleT = "|-";
            chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                       
        }
    }

    else if (ijDif == 1 && jkDif == N_-1 && (j%N_==1)){
    // std::cout << "-| : this is the case for the Nth site where we wrapping" << std::endl;
    // (i,j,k) : (Sm,Sp,Sz) +VE AND NOTE THAT THIS HAS ACQUIRED A RELATIVE MINUS SIGN
    std::array<std::string,3> basePositive = {Sm,Sp,Sz};
    for (std::array<std::string,3> a : chiralCombinatorics)
    {
        // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
        std::complex<double> sign =  signPermutation(a,basePositive);
        std::string triangleT = "-|";
        chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                    
    }
    }
    else if (ijDif == N_ && jkDif == N_-1 && (j%N_==1)){
        // std::cout << "_| : this is the case for 1st site where we have wrapping" << std::endl;
        // (i,j,k) : (Sp,Sz,Sm) +VE and NOTE THAT WE HAVE ACUIRED A MINUS SIGN 
        std::array<std::string,3> basePositive = {Sp,Sz,Sm};
        for (std::array<std::string,3> a : chiralCombinatorics)
        {
            // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
            std::complex<double> sign =  signPermutation(a,basePositive);
            std::string triangleT = "_|";
            chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                 
        }

    } 

    //THIS DEFINES ALL THE TRIANGLES WHICH ARE OF THE LHS SORT 
    else if (ijDif == 1 && jkDif == N_){
        // std::cout << "|_ : this is the case for all which are 2 to N around circumference sites" << std::endl; 
        // (i,j,k) : (Sp,Sz,Sm) +VE
        std::array<std::string,3> basePositive = {Sp,Sz,Sm};
        for (std::array<std::string,3> a : chiralCombinatorics)
        {
            // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
            std::complex<double> sign =  signPermutation(a,basePositive);
            std::string triangleT = "|_";
            chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                  
        }

    }
    else if (ijDif == 1 && jkDif == N_-1){
        // std::cout << "|- : this is the case for all which are 1 to N-1 around circumference sites" << std::endl; 
        // (i,j,k) : (Sz,Sm,Sp) +VE
        std::array<std::string,3> basePositive = {Sz,Sm,Sp};
        for (std::array<std::string,3> a : chiralCombinatorics)
        {
            // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
            std::complex<double> sign =  signPermutation(a,basePositive);
            std::string triangleT = "|-";
            chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                      
        }
    }
    //THIS DEFINES ALL THE TRIANGLES WHICH ARE OF THE RHS SORT 
    else if (ijDif == N_-1 && jkDif == 1){
        // std::cout << "_| : this is the case for all which are 2 to N around circumference sites" << std::endl;
        // (i,j,k) : (Sp,Sm,Sz) +VE 
        std::array<std::string,3> basePositive = {Sp,Sm,Sz};
        for (std::array<std::string,3> a : chiralCombinatorics)
        {
            // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
            std::complex<double> sign =  signPermutation(a,basePositive);
            std::string triangleT = "_|";
            chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                       
        }
    }
    else if (ijDif == N_ && jkDif == 1){
        // std::cout << "-| : this is the case for all which are 1 to N-1 around circumference sites" << std::endl;
        // (i,j,k) : (Sm,Sz,Sp) +VE 
        std::array<std::string,3> basePositive = {Sm,Sz,Sp};
        for (std::array<std::string,3> a : chiralCombinatorics)
        {
            // WE GET THE SIGN AND THEN WE ADD THIS TO A VECTOR OF ALL CHIRAL TERMS
            std::complex<double> sign =  signPermutation(a,basePositive);
            std::string triangleT = "-|";
            chiralTerms.push_back(chiralObject({i,j,k},a,triangleT,sign*JC_));                       
        }
    }

    else{
        // std::cout << "Check code again lol" << std::endl; 
        // DEFAULTED TO ZERO AND SO CANNOT CONTRIBUTE IN ANY CIRCUMSTANCE
        std::complex<double> theBadZero = 0;
        chiralTerms.push_back(chiralObject({-1,-1,-1},{"no","no","no"},"_",theBadZero)); 
    }

    return chiralTerms; 

    // AT THIS POINT WE HAVE CLASSIFIED ALL THE DIFFERENT POTENTIAL 'TRIANGLES' ON A 2D CYLINDER AND MAPPED THEM TO A 1D MPS 
    // AND WE HAVE NOTE USED MODULUS % SINCE WE FEEL THAT THIS SHOULD CAUSE AN ERROR 
}


// AT THIS POINT WE HAVE CLASSIFIED ALL OF THE DIFFERENT PLAQUETTES AND SO WE NOW HAVE THE TASK OF CREATING A FUNCTION WHICH TAKES
// A LIST OF CHIRAL TERMS AND THEN RETURNS WHERE WE SHOULD PLACE THE MIDDLE OPERATOR IN THE MPO AND THE SIGN - THIS IS EXACTLY THE 
// VEIN AS WHAT WE HAD BEFORE FOR THE NNN BUT NOW WE HAVE INTRODUCED THESE CHIRAL TERMS AND SO NOW WE WILL SEPARATE THE TWO INDEX CONTAINERS 
// SO THAT WE CAN 'TURN THE CHIRAL TERMS OFF' MORE EASILY FOR DEBUGGING

// WE START BY DEFINING A CLASS FOR THE NEW CONTAINERS WHICH WILL BE DEFINED IN THE SAME WAY AS IN THE CHIRALFINDER.H BUT AS A NEW OBJECT BECAUSE 
// THESE FILES ARE NOT GOING TO BE LINKED AND THE OBJECTS NEED TO BE WITHIN SCOPE

class chiralCont{
    
    public: 

        int row; 
        int col;
        std::string opType;
        std::complex<double> cConst;  
        chiralCont() : row(0),
                    col(0),  
                    cConst(0)
                    {};
        chiralCont(int, int, std::string, std::complex<double>); 
};

inline chiralCont:: 
chiralCont(int row_, int col_, std::string opType_,std::complex<double> cConst_)
        : row(row_), 
          col(col_), 
          opType(opType_), 
          cConst(cConst_) {};



// THIS MATCHES THE FIRST OPERATOR WITH A STANDARD SET OF OPERATORS AND GIVES THE STARTING SECTOR
// SECTOR 1 MEANS THAT WE START WITH AN SZ OPERATOR, SECTOR 2 MEANS THAT WE STARTED WITH AN SM, 
// SECTOR 3 MEANS THAT WE START WITH AN SP OPERATOR 
int checkOperatorChiral(std::string T, std::array<std::string,3> S){

     int counter=1;

     while(T != S[counter-1]){
         counter +=1;
     }

     return counter; 

}


//n REPRESENTS THE POSITION IN THE CIRCUMFERENCE AS THIS WILL RETURN A VECTOR CONTAINING ALL THE RELEVANT 
// CHIRAL TERMS FOR THE MPO FOR SITE n 
std::vector<chiralCont> chiralHelperFunction(int const& n, const std::complex<double>& JC_, const int& N_)
{
    std::vector<chiralCont> siteChiralTerms;
    std::vector<chiralObject> bigOlVector; 
    const std::array<std::string,3> horizontalSectorOp{Sz,Sm,Sp};
    const std::array<std::string,3> verticalSectorOp{Sz,Sp,Sm};
    int sGap = 3*N_-1;

    int nmod1 = n%N_;
    int nmod2; 
    if(nmod1==0){
        nmod1 = N_;
        nmod2 = nmod1 + N_;
    }
    if(nmod1==1){
        nmod1 = N_+1;
        nmod2 = nmod1 + N_;
    }
    else{
        nmod2 = nmod1 + N_;
    }
    
    for (int i = 1; i <= 2*N_ ; i++){
            for (int k = 1; k <= 2*N_ ; k++){

                std::vector<chiralObject> chiralTerms1 = triangular(i,nmod1,k,JC_,N_);
                std::vector<chiralObject> chiralTerms2 = triangular(i,nmod2,k,JC_,N_);
                
                
                for(chiralObject i : chiralTerms1){

                    if(i.sites[0] != -1){
                        bigOlVector.push_back(i);
                    }

                }

                for(chiralObject i : chiralTerms2){

                    if(i.sites[0] != -1){
                        bigOlVector.push_back(i);
                    }

                }

            }
        }

    // NOW THE BIGOLVECTOR IS POPULATED WITH THE CORRECT CHIRAL TERMS FOR A GIVEN ORDER - WE HAD SOMETHING DIFFERENT IN THE OTHER FILE 
    // WHICH COULD POTENTIALLY BE THE BUG WITH THE OTHER CODE. WE CHECK WHAT THE VALUE OF THE FIRST OPERATOR IS IN THE CHIRALOBJECT

    for (chiralObject T : bigOlVector){

        //THE CHIRAL TERM IS FULLY DEFINED BY THE CENTRAL TERM WHICH DOES THE HEAVY LIFTING - WE CAN FULLY DEFINE THE PLACEMENT IF WE KNOW THE 
        //STARTING AND FINAL SECTORS (AND CHIRAL TERM) AND THE SPACING WHICH DEFINES THE GAPS BETWEEN SUCCESSIVE OPERATOR CHAIN. WE ARE DEFINING THE 
        //STARTING SECTOR AS THE FIRST TERM IN THE NORMAL (NNN) SECTOR AND SO MULTIPLY THE SECTOR BY THE GAP BETWEEN SECTORS. THE GAP BETWEEN SECTORS IS 3N 
        //AND WE HAVE 2 STATES AT THE START. THE ROW FOR THE MIDDLE OPERATOR IS DETERMINED BY THE STARTING SECTOR - (IJDIF-1). 
        //THE COLUMN FOR THE MIDDLE OPERATOR IS DETERMINED BY THE SECTOR THE ENDING SECTOR AND THE FIRST TERM IN THIS AND IN THE CASE THAT WE HAVE 2N-1 IDENTITIES 
        //ALONG THE OFF-DIAGONAL A GAP THAT NOTHING CAN GET THROUGH AND THEN WE HAVE N-1 IDENTITIES AND THE FINAL TERM. THUS WE WANT TO SUBTRACT (2N+N-jkDIF) 

        int startingSector = (sGap*checkOperatorChiral(T.chiralTerm[0], horizontalSectorOp)) + 2 ; 
        std::string middleOperator = T.chiralTerm[1];
        int endingSector = (sGap*checkOperatorChiral(T.chiralTerm[2], verticalSectorOp)) + 2;

        int ijDif = T.sites[1] - T.sites[0]; 
        int jkDif = T.sites[2] - T.sites[1];

        int row = startingSector - (ijDif -1); 
        int col = endingSector - (sGap - jkDif);

        std::complex<double> imag (0,1); 

        siteChiralTerms.push_back(chiralCont(row, col, middleOperator, imag * T.cConst/2));

        //FINISHING STATE FOR THE CHIRAL TERMS IS DEFINED IN THE DRESSEDNEARESTNEIGHBOUR.H FILE 


    }

    return siteChiralTerms;
}











#endif