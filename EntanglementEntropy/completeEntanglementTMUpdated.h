#ifndef __ENTANGLEMENTSPECTRUM 
#define __ENTANGLEMENTSPECTRUM

// CODE GENERATES THE STRUCTURE REQUIRED TO CALCULATE THE R MATRIX WHICH IS THE FIXED POINT OF THE MIXED TRANSFER MATRIX
// WE SHALL CONSTRUCT THIS FIXED POINT BY CONSIDER THE RIGHT MAP OF THE MIXED TRANSFER MATRIX SINCE WE CAN USE THE RIGHT 
// ORTHOGONALITY CONDITION

#include <complex>
#include <array> 
#include <iostream> 
#include <vector>
#include <cmath>

using namespace itensor;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// SHOWING THAT THE TRANSLATION OPERATOR IS A GOOD QUANTUM NUMBER //


// TAKES IN A CONSTANT REFERENCE OF THE MPS AND THEN THE NUMBER THAT WE 
// WANT TO TRANSLATE THE MPO BY
// THIS RETURNS THE LIST OF INDICES WHICH WE NEED TO CHANGE 
// THIS IS BASED ON THE CONSTRUCTIUON IN VIDAL WE WILL ALSO TRY FULLY SWAPPING 
// ALL THE POSITIONS HOWEVER THIS COULD BE PROBLEMATIC AS WE ARE IN THE GAUGE

std::vector<std::array<Index,2>> translateIndices(MPS const& psi,int T){

    //GET THE LENGTH OF THE CYLINDER 
    int NuC = length(psi)/2; 
    
    // DEFINE THE STORAGE VECTOR
    std::vector<std::array<Index,2>> translationList;
    

    for (int i=1; i <= NuC; i++){

        //SETTING AS MODULUS FOR TRANSFORMATIONS GREATER THAN THE CYLINDER WIDTH
        int iprime = (T+i)%NuC;
        if ((iprime)==0){
            iprime = NuC;
        }

        //WE NOW HAVE TO DO THE RELABELLING OF THE SITE LABELS
        // DEFINE THE OLD INDEX LABEL
        Index oldIndex; 
        for(auto& index : inds(psi.A(i))){

            if(hasTags(index,"Site")){
                oldIndex = index;
            }
        }
        Index newIndex;
        for(auto& index : inds(psi.A(iprime))){

            if(hasTags(index,"Site")){
                newIndex = index;
            }
        }

        translationList.push_back({oldIndex, dag(newIndex)});

    }

        // WE SHOULD NOW HAVE THE CORRECT TAGS FOR THE OLD AND NEW INDEX

        for (int i=NuC+1; i <= 2*NuC; i++){

        // std::cout << "------- ALL INDEX : " << i << " -----------" << std::endl;
        // for(auto& i : inds(psi(i))){
        //     println(i);
        // }

        //SETTING AS MODULUS FOR TRANSFORMATIONS GREATER THAN THE CYLINDER WIDTH
        int iprime = (T+i)%NuC;
        if ((iprime)==0){
            iprime = 2*NuC;
        }
        else{
            iprime += NuC;
        }

        //WE NOW HAVE TO DO THE RELABELLING OF THE SITE LABELS
        // DEFINE THE OLD INDEX LABEL
        Index oldIndex; 
        for(auto& index : inds(psi.A(i))){

            if(hasTags(index,"Site")){
                oldIndex = index;
                // std::cout << "------- OLD INDEX : " << i << " -----------" << std::endl;
                // println(oldIndex);
            }
        }
        Index newIndex;
        for(auto& index : inds(psi.A(iprime))){

            if(hasTags(index,"Site")){
                newIndex = index;
                // std::cout << "------- NEW INDEX : " << i << " -----------" << std::endl;
                // println(dag(newIndex));
            }
        }

        // WE SHOULD NOW HAVE THE CORRECT TAGS FOR THE OLD AND NEW INDEX

        translationList.push_back({oldIndex, dag(newIndex)});
        
    }

    return translationList; 

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// FIRST WE WANT TO CHECK THAT THE RIGHT ORTHOGONALITY CONDITION IS SATISFIED FOR OUR TENSORS BEFORE WE CAN CONTINUE
// IF IN THE ORTHOGONAL GAUGE WE SHOULD GET SOMETHING WHICH CONTRACTS TO THE IDENTITY SINCE WE HAVE PRIMED THE LEFT INDICES 
// AND THE THE RIGHT INDICES ARE UNPRIMED AND THUS CONTRACTED OVER 
// THIS CORRESPONDS TO ------0--------|
//                           |        |  = IDENTITY IFF IN THE RIGHT ORTHOGONALITY GAUGE
//                     ------0--------| 
// THIS PROCESS IS INCREDIBLY COMPUTATIONALLY INTENSIVE AS IT SCALES WITH BOND DIMENSION QUARTICALLY THUS WE LIMIT TO MPS WHICH HAVE 
// LOW BOND DIMENSION

void checkRightOrthogonality(MPS const& psi){

    int NuC = length(psi); 

    if(dim(inds(psi(0))[0]) < 15){
        for (auto ni=1; ni <= NuC; ++ni){
            std::cout << "--------------- RIGHT ORTHOGONALITY OF: " << ni << "--------------" << std::endl; 
            auto pa =psi.A(ni);
            auto indc = commonIndex(pa,psi.A(ni-1));
            PrintData(pa*prime(dag(pa),indc));
        } 
    }
}

void checkLeftOrthogonality(MPS const& psi){

    int NuC = length(psi); 

    if(dim(inds(psi(0))[0]) < 15){
        for (auto ni=1; ni <= NuC; ++ni){
            std::cout << "--------------- LEFT ORTHOGONALITY OF: " << ni << "--------------" << std::endl; 
            auto pa = psi.A(ni);
            auto indc = commonIndex(pa,psi.A(ni+1));
            PrintData(pa*prime(dag(pa),indc));
        } 
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// THIS FUNCTION TAKES IN AN MPS BY REFERENCE AND SUCCESSFULLY COMPUTES SVD TO TURN INTO THE RIGHT ORTHOGONALITY CONDITION WITH THE SCHMIDT 
// VALUES ON THE RHS, THIS SHOULD BE OBVIOUSLY COMPLETELY SINGULAR. 

std::vector<ITensor> orthogonalRight(MPS const& psi){

    int NuC1 = length(psi)/2; 
    int NuC2 = NuC1 + 1; 

    std::vector<ITensor> rightTensors;

    // THIS CONTAINS THE STATE WHICH IS BEING MADE RIGHT ORTHOGONAL
    ITensor Z = psi(NuC2); 

    for (int i = 0; i <= NuC1; i++){

        int N0 = NuC1 - i; 

        ITensor Q = psi(N0) * Z; 

        Index l1 = uniqueIndex(psi(N0),Z,"Link");
        Index r2 = uniqueIndex(Z,psi(N0),"Link");
        Index s1; 
        Index s2;
        
        for (auto& i : inds(psi(N0))){
            if (hasTags(i,"Site")){
                s1 = i;
            }
        }

        for (auto& i : inds(Z)){
            if (hasTags(i,"Site")){
                s2 = i;
            }
        }

        auto [U,S,V] = svd(Q,{l1,s1});

        rightTensors.push_back(V); 

        // NOW CONSTRUCT THE NEW LEFT STATE FOR THE NEXT ITERATION. 

        Z = U * S;

    }

    rightTensors.push_back(Z);




    return rightTensors; 

   
    
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// WE NOTE THAT WE ONLY CARE ABOUT THE RIGHT EIGENVECTOR OF THE MIXED TRANSFER MATRIX SINCE WE WANT TO BE IN THE RIGHT ORTHOGONALITY 
// CONDITION HOWEVER HOW CAN WE BE IN THE RIGHT ORTHOGONALITY CONDITION IF WE NEVER MULTIPLY IN THE ORTHOGONALITY CENTRE


class rightMap
  {
  std::vector<ITensor> const& A_;
  mutable long size_;
  

  public:
  int NuC_;
  int T_;
  std::vector<std::array<Index,2>> List_;


  rightMap(std::vector<ITensor> const& A, int T, std::vector<std::array<Index,2>> list)
    : A_(A), 
      T_(T),
      List_(list)

    {   
    size_ = 1;

    NuC_ = (A_.size()-1)/2;

    Index indexLOut = uniqueIndex(A_[0],A_[1],"Link");
    Index indexLIn = uniqueIndex(dag(prime(A_[0])),dag(prime(A_[1])),"Link");
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    auto [L,indexComL] = combiner(indexLOut,indexLIn);

    size_ = dim(indexComL);

    }
    

  void
  product(ITensor const& x, ITensor& b) const
    {

    Index indexLIn = inds(A_[1])[2];
    Index indexLOut = inds(dag(prime(A_[1],"Link")))[2];
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    auto [L,indexComL] = combiner(indexLOut,indexLIn);
    auto [R,indexComR] = combiner(indexROut, indexRIn);

    // Index indexX;

    // for (Index i : inds(x)){
    //     indexX = i;
    // }
    // ONE SINGULAR INDEX 
    // b = x * delta(dag(indexX),dag(indexComR));
    // b *= R; 

    // TWO INDICES 

    b = x * delta(dag(inds(x)[0]),dag(indexROut)) * delta(dag(inds(x)[1]),dag(indexRIn));

    b *= A_[2*NuC_];
    b *= (dag(prime(A_[2*NuC_],"Link"))*delta(List_[2*NuC_-1][0],List_[2*NuC_-1][1]));


    for (int i = 1; i <  2*NuC_; i++){
        
        b *= A_[2*NuC_-i];
        b *= (dag(prime(A_[2*NuC_-i],"Link"))*delta(List_[2*NuC_-i-1][0],List_[2*NuC_-i-1][1]));

    }


    b *= delta(inds(x)[0],dag(indexLOut));
    b *= delta(inds(x)[1],dag(indexLIn));
    
    // b *= L; 

    // b*=delta(dag(indexComL),indexX);
    }

  long
  size() const
    {
    return size_;
    }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// RIGHT EIGENVECTOR OF THE TRANSFER MATRIX /////////////////////////////////////////////////////
////////////////////////////////////////////////// IN THE RIGHT ORTHOGONAL GAUGE/RIGHT ISOMETRIC ////////////////////////////////////////////////
////////////////////////////////////////////////// ACCORDING TO THE ORTHOGONALITY CONDITIONS WE  ///////////////////////////////////////////////
////////////////////////////////////////////////// SHOULD GET T_R |I> = |I>  2019 Zauner-Stauber ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class rightTransferMatrix
  {
  std::vector<ITensor> const& A_;
  mutable long size_;
  

  public:
  int NuC_;
  std::vector<std::array<Index,2>> List_;


  rightTransferMatrix(std::vector<ITensor> const& A)
    : A_(A)

    {   
    size_ = 1;

    NuC_ = (A_.size()-1)/2;

    Index indexLOut = uniqueIndex(A_[0],A_[1],"Link");
    Index indexLIn = uniqueIndex(dag(prime(A_[0])),dag(prime(A_[1])),"Link");
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    auto [L,indexComL] = combiner(indexLOut,indexLIn);

    size_ = dim(indexComL);

    }
    

  void
  product(ITensor const& x, ITensor& b) const
    {

 
    Index indexLIn = inds(A_[1])[2];
    Index indexLOut = inds(dag(prime(A_[1],"Link")))[2];
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    // auto [L,indexComL] = combiner(indexLOut,indexLIn);
    // auto [R,indexComR] = combiner(indexROut, indexRIn);

    // Index indexX;

    // for (Index i : inds(x)){
    //     indexX = i;
    // }
    // ONE SINGULAR INDEX 
    // b = x * delta(dag(indexX),dag(indexComR));
    // b *= R; 

    // TWO INDICES  
    b = x * delta(dag(inds(x)[0]),dag(indexROut)) * delta(dag(inds(x)[1]),dag(indexRIn)); 
    b *= A_[2*NuC_];
    b *= dag(prime(A_[2*NuC_],"Link"));


    for (int i = 1; i <  2*NuC_; i++){
        
        b *= A_[2*NuC_-i];
        b *= dag(prime(A_[2*NuC_-i],"Link"));
    }
    


    b *= delta(inds(x)[0],dag(indexLOut));
    b *= delta(inds(x)[1],dag(indexLIn));
    
    // b *= L; 

    // b*=delta(dag(indexComL),indexX);
    }

  long
  size() const
    {
    return size_;
    }

};

// PROJECTING OUT THE DOMINANT SUBSPACE OF THE EIGENVALUE SPACE 
// L INDICATES THE LEFT EIGENVECTOR OF THE TRANSFER MATRIX THEREFORE GOES ON THE RHS 
// R INDICATES THE RIGHT EIGENVECTOR OF THE TRANSFER MATRIX AND THEREFORE GOES ON THE LHS

class rightTransferMatrixDom
  {
  std::vector<ITensor> const& A_;
  mutable long size_;
  ITensor L_; 
  ITensor R_; 
  std::complex<double> const& lambda_;


  public:
  int NuC_;
  std::vector<std::array<Index,2>> List_;


  rightTransferMatrixDom(std::vector<ITensor> const& A, 
                                    ITensor L, 
                                    ITensor R, 
                                    std::complex<double> const& lambda)
    : A_(A), 
      L_(L),
      R_(R),
      lambda_(lambda)

    {   
    size_ = 1;

    NuC_ = (A_.size()-1)/2;

    Index indexLOut = uniqueIndex(A_[0],A_[1],"Link");
    Index indexLIn = uniqueIndex(dag(prime(A_[0])),dag(prime(A_[1])),"Link");
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    auto [left,indexComL] = combiner(indexLOut,indexLIn);

    size_ = dim(indexComL);

    }
    

  void
  product(ITensor const& x, ITensor& b) const
    {

 
    Index indexLIn = inds(A_[1])[2];
    Index indexLOut = inds(dag(prime(A_[1],"Link")))[2];
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    // auto [L,indexComL] = combiner(indexLOut,indexLIn);
    // auto [R,indexComR] = combiner(indexROut, indexRIn);

    // Index indexX;

    // for (Index i : inds(x)){
    //     indexX = i;
    // }
    // ONE SINGULAR INDEX 
    // b = x * delta(dag(indexX),dag(indexComR));
    // b *= R; 

    // TWO INDICES  
    b = x * delta(dag(inds(x)[0]),dag(indexROut)) * delta(dag(inds(x)[1]),dag(indexRIn)); 
    b *= A_[2*NuC_];
    b *= dag(prime(A_[2*NuC_],"Link"));


    for (int i = 1; i <  2*NuC_; i++){
        
        b *= A_[2*NuC_-i];
        b *= dag(prime(A_[2*NuC_-i],"Link"));
    }
    


    b *= delta(inds(x)[0],dag(indexLOut));
    b *= delta(inds(x)[1],dag(indexLIn));

    // WE FIRST WANT TO GET THE INDICES OF THE R_ AND L_ VECTORS SUCH THAT 

    ITensor c = dag(R_) * delta(dag(inds(x)[0]), inds(R_)[0]) * delta(dag(inds(x)[1]), inds(R_)[1]); 
    
    c *= x; 
    c *= lambda_ * R_ * delta(inds(x)[0],dag(inds(R_)[0])) * delta(inds(x)[1],dag(inds(R_)[1]));;

    b -= c;
    
    // b *= L; 

    // b*=delta(dag(indexComL),indexX);
    }

  long
  size() const
    {
    return size_;
    }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// LEFT EIGENVECTOR OF THE TRANSFER MATRIX /////////////////////////////////////////////////////
////////////////////////////////////////////////// IN THE RIGHT ORTHOGONAL GAUGE/RIGHT ISOMETRIC ////////////////////////////////////////////////
////////////////////////////////////////////////// ACCORDING TO THE ORTHOGONALITY CONDITIONS WE  ///////////////////////////////////////////////
////////////////////////////////////////////////// SHOULD GET <L|T_R  = <L|  2019 Zauner-Stauber ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class leftTransferMatrix
  {
  std::vector<ITensor> const& A_;
  mutable long size_;
  

  public:
  int NuC_;


  leftTransferMatrix(std::vector<ITensor> const& A)
    : A_(A)

    {   
    size_ = 1;

    NuC_ = (A_.size()-1)/2;

    Index indexLOut = uniqueIndex(A_[0],A_[1],"Link");
    Index indexLIn = uniqueIndex(dag(prime(A_[0])),dag(prime(A_[1])),"Link");
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    auto [L,indexComL] = combiner(indexLOut,indexLIn);

    size_ = dim(indexComL);

    }
    

  void
  product(ITensor const& x, ITensor& b) const
    {


    Index indexLIn = inds(A_[1])[2];
    Index indexLOut = inds(dag(prime(A_[1],"Link")))[2];
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));


    auto [L,indexComL] = combiner(indexLOut,indexLIn);
    auto [R,indexComR] = combiner(indexROut, indexRIn);

    // Index indexX;

    // for (Index i : inds(x)){
    //     indexX = i;
    // }
    // ONE SINGULAR INDEX 
    // b = x * delta(dag(indexX),dag(indexComR));
    // b *= R; 

    // TWO INDICES 

    b = x * delta(dag(inds(x)[0]),dag(indexLOut)) * delta(dag(inds(x)[1]),dag(indexLIn));

    b *= A_[1];
    b *= dag(prime(A_[1],"Link"));


    for (int i = 2; i <=  2*NuC_; i++){
        
        b *= A_[i];
        b *= dag(prime(A_[i],"Link"));
    }
    

    b *= delta(inds(x)[0],dag(indexRIn));
    b *= delta(inds(x)[1],dag(indexROut));
    
    // b *= L; 

    // b*=delta(dag(indexComL),indexX);
    }

  long
  size() const
    {
    return size_;
    }

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// CHECKS THE LEFT EIGENVECTORS AND SEES WHETHER /////////////////////////////////////////////////////
////////////////////////////////////////////////// THEY ARE TRULY A GOOD MATCH FOR THE REGULAR ////////////////////////////////////////////////
////////////////////////////////////////////////// EIGENVALUE EQUATION IN WHICH CASE WE ARE CALCULATING  ///////////////////////////////////////////////
////////////////////////////////////////////////// SHOULD GET <L|T_R  = <L|  2019 Zauner-Stauber ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// WE WANT TO TAKE IN A VECTOR WHICH HAS TWO INDICES, MULTIPLY THIS ONTO THE OPERATOR AND SUBTRACT WITH THE EIGENVECTOR MULTIPLIED BY THE EIGENVALUE
// IN THIS CASE WE ONLY NEED TO CONSIDER THE ACT OF MULTIPLYING THE VECTOR ONTO THE OPERATOR AND KEEPING THE PREVIOUS INDICES SO THAT WE CAN SUBTRACT 
// THE ORIGINAL VECTOR IN THE NORM CALCULATIONS. THIS SHOULD THEREFORE RETURN A VECTOR 

ITensor leftEigenCheck(ITensor const& x, std::vector<ITensor> const& A_){ 

    int NuC_ = (A_.size()-1)/2;

    Index indexLIn = inds(A_[1])[2];
    Index indexLOut = inds(dag(prime(A_[1],"Link")))[2];
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    // TWO INDICES 

    ITensor b = x * delta(dag(inds(x)[0]),dag(indexLIn)) * delta(dag(inds(x)[1]),dag(indexLOut));
    b *= A_[1];
    b *= dag(prime(A_[1],"Link"));

    for (int i = 2; i <=  2*NuC_; i++){
        
        b *= A_[i];
        b *= dag(prime(A_[i],"Link"));
    }
    

    b *= delta(inds(x)[0],dag(indexRIn));
    b *= delta(inds(x)[1],dag(indexROut));

    return b; 


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// CHECKS THE RIGHT EIGENVECTORS AND SEES WHETHER /////////////////////////////////////////////////////
////////////////////////////////////////////////// THEY ARE TRULY A GOOD MATCH FOR THE REGULAR ////////////////////////////////////////////////
////////////////////////////////////////////////// EIGENVALUE EQUATION IN WHICH CASE WE ARE CALCULATING  ///////////////////////////////////////////////
////////////////////////////////////////////////// SHOULD GET T_R|I>  = |I>  2019 Zauner-Stauber ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ITensor rightEigenCheck(ITensor const& x, std::vector<ITensor> const& A_){ 

    int NuC_ = (A_.size()-1)/2;

    Index indexLIn = inds(A_[1])[2];
    Index indexLOut = inds(dag(prime(A_[1],"Link")))[2];
    Index indexROut = uniqueIndex(A_[2*NuC_],A_[2*NuC_-1],"Link");
    Index indexRIn = dag(prime(indexROut));

    ITensor b = x * delta(dag(inds(x)[0]),dag(indexROut)) * delta(dag(inds(x)[1]),dag(indexRIn));
    b *= A_[2*NuC_];
    b *= dag(prime(A_[2*NuC_],"Link"));


    for (int i = 1; i <  2*NuC_; i++){
        
        b *= A_[2*NuC_-i];
        b *= dag(prime(A_[2*NuC_-i],"Link"));
    }
    
    b *= delta(inds(x)[0],dag(indexLOut));
    b *= delta(inds(x)[1],dag(indexLIn));

    return b; 

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// CHECKS THE AFTER ROTATION THE EIGENVECTORS ///////////////////////////////////////////////////
////////////////////////////////////////////////// ARE TRULY A GOOD MATCH      //////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




ITensor eigenCheck(ITensor x, MPS const& psi, int T){

    int NuC = length(psi)/2;
    std::vector<std::array<Index,2>> list = translateIndices(psi,T);

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

    Index indexX;

    x *= delta(dag(indexX),dag(indexCom));

    ITensor b = C*x; 

    b *= psi(1);
    b *= dag(prime(psi(1),"Link"))*delta(list[0][0],list[0][1]);

    // IF WE ARE DOING THE ENTIRE UNIT CELL *2 THEN UNCOMMENT THESE SETS OF LINES

    for (int i = 2; i <= 2*NuC; i++){
        
        b *= psi(i);
        b *= (dag(prime(psi(i),"Link"))*delta(list[i-1][0],list[i-1][1]));
    }

    b *= dag(C); 
    b *= delta(indexCom,indexX);

    // IF ONLY USING THE FIRST UNIT CELL THEN UNCOMMENT THIS SET OF LINES

    // for (int i = 1; i <=  NuC; i++){

    //     b *= psi(i);
    //     b *= (dag(prime(psi(i),"Link"))*delta(list[i-1][0],list[i-1][1]));
    // }

    // b *= E; 
    // b *= delta(dag(indexCom1),indexX);
    

    return b; 
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// CHECK IF THE R MATRICES WHICH HAVE BEEN CALCU /////////////////////////////////////////////////////
////////////////////////////////////////////////// ILATED DO INDEED RELATE TO THE PERMUTATION OF ////////////////////////////////////////////////
////////////////////////////////////////////////// THE SITES BY CALCULATING THE OVERLAP          ///////////////////////////////////////////////
////////////////////////////////////////////////// <R'PSI R,PSI>             2019 Zauner-Stauber ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::complex<double> checkRRep(MPS const& psi, ITensor R){

    // WE START BY GETTING THE MPS AND MULTIPLYING THE END OF THE WAVEFUNCTION 
    // IN THE RIGHT ORTHOGONALITY CONDITION AND SO ALL THE RHS WOULD CONTRACT TO THE IDENTITY

    int NuC = length(psi)/2;

    // EVERYTHING FROM 1 - N IS IN THE RIGHT ORTHOGONALITY GAUGE 
    auto oi = uniqueIndex(psi(1),psi(2),"Link");
    auto ui = uniqueIndex(psi(NuC),psi(NuC-1),"Link");

    auto indexRIn = inds(R)[0]; 
    auto indexROut = inds(R)[1];

    ITensor x = (R * delta(dag(indexRIn), prime(dag(oi)))) * prime(psi(1),oi);
    x *= delta(dag(indexROut),prime(oi));

    ITensor y = dag(prime(psi(1),"Link"));

    ITensor norm = x*y;

    
    for (int i = 2; i < NuC; i++){

        norm *= psi(i);
        norm *= (dag(prime(psi(i),"Link")));

    }

    norm *= prime(psi(NuC),ui);
    norm *= (dag(R) * delta(indexROut, prime(dag(ui))));
    norm *= delta(indexRIn,prime(ui));
    norm *= (dag(prime(psi(NuC),"Link")));

    std::complex<double> normVal = eltC(norm);

    return normVal;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// CHECK WHETHER ROTATIONS OF THE TRANSFER ////////////////////////////
//////////////////////////////////// MATRIX IS A GOOD QUANTUM NUMBER         ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::complex<double> overlapR(MPS const& psi,std::vector<std::array<Index,2>> list){

    // std::cout << "------------------------------------------" << std::endl;

    // WE START BY GETTING THE MPS AND MULTIPLYING THE END OF THE WAVEFUNCTION 
    // IN THE RIGHT ORTHOGONALITY CONDITION AND SO ALL THE RHS WOULD CONTRACT TO THE IDENTITY

    int NuC = length(psi)/2;

    // EVERYTHING FROM 1 - N IS IN THE RIGHT ORTHOGONALITY GAUGE 
    auto oi = uniqueIndex(psi(1),psi(2),"Link");
    auto ui = uniqueIndex(psi(NuC),psi(NuC-1),"Link");

    ITensor x = prime(psi(1),oi);

    ITensor y = dag(prime(psi(1),"Link"));

    ITensor norm = x*(y*delta(list[0][0],list[0][1]));
    
    for (int i = 2; i < NuC; i++){

        norm *= psi(i);
        norm *= (dag(prime(psi(i),"Link"))*delta(list[i-1][0],list[i-1][1]));

    }

    norm *= prime(psi(NuC),ui);
    norm *= (dag(prime(psi(NuC),"Link"))*delta(list[NuC-1][0],list[NuC-1][1]));

    std::complex<double> normVal = eltC(norm);

    return normVal;

}

std::complex<double> overlapR(std::vector<ITensor> const& psi,std::vector<std::array<Index,2>> list){

    // std::cout << "------------------------------------------" << std::endl;

    // WE START BY GETTING THE MPS AND MULTIPLYING THE END OF THE WAVEFUNCTION 
    // IN THE RIGHT ORTHOGONALITY CONDITION AND SO ALL THE RHS WOULD CONTRACT TO THE IDENTITY

    int NuC = (psi.size() - 1)/2;

    // EVERYTHING FROM 1 - N IS IN THE RIGHT ORTHOGONALITY GAUGE 
    auto oi = uniqueIndex(psi[1],psi[2],"Link");
    auto ui = uniqueIndex(psi[NuC],psi[NuC-1],"Link");

    ITensor x = prime(psi[1],oi);

    ITensor y = dag(prime(psi[1],"Link"));

    ITensor norm = x*(y*delta(list[0][0],list[0][1]));
    
    for (int i = 2; i < NuC; i++){

        norm *= psi[i];
        norm *= (dag(prime(psi[i],"Link"))*delta(list[i-1][0],list[i-1][1]));

    }

    norm *= prime(psi[NuC],ui);
    norm *= (dag(prime(psi[NuC],"Link"))*delta(list[NuC-1][0],list[NuC-1][1]));

    std::complex<double> normVal = eltC(norm);

    return normVal;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CHECK THE NORMALISATION OF THE TRANSFER MATRICES ///////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


std::complex<double> normalisation(MPS const& psi){

    // std::cout << "------------------------------------------" << std::endl;

    // WE START BY GETTING THE MPS AND MULTIPLYING THE END OF THE WAVEFUNCTION 
    // IN THE RIGHT ORTHOGONALITY CONDITION AND SO ALL THE RHS WOULD CONTRACT TO THE IDENTITY

    int NuC = length(psi)/2;

    // EVERYTHING FROM 1 - N IS IN THE RIGHT ORTHOGONALITY GAUGE 
    auto oi = uniqueIndex(psi(1),psi(2),"Link");
    auto ui = uniqueIndex(psi(NuC),psi(NuC-1),"Link");

    ITensor x = prime(psi(1),oi);

    ITensor y = dag(prime(psi(1),"Link"));

    ITensor norm = x*y;
    
    for (int i = 2; i < NuC; i++){

        norm *= psi(i);
        norm *= (dag(prime(psi(i),"Link")));

    }

    norm *= prime(psi(NuC),ui);
    norm *= (dag(prime(psi(NuC),"Link")));

    std::complex<double> normVal = eltC(norm);

    return normVal;

}


std::complex<double> normalisation(std::vector<ITensor> const& psi){

    // std::cout << "------------------------------------------" << std::endl;

    // WE START BY GETTING THE MPS AND MULTIPLYING THE END OF THE WAVEFUNCTION 
    // IN THE RIGHT ORTHOGONALITY CONDITION AND SO ALL THE RHS WOULD CONTRACT TO THE IDENTITY

    int NuC = (psi.size() - 1)/2;

    // EVERYTHING FROM 1 - N IS IN THE RIGHT ORTHOGONALITY GAUGE 
    auto oi = uniqueIndex(psi[1],psi[2],"Link");
    auto ui = uniqueIndex(psi[NuC],psi[NuC-1],"Link");

    ITensor x = prime(psi[1],oi);

    ITensor y = dag(prime(psi[1],"Link"));

    ITensor norm = x*y;
    
    for (int i = 2; i < NuC; i++){

        norm *= psi[i];
        norm *= (dag(prime(psi[i],"Link")));

    }

    norm *= prime(psi[NuC],ui);
    norm *= (dag(prime(psi[NuC],"Link")));

    std::complex<double> normVal = eltC(norm);

    return normVal;

}



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// COMPUTE THE RHO MATRIX WE CAN  ////////////////////////////////////
/////////////////////////////////////ADD SINGULAR VALUES IN ALONG THE //////////////////////////////////
/////////////////////////////////////SCHMIDT CUT IF WE NEED EASILY/////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



ITensor rhoTildeConstruct(std::vector<ITensor> const& psi){

    // WE EXPLICITELY COMPUTE THE DIAGONALISATION OF THE ITENSOR

    int NuC = (psi.size() - 1)/2;

    // EVERYTHING FROM 1 - N IS IN THE RIGHT ORTHOGONALITY GAUGE 
    auto oi = uniqueIndex(psi[1],psi[2],"Link");
    auto li = uniqueIndex(psi[0], psi[1],"Link");
    auto ui = uniqueIndex(psi[2*NuC],psi[2*NuC-1],"Link");

    ITensor x = prime(psi[2*NuC],ui);
    ITensor y = dag(prime(psi[2*NuC],"Link"));

    ITensor rhoTilde = x*y;
    
    for (int i = 2*NuC-1; i >= 1; i--){

        rhoTilde *= psi[i];
        rhoTilde *= (dag(prime(psi[i],"Link")));

    }

    rhoTilde *= psi[0];

    for (auto& i : inds(rhoTilde)){
        println(i);
    }

    for (auto& i : inds(prime(dag(psi[0])))){
        println(i);
    }

    rhoTilde *= prime(dag(psi[0]));

    PrintData(rhoTilde);

    return rhoTilde;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// DIAGONALISE RHO_D USING HERMITIAN ///////////////////////////////
///////////////////////////////////// SOLVER                         //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<ITensor> exactDiagonalisationH(ITensor const& rhoTilde){

    auto [Q,D] = diagHermitian(rhoTilde);

    std::vector<ITensor> diagonalisationTensors{ Q,D };

    return diagonalisationTensors;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////WE HAVE AN EIGENVALUE EQUATION FOR THIS WHICH GOES SOMETHING LIKE THIS//////////////////////////
/////////////////////////////////// R|P> = exp(iky)|P> SINCE THE EIGENVECTORS OF THE REDUCED DENSITY MATRIX///////////////////////
/////////////////////////////// SIMULTANEOUSLY DIAGONALISE THE ROTATION OPERATOR AROUND THE CYLINDER ///////////////////////////////
///////////////////////////////////// SOLVER ////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void momentumResolvedEigenSolver(ITensor R, ITensor Q, ITensor D, std::string path){
    
    std::ofstream myfile1; 
    myfile1.open(path + "/momentumResolvedEigenvalues.txt"); 

    // GET THE BLOCK SIZES OF THE MATRICES THIS NEEDS USER INPUTTTTT 

    std::vector<int> blocksizes; 
    for (int i = 1; i <= inds(R)[0].nblock(); i++){
        blocksizes.push_back(blocksize(inds(R)[0],i));
    }
    
    int startingP = 1;
    /// WE START BY NORMALISING THE COLUMNS OF R 
    for (int max : blocksizes){

        for (int i = startingP; i < max+startingP; i++){

            ITensor ri = R*dag(setElt((inds(R)[1]=i)));
            std::complex<double> normalisation = norm(ri); 
            ri /= normalisation;
            
            for (int j=startingP; j < max+startingP; j++){
                R.set(inds(R)[0]=j,inds(R)[1]=i,eltC(ri,inds(ri)[0]=j));
            }
            
        }

        startingP += max; 

    }

    for (int i = 1; i <= maxDim(R); i++){

        ITensor q1 = Q*dag(setElt((inds(Q)[0]=i)));
        std::complex<double> normalisation = norm(q1); 
        q1 /= normalisation;

        myfile1 << eltC(D,inds(D)[0]=i,inds(D)[1]=i) << ",";

        myfile1 <<   eltC(dag(q1)
                        *delta(dag(inds(dag(q1))[0]),dag(inds(R)[1]))
                        *(R 
                        *delta(inds(dag(q1))[0],dag(inds(R)[0]))
                        *q1)) << std::endl;
    }

    myfile1.close();


}



#endif


