#include "itensor/all.h"
#include <iostream>


using namespace itensor; 

// IN THIS CODE WE CALCULATE THE TRANSFER MATRIX BY THE ARNOLDI ITERATIVE SOLVER WHICH 
// FINDS DOMINANT AND SUB-DOMINANT EIGENVALUES. WE START BY CONTRACTING TRANSFER MATRIX 
// ONTO SOME RANDOM ITENSOR X WHICH CONSERVES SZ. 

class leftMap
  {
  MPS const& A_;
  mutable long size_;

  public:
  int NuC_;


  leftMap(MPS const& A)
    : A_(A)

    {
    size_ = 1;

    NuC_ = length(A_)/2;

    ITensor sdag = dag(A_(1)); 

    Index index1 = inds(A_(1))[0];
    Index index2 = inds(A_(1))[1];
    Index index3 = inds(A_(1))[2];

    Index index1dag = inds(prime(sdag))[0];
    Index index2dag = inds(prime(sdag))[1];
    Index index3dag = inds(prime(sdag))[2];

    Index rIndex = rightLinkIndex(A_,NuC_);
    Index rIndexdag = rightLinkIndex(dag(prime(A_,"Link")),NuC_);

    auto [E,indexCom1] = combiner(rIndex, rIndexdag);
    auto [C,indexCom] = combiner(index3,index3dag);


    size_ = dim(indexCom);

    }
    

  void
  product(ITensor const& x, ITensor& b) const
    {

    ITensor sdag = dag(A_(1)); 

    Index index1 = inds(A_(1))[0];
    Index index2 = inds(A_(1))[1];
    Index index3 = inds(A_(1))[2];

    Index index1dag = inds(prime(sdag))[0];
    Index index2dag = inds(prime(sdag))[1];
    Index index3dag = inds(prime(sdag))[2];

    Index rIndex = rightLinkIndex(A_,NuC_);
    Index rIndexdag = rightLinkIndex(dag(prime(A_,"Link")),NuC_);

    auto [E,indexCom1] = combiner(rIndex, rIndexdag);
    auto [C,indexCom] = combiner(index3,index3dag);

    /// THIS IS IN ESSENCE THE CODE THAT WE WILL BE IMPLEMENTING

    Index indexX;

    for (Index i : inds(x)){
        indexX = i;
    }

    b = x * delta(dag(indexX),dag(indexCom));

    b = C*b; 
    b *= A_(1);

    b *= dag(prime(A_(1),"Link"));

    for (int i = 2; i <=  2*NuC_; i++){
        b *= A_(i);
        b *= dag(prime(A_(i),"Link"));
    }

    b*= dag(C); 
    b*=delta(indexCom,indexX);

    }

  long
  size() const
    {
    return size_;
    }

};


class rightMap
  {
  MPS const& A_;
  mutable long size_;

  public:
  int NuC_;


  rightMap(MPS const& A)
    : A_(A)
 

    {
    size_ = 1;
    NuC_ = length(A_)/2; 

    ITensor sdag = dag(A_(1)); 

    Index index1 = inds(A_(2))[0];
    Index index2 = inds(A_(1))[1];
    Index index3 = inds(A_(1))[2];

    Index index1dag = inds(prime(sdag))[0];
    Index index2dag = inds(prime(sdag))[1];
    Index index3dag = inds(prime(sdag))[2];

    auto [C,indexCom] = combiner(index3,index3dag);


    size_ = dim(indexCom);

    }
    

  void
  product(ITensor const& x, ITensor& b) const
    {

    ITensor sdag = dag(A_(1)); 

    Index index1 = inds(A_(1))[0];
    Index index2 = inds(A_(1))[1];
    Index index3 = inds(A_(1))[2];

    Index index1dag = inds(prime(sdag))[0];
    Index index2dag = inds(prime(sdag))[1];
    Index index3dag = inds(prime(sdag))[2];

    auto [C,indexCom] = combiner(index3,index3dag);

    /// THIS IS IN ESSENCE THE CODE THAT WE WILL BE IMPLEMENTING

   Index indexX;

    for (Index i : inds(x)){
        indexX = i;
    }


    b = x * delta(dag(indexX),indexCom);


    b = dag(C)*b; 

    b *= A_(2*NuC_);

    for (int i = 1; i <  2*NuC_; i++){
        b *= A_(2*NuC_-i);
    }

    b *= dag(prime(A_(2*NuC_),"Link"));

    for (int i = 1; i <  2*NuC_; i++){
        b *= dag(prime(A_(2*NuC_-i),"Link"));
    }

    b*= C; 
    
    b*=delta(dag(indexCom),indexX);
    }

  long
  size() const
    {
    return size_;
    }

};


ITensor eigenCheck(ITensor x, MPS const& psi){

    int NuC = length(psi)/2;

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

    /// THIS IS IN ESSENCE THE CODE THAT WE WILL BE IMPLEMENTING

    Index indexX;

    for (Index i : inds(x)){
        indexX = i;
    }

    x *= delta(dag(indexX),dag(indexCom));

    ITensor b = C*x; 

    b *= psi(1);
    b *= dag(prime(psi(1),"Link"));

    for (int i = 2; i <= 2*NuC; i++){
        
        b *= psi(i);
        b *= dag(prime(psi(i),"Link"));
    }

    b *= dag(C);     
    b *= delta(indexCom,indexX);

    return b; 
}