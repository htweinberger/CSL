
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
********************************************************************
////////////////////////////////////////////////////////////////////
///////////////////////////// READ ME //////////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CSL master's project code using the ITensor C++ library and iDMRG for
a J1-J2 Heisenberg interaction with chiral perturbation Hamiltonian 
over all plaquettes. 

Need to link the iDMRG.h file using MakeFile which is included in the 
GitHub. Specific flags may depend on usage. Note that the iDMRG.h code 
is ITENSOR code thanks to "The ITensor Software Library for Tensor Network 
Calculations", Matthew Fishman, Steven R. White, E. Miles Stoudenmire, 
arxiv:2007.14822. 

USER NOTE : The MakeFiles require specification of the exact location for 
both the idmrg.h header file and itensor library. This requires 
updating each MakeFile. There should only be  one MakeFile per 
subdirectory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
********************************************************************
////////////////////////////////////////////////////////////////////
/////////////////////////// runAuto.cc /////////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runAuto.cc file takes in one command line argument to open a .txt
parameter file. All of the parameters in the file need to be specified 
and an example of this is included on the GitHub as the ExampleSimulation.txt. 
This file specifes both the initialisation parameters for the MPO object 
as well as the iDMRG sweep object. This has to be in a folder with the 
same name i.e. /SavedFiles/ExampleSimulations/ExampleSimulations.txt

PULL : 15/09/2022 Boolean operators to specify a random initial setup 
and infinite system are now parameters in the command line argument 
which are {yes,no}

PULL : 15/09/2022 the randomGen.h file now takes in c string arrays 
for Fisher-Yates algorithm to randomly assign {"Up","Down"} to the 
site in the 1D lattice. SU(2) symmetry means that the overall spin-1/2
projected on the z-axis must be zero.

Dependencies .h files : 
 -- "MPOconstruction.h" contains the files for the construction of the MPO 
 -- "dmrgCustomObserver.h" to allow for printing of intermediate values 
 -- "randomGen.h" in the initialisation of the lattice 

NOTE: the savefolder has to be already exist before the iDMRG
simulation. Saves the ground state solution in this directory for 
subsequent analysis. 

Outputs in .../SavedFiles/YOURSIMULATION directory: 
  -- correlationsSz.txt (for analysis) for <psi|Sz_1 Sz_j|psi> 
  -- [filename]Intermediate (for analysis) 
  -- psiFinalFile (for reading) 
  -- psiInitialFile (for reading) 
  -- sitesInitialFile (for reading) 
  -- sitesFinalFile (for reading) 
  -- initialStateFile.txt (for analysis) 

********************************************************************
////////////////////////////////////////////////////////////////////
//////////////////////// MPOConstruction.h//////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************

This file is based on the example file on the iTensor github however
is adapted for our CSL use case. Importantly we have a 2D system and 
therefore have added complexity keeping track of the elements in the 
infinite MPO representation. This is done using helper functions, many
of which have similar functionality and therefore should ideally use 
inheritance which is not implemented. 

Each header file in MPOConstuction.h relates to different helper functions
 -- dressedChiralFinder.h (Chiral terms) 
 -- dressedNNObject.h (NN terms, NNN terms and identities)

Each of these deals with the 2D geometry of the system and wrapping and 
contains functions which return a structure which has:
 -- row, col (of the MPO matrix) for a given site in the unit cell as well 
as the operator ("Sz","Sp","Sm","Id") and the coupling value. 

These objects are used by MPOConstruction.h to construct the correct MPO 
matrix.

********************************************************************
////////////////////////////////////////////////////////////////////
//////////////////////////// MakeFile //////////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************

MakeFile should be in the same directory as the RunAuto.cc and .h files 
for linking. This links to the ITensor Library directory (stored locally) 
and the directory containing the iDMRG source files. Flags should be 
correctly set for compiling the code. Once compiled use the executible 
file with the parameter.txt file as the command line argument. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
********************************************************************
////////////////////////////////////////////////////////////////////
////////////////////////// Analysis ////////////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Read in the psiFinalFile and sitesFinalFile. Again each analysis file 
needs to be suitably compiled and linked correctly using a MakeFile.
Each needs to be compiled once and then takes the parameter.txt file  
as a command line argument. Looks for files in "...SavedFiles/SAVEFILENAME/..."

There are three main separate streams of analysis: 

********************************************************************
////////////////////////////////////////////////////////////////////
////////////////////////// Correlations ////////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************

 -- Correlations (Sz_iSz_j,Sp_iSm_j,Sp_jSm_i) for i-j on lattice and 
calculates pairwise correlations as a std::vector. This automatically
saves ".../SavedFiles/YOURSIMULATION/correlations.txt" containing the different 
correlations. Requires user input of exact path in correlations.cc.
correlations.h contains pairwiseCorrelations which returns a std::vector
containing a class with correlation data. 
	- correlations.cc depends on correlations.h 


Takes path to ".../SavedFiles/YOURSIMULATION" as a variable defined in the 
code and therefore needs user edit. 

********************************************************************
////////////////////////////////////////////////////////////////////
//////////////////////// Transfer Matrix ///////////////////////////
////////////////////////////////////////////////////////////////////
********************************************************************

-- transferMatrix.cc calculates the leading and sub-dominant eigenvalues 
for twice the unit cell all in the left orthogonal gauge. This returns 
the left and right, dominant and sub-dominant eigenvalues of the 
transfer matrix in the right orthogonal gauge. arnoldiHelpers.h contains 
a map  T: V_{ij} --> V_{ij} for the arnoldi iterative method and therefore 
efficiently calculates the transfer matrix contracted onto a vector. 
	- transferMatrix.cc depends on arnoldoHelpers.h 

Potential issues with regauging and so do not trust this. Moreover the 
ITensor arnoldi iterative code needs improving. Discrepancy with 
subleading eigenvalue and decay in Sz correlations.

Takes path to ".../SavedFiles/YOURSIMULATION" as a variable defined in the code 
and therefore needs user edit. 

********************************************************************
////////////////////////////////////////////////////////////////////
///////////////// K Resolved Entanglement Spectrum /////////////////
////////////////////////////////////////////////////////////////////
********************************************************************


********************************************************************
**************** PULL : RTENSOR FOLDER 16/09/2022 ******************
********************************************************************

I had not linked the files which compute the leading eigenvector of the 
permuted transfer matrix which is needed to calculate the momentum resolved 
entanglement spectrum. This is now on the GitHub  this is contained within 
the directory RMatrix which computes and saves the left stationary point 
to the permuted transfer matrix and this is a representation of the
rotation operator on the bond index level. 

Again we have a specific MakeFile to compile and link this .cc file and its 
dependencies. 

completeEntanglement.h is a header file containing all important functions and 
classes for the 

This saves a file in the directory .../SavedFiles/YOURSIMULATION/RTensorFile
containing a ITensor object which can be read by kEigensolver 

We also extract Dense representation of the QNMPS as a .txt file saving the 
values of the RTensor matrix at .../SavedFiles/YOURSIMULATION/RMatrix.txt

We also calculate the overlap of <psi|R|psi> and save this in the file 
.../SavedFiles/YOURSIMULATION/rotationQN.txt. This contains normalisation 
information as well as some other metrics to quantify how good the 
R matrix is. 

********************************************************************
************************ kEigensolver.cc ****************************
********************************************************************

-- kEigensolver.cc takes in the /RTensorFile which is at ".../SavedFiles/
YOURSIMULATION/RTensorfile and takes the parameter.txt file as a command line 
argument to specify .../SavedFiles/YOURSIMULATION directory. 

We have an eigenvalue equation R|P> = exp(iky)|P> where |p> are the eigenvectors
of the transfer matrix since [R,H]=0, we can simultaneously diagonalise the 
eigenvectors of the Hamiltonian and Rotation operator therefore R is a good 
quantum number. 

(1) Compute the reduced density operator using the Schmidt values at the cut 
 - This is done by rhoTildeConstruct which takes in the MPS 

(2) Diagonalise the density operator which is closed at the cut and find its 
eigenvectors
 - exactDiagonalisation takes in an the density operator which is closed 
on the LHS at the cut and returns its eigenvectors and eigenvalues as a vector
 - The eigenvalues are the singular values of the cut and since this is a 
global topological feature, the matrix of Schmidt values is the LHS eigenvector 
of the Transfer matrix with unit eigenvalue. 
 
(3) momentumResolvedEigensolver takes rotationVector
 - Normalises each column of R since we expect R to be unitary R. 
 - Compute <P|R|P> where |P> is an eigenvector of the R 
 - Outputs (Eigenvalue), (Real{<P|R|P>},Im{<P|R|P>}) \n 
in "SavedFiles/YOURSIMULATION/momentumResolvedEigenvalues.txt"

This is the important file containing the momentumResovledEigenvalues and
can use desired output. 

Takes path .../SavedFiles/YOURSIMULATION as a variable defined in the code 
and therefore needs user edit. 

********************************************************************
////////////////////////////////////////////////////////////////////
//////////////////completeEntanglementTMUpdated.h //////////////////
////////////////////////////////////////////////////////////////////
********************************************************************

This header file is quite messy and needs a lot of cleaning. A lot of consists
of unnecessary debugging code. 

-- completeEntanglementTMUpdated.h is a header file containing all important 
functions and classes related to the kEigensover.cc 

(-) translateIndices function which takes in &MPS and number of sites translated 
T and returns an index mapping F : i --> i' for the real index on a unit cell 
which is used for the contraction of the permuted transfer matrix. 
Returns std::vector containing (i,i'). 
(-) checkRightOrthogonality(const &MPS) 
(-) checkLeftOrthogonality(const &MPS)                
(-) orthogonalRight mapping : MPS --> list(ITensor) in right orthogonal gauge  
(-) rightMap for arnoldi Iteration  T: V_{ij} --> V_{ij} to find the right 
eigenvector of a vector of ITensors which make up the unit cell
(-) rhoTildeConstuct takes in psi in the right orthogonal gauge and constructs 
the reduced density operator on the bond level and is closed on the LHS of the cut. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






