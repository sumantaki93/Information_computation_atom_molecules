This project demonstrates a basic idea of how the Gaussian-type orbitals (GTO) are algorithmically implemented inside any Hartee-Fock (HF) or Density Functional theory (DFT) code. In brief, this code reconstructs atomic/molecular orbitals in both position and momentum space to compute all-electron density and its analytical gradient by accessing coefficient matrices, atomic positions, GTO details, and the atomic numbers for a given atomic or molecular systems, from any standard DFT/HF package. Since the computational expenses gradually increase as the system size grows, that bottleneck has been tried to overcome by implementing a multiprocessing environment. So, this is a completely parallel FORTRAN 90 code. Finally, by performing a couple of numerical integrations on the parallel computed electron density one can obtain different information values (Shannon, fisher, etc.), various moments of densities and complexities.     

As a clear instruction, this code provides no further scripts to generate the input data set. One may integrate any standard script to provide the following input files for simplicity. Here 'cmatup.dat' and/or 'cmatdown.dat' are coefficient matrices (from the DFT/HF calculations), 'geo.dat' is the atomic positions, 'lmn.dat' is the information about atomic orbitals for a given basis set, 'alpha beta method.dat' is the file which contains the number of up or down spins present in the system, 'noa npg noc norb.dat' file stands for the number of atoms, number of primitive gaussian, number of contructed gaussian and the number of occupied orbitals. Furthermore, 'grid point.dat' represents the number of points that need to be framed along r, theta, and phi directions (cause the scheme has been implemented in the spherical polar grid). 'npg1.dat' contains the frequency of the atomic orbitals.  'exponent.dat' and 'precoef.dat' are the GTO data available from the choice of GTO-basis (one may follow the basis set exchange website). In the file 'CPU.txt', you should specify the actual number of processors utilized during the computation and the maximum number of processors accommodated for a given task. All the data you need to provide into the **'main_cal.f90'** routine (inside the source directory). For more information follow https://doi.org/10.1016/j.comptc.2020.112801                    
 


Some calculations conducted using the DFT-B3LYP method and 6-311G** basis for a couple of atoms and molecules are attached as examples. To execute information and complexity for the chemical systems follow the given instructions. 
  
Input and output files
-----------------------------
1) cmatup.dat** (AO/MO, DFT/HF)  
2) cmatdw.dat** (AO/MO, DFT/HF)  
3) geo.dat  
4) lmn.dat  
5) alpha beta method.dat (Put method=1 or 2 for spin-restricted or unrestricted calculations)
6) noa npg noc norb.dat  
7) grid point.dat  
8) atomic number.dat  
9) npg1.dat  
10) exponent.dat  
11) precoef.dat  
12) CPU.text

 Then given code will generate the two following files as the final outputs   
 1) 'information.dat'  
 2) 'complexity.dat'     
 3) 'Moments.dat'
  
  
**For unrestricted calculations both files (up and down) are necessary, but for restricted cases, one can use any one of them (up or down).    
 
   

Commands to run  
-------------------------
Install Open MP and Fortran compiler (if necessary)

For atomic systems    
gfortran -O3 -fopenmp -o emd_atom_parallel psix.f90 basisfn.f90 cinorm.f90 factor.f90 gauss legendre.f90 main cal.f90  
./emd_atom_parallel

For  molecules   
gfortran -O3 -fopenmp -o emd_mol_parallel psix.f90 basisfn.f90 cinorm.f90 factor.f90 gauss legendre.f90 main cal.f90  
./emd_mol_parallel

 




























 





