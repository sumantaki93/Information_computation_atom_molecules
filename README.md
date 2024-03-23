This project demonstrates a basic idea of how the Gaussian-type orbitals (GTO) are algorithmically implemented inside any Hartee-Fock (HF) or Density Functional theory (DFT) code. In brief, this code reconstructs atomic/molecular orbitals in both position and momentum space to compute all-electron density and its analytical gradient by accessing coefficient matrices, atomic positions, GTO details, and the atomic numbers for a given atomic or molecular systems, from any standard DFT/HF package. Since the computational expenses gradually increase as the system size grows, that bottleneck has been tried to overcome by implementing a multiprocessing environment. So, this is a completely parallel FORTRAN 90 code. Finally, by performing a couple of numerical integrations on the parallel computed electron density one can obtain different information values (Shannon, fisher, etc.)     

As a clear instruction, this code provides no further scripts to generate the input data set. One may integrate any standard script to provide the following input files for simplicity. Here 'cmatup.dat' and/or 'cmatdown.dat' are coefficient matrices (from the DFT/HF calculations), 'geo.dat' is the atomic positions, 'lmn.dat' is the information about atomic orbitals for a given basis set, 'alpha beta method.dat' is the file which contains the number of up or down spins present in the system, 'noa npg noc norb.dat' file stands for the number of atoms, number of primitive gaussian, number of contructed gaussian and the number of occupied orbitals. Furthermore, 'grid point.dat' represents the number of points that need to be framed along r, theta, and phi directions (cause the scheme has been implemented in the spherical polar grid). 'npg1.dat' contains the frequency of the atomic orbitals.  'exponent.dat' and 'precoef.dat' are the GTO data available from the choice of GTO-basis (One may follow the Basis set exchange website). In the file 'CPU.txt', you should specify both the actual number of processors utilized during the computation and the maximum number of processors that can be accommodated for a given task. All the data you need to provide into the 'main_cal.f90' routine (inside the source directory). For more information follow https://doi.org/10.1016/j.comptc.2020.112801                    
 


One may use this script to execute information and complexity for atoms or molecules (As an example {B3LYP,6-311G**} basis and functional are included for a test run)  
  
Input files are as follows  

1) cmatup.dat** (AO/MO, DFT/HF)  
2) cmatdw.dat** (AO/MO, DFT/HF)  
3) geo.dat  
4) lmn.dat  
5) alpha beta method.dat  
6) noa npg noc norb.dat  
7) grid point.dat  
8) atomic number.dat  
9) npg1.dat  
10) exponent.dat  
11) precoef.dat  
12) CPU.text  
  
  
  
**For unrestricted calculations both files (up and down) are necessary, but for restricted cases, one can use any one of them (up or down).    
 
Install Open MP, then provide the Num of CPUs and the maximum number of CPUs   

Commands to run  
-------------------------
gfortran -O3 -fopenmp -o emd atom parallel psix.f90 basisfn.f90 cinorm.f90 factor.f90 gauss legendre.f90 main cal.f90  
./emd atom parallel



 




























 





