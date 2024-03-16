
#  One may use this script to execute information and complexity for atoms or molecules form  (As an example {B3LYP,6-311G**} basis and functional are included for test run). 
# Needed Input files:: 1) cmatup.dat (AO/MO, DFT/HF) 
#                       2) cmatdw.dat (AO/MO, DFT/HF)               !! FOR UHF CALCULATION BOTH FILES ARE NECESSARY. BUT FOR RHF  
#                                                                   !! ONE MAY USE ONLY ONE OF THEM. BUT ANOTHER FILE SHOULD BE AVAILABLE IN THE SOURCE CODE DIRECTORY FOR SMOOTH RUN.  
#
#                       3)geo.dat
#                       4)lmn.dat
#                       5)alpha beta method.dat
#                       6)noa npg noc norb.dat
#                       7)grid point.dat
#                       8)atomic number.dat
#                       9)npg1.dat
#                      10)exponent.dat
#                      11)precoef.dat
#                      12)npg1.dat                      

# For 5 to 9 you have to put explicit values. For 'noc' and 'npg' you may see 'npgnoc', which will appear in the 'Input tool' Directory.   



 
#gfortran -O3 -fopenmp -o emd atom parallel psix.f90 basisfn.f90 cinorm.f90 factor.f90 gauss legendre.f90 main cal.f90 
#./emd atom parallel

# For more information follow https://doi.org/10.1016/j.comptc.2020.112801
 




























 





