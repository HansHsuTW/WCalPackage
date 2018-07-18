# TBHmftn
This script will convert the case_hr.dat generated from the wannier90 to a packed file named ftn58sparse.mat.

## Getting Strated 
1) Prepare the structure input file (st_input.txt) for your Wannier tight-binding model. The example for SnTe is given below
```
SnTe_bulk
6.38865 6.38865 6.38865
0.0   0.5   0.5
0.5   0.0   0.5
0.5   0.5   0.0
2
Sn 2 3 4
Te 2 3 4
1 1
 0.00000   0.00000   0.00000	3
-0.50000   0.50000   0.50000	3
```
