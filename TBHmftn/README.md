# TBHmftn
This script will convert the case_hr.dat generated from the wannier90 to a packed file named ftn58sparse.mat.

## Getting Strated 
1) Prepare the structure input file (st_input.txt) for your Wannier tight-binding model. An example for SnTe is given below
```
SnTe_bulk                        % Name for the model
6.38865 6.38865 6.38865          % Lattice constant  
0.0   0.5   0.5                  % Bravais lattice vectors
0.5   0.0   0.5
0.5   0.5   0.0
2                                % Number of element species
Sn 2 3 4                         % (1st type) with the orbit index
Te 2 3 4                         % (2nd type)
1 1                              % Number of atom for each species, for exampel, 1 Sn and 1 Te. 
 0.00000   0.00000   0.00000	3   % Fractional coordinates for each atom followed by the number of orbits. 
-0.50000   0.50000   0.50000	3
```
