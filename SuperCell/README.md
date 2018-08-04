# SuperCell
The program can generate an arbitrary supercell from a TB model with new unit vectors constructed from the linear combination of the original unit vecotrs. 

## Getting Started 
1) Prepare the input file (input.txt) for the program "SuperCell.m"
```
SuperCell

wcal.vec1   = [2 0 0]                 % Transfromation matrix 
wcal.vec2   = [0 2 0]
wcal.vec3   = [0 0 1]
wcal.file   = ftn58sparse_cut_p.mat   % Input TB model 

endSuperCell
```
2) Run "SuperCell.m" to get the necessary infomation to construct the supercell TB model from the original TB model
```
>> SuperCell.m
```
After running this, a structure files "Unit.vasp" is generated. You can check "Unit.vasp" to see whether this is the supercell you want. 

3) Run "SuperTBHmftn.m" to get the supercell TB model. 
```
>> SuperTBHmftn.m
```
4) Check the band structure by running "band_ftn.m"
```
>> band_ftn.m
```

5) Change the "InFile" in "band_ftn.m" to "ftn58sparse_cut_cdw.mat" to compare the band structure w/ CDW and w/o CDW (supercell constructed from the prinstine monolayer TiSe2 ).
```
>> band_ftn.m
```
