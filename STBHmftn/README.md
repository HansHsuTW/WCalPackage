# STBHmftn
This program can create a finite slab tight-binding model with the (hkl) surface from the bulk tight-binding model. 

## Getting Started 
1) Prepare the input file (input.txt) for the program "SlabCell.m"
```
SlabCell

wcal.nL     = 40                  % number of layers stacking along [hkl]
wcal.isCut  = 0                   % isCut = 1 to prepare the desired termination for the finite slab
wcal.LB     = [-0.1 0.5]          % if isCut = 1, LB is for the lower and upper boundaries. 
wcal.hkl    = [0 0 1]             % (hkl)
wcal.isConv = 0                   % For the case when the conventional unit vectors do not coincde with the Cartesian xyz
wcal.ref    = ftn58sparse.mat     % name for the bulk TB model

endSlabCell
```
2) Run "SlabCell.m" to get the necessary infomation to construct the slab TB model from the bulk TB model
```
>> SlabCell.m
```
After running this, two structure files "Unit.vasp" and "Slab.vasp" are generated. You can check "Slab.vasp" to see whether this is the structure you want. 

3) Run "STBHmftn.m" to get the slab TB model. 
```
>> STBHmftn.m
```
