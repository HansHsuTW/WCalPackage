# Unfold
A sub-program to calculate unfolding band structure.

# Getting Started 
run the following script on Matlab to get unfolding wavefunction weight at each k points,
```
>> Bandunfold.m
```

run "Plot_unfold.m" to plot the wavefunction weight. 
```
>> Plot_unfold.m
```

Change the "InFile" in "Bandunfold.m" from "ftn58sparse_2x2_pristine.mat" to "ftn58sparse_cut_cdw.mat" for comparing the unfolding band structure w/o CDW and w/ CDW (supercell constructed from the prinstine monolayer TiSe2 ), and run the above steps.
```
>> Bandunfold.m
```
```
>> Plot_unfold.m
```
