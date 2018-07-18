# GreenSur
This program can calculate the surface states on a semi-inifite thin film by using the Green's function method. 

## Getting Started 
From the step in "STBHmftn", you will know that the maximum hopping distace for the (001) surface is 7. Thus you should consturct the finite slab model with 7x2 = 14 layers. (for the reason please referring to the method paper)
1) Prepare the input file (input.txt) for the program "GreenSur.m". 
```
GreenSur

wcal.Ni     = 13                      % number of iterations
wcal.nn     = 0.0008                  % broading facor in the Green's function
wcal.Nk     = 50                      % number of k points along the k path 
wcal.kp1    = [0.0 0.0]               % the starting k point for the k path
wcal.kp2    = [0.5 0.0]               % the ending k point for the k path

wcal.Ef     = 5.2503                  % the Fermi energy 

wcal.Emu    = [-1 1]                  % the energy window
wcal.nE     = 50                      % number of energy points 

wcal.ncpu   = 2                       % number of cups for the paralle running
wcal.file   = Sftn58sparse_16L.mat    % the slab model gotten from "STBHmftn"

wcal.isPot  = 1                       % isPot = 1 adding an extra layer on top
wcal.PeV    = 0.1                     % the onsite energy added on the extra layer 
wcal.dis    = 0.86                    % to specify how many atom layers to be modified 

endGreenSur
```
2) Run "GreenSur.m" to get the Ek dispersion of the semi-infinite thin film with (hkl) surface
```
>> GreenSur.m
```
3) Run "PlotEK.m" to plot the spectrum weight
```
>> PlotEk.m
```
