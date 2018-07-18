%%% Program for Calculating Z2 via Wilson Loop %%%
%%% ---------------------------------------- %%%
clear all;
load ftn58sparse_SnTe.mat

%%-- Inital Setup --%%
nkx    = 100;
nky    = 150;
ocnorb = 6;

%%-- Actualy run --%%
WL(ftn58sparse,ocnorb,nkx,nky);