%%% Calculate Maximum hopping distance %%%
load ftn58sparse.mat
BRori  = ftn58sparse.BR;

load slab_info.mat
BRslab = BR;

ddori  = ftn58sparse.dd;
ddori  = BRori'*ddori';
ddslab = inv(BR')*ddori;

max(ddslab')
