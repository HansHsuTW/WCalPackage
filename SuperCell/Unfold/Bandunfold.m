%%%           Band Unfolding Program           %%%
%%% ------------------------------------------ %%%
%%% Reference: PhysRevLett.104.216401          %%%
%%% ------------------------------------------ %%%
%%%  1 Nov created by Wing Chi                 %%%
%%% 15 Nov modified by Hans                    %%%
%%% ------------------------------------------ %%%
clear all

%% Inputs %%
equi_atom = [1 2 3 4;
             5 6 7 8;
             9 10 11 12];
InFile    = 'ftn58sparse_2x2_pristine.mat';   
latss     = [2 0 0;0 2 0;0 0 1];

%% k-paths %%
nk = 50;
p1 = [linspace(1/3,0,nk+1)' linspace(1/3,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p2 = [linspace(0.0,0.5,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p3 = [linspace(0.5,1/3,nk)' linspace(0.0,1/3,nk)' linspace(0.0,0.0,nk)'];

kpt_uc = [p1(1:nk,:);p2(1:nk,:);p3(1:end,:)];

%% Loading
instruct = load(InFile);
if isfield(instruct,'ftn58sparse')
    ftn58sparse = instruct.ftn58sparse;
    ftn58sparse.latss = latss;
elseif isfield(instruct,'SPftn58sparse')
    ftn58sparse = instruct.SPftn58sparse;
else
    ftn58sparse = instruct.Sftn58sparse;
end

Orbitps   = ftn58sparse.Orbitps;
norbss    = ftn58sparse.norb;
ijss      = ftn58sparse.ij;
ddss      = ftn58sparse.dd;
ttss      = ftn58sparse.tt;
BR        = ftn58sparse.BR;
latss     = ftn58sparse.latss;

Orbitpsss = Orbitps;

%--- 1/2 cell ---%
%--- up spin ---%
for iat=1:size(equi_atom,2)-1
    for ispe=1:size(equi_atom,1)
        iat1 = equi_atom(ispe,2:end);%9
        at1 = find(Orbitps(1:norbss/2,3)==iat1(iat));
        at2 = find(Orbitps(1:norbss/2,3)==equi_atom(ispe,1));
        Orbitpsss(at1,3) = Orbitps(at2,3);  
    end
end
%--- down spin ---%
for iat=1:size(equi_atom,2)-1
    for ispe=1:size(equi_atom,1)
        iat1 = equi_atom(ispe,2:end);%9
        at1 = find(Orbitps(1+norbss/2:end,3)==iat1(iat));
        at2 = find(Orbitps(1+norbss/2:end,3)==equi_atom(ispe,1));
        Orbitpsss(at1,3) = Orbitps(at2+norbss/2,3);  
    end
end

%% generate supercell klist corresponding to the unit cell k
nks    = size(kpt_uc,1);
% --- calculate the kpoints w.r.t. super cell reciprocal vector
kpt_sc = (latss*(kpt_uc(:,1:3))')';

% --- find the supercell orbit index which arise from the same unit cell orbital
rep_orb = unique(Orbitpsss(:,2));
ncell   = norbss/length(rep_orb);
norb_uc = length(rep_orb);
for ii=1:norb_uc
    map_orb(ii,:) = find(Orbitpsss(:,2)==rep_orb(ii));
end

%% calculate eigenvalues and eigenfunctions of supercell H in the sc klist %%
eigvec = cell(nks,1);
Ek     = zeros(nks,norbss);
partial_weight=zeros(nks,norbss,norb_uc);
Weight = zeros(nks,norbss);

tic
for ik=1:nks
    time_init=tic;
    kcolumnvec = kpt_sc(ik,1:3)'*2*pi;
    Hsparse    = sparse(ijss(:,1),ijss(:,2),exp(1i*ddss*kcolumnvec).*ttss,norbss,norbss);
    HH         = full(Hsparse);
    HH         = (HH+HH')/2;
    [vec, Etmp]= eig(HH);
    Ek(ik,:)   = diag(Etmp);
    eigvec{ik,1} = vec;
     
    kcolvec = kpt_uc(ik,:)'*2*pi;
    for iband=1:norbss
        Evec = eigvec{ik,1}(:,iband);
        rvec = Orbitpsss(:,4:6)*latss;
        tmp  = exp(-1*1i*rvec*kcolvec).*Evec;
        for ii=1:norb_uc
            partial_weight(ik,iband,ii) = abs(sum(tmp(map_orb(ii,:))))^2;
        end
        Weight(ik,iband) = sum(partial_weight(ik,iband,1:end));
%         Weight(ik,iband,:) = partial_weight(ik,iband,1:end);
    end
     
     if mod(ik,10)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
     end
end

clear vec Etmp
toc

save unfold.mat Weight norbss nks Ek