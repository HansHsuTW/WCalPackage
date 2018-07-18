%%% This code will transfer ftn58 to a standard format (ftn58sparse) %%%
%%% ---------------------------------------------------------------- %%%
%%% Necessary Input files:                                           %%%
%%% 1) Wannier interpolation Hamiltonian: case_hr.dat                %%%
%%% 2) Crystal information: st_input                                 %%%
%%% ---------------------------------------------------------------- %%%
%%% Input data description:                                          %%%
%%% isEc ==> cutting of hopping strength                             %%%
%%%          (the magnitude is assigned by "ecut")                   %%%
%%% isSW ==> orbital  switcher                                       %%%
%%%          (this will re-arrange the orbital id by following the   %%%
%%%           standard form, see the table below)                    %%%
%%% isSO ==> option of various spin-orbital types                    %%%
%%%          ("0" is for non-SOC; "1" is for DFT-SOC;                %%%
%%%           "2" is for SOC added manually)                         %%%
%%% ---------------------------------------------------------------- %%%

%%% ---------------------------------------------------------------------- %%%
%%%                     Order of atomic orbitals: 
%%%    id   Wien2k       QE(w90)      OpenMX       TB2KKR
%%% s   1   s            s            s            s
%%% p   2   px           pz           px           px
%%%     3   py           px           py           py
%%%     4   pz           py           pz           pz
%%% d   5   dxy          dz2          dz2          dz2
%%%     6   dyz          dxz          d(x2-y2)     dxz
%%%     7   dxz          dyz          dxy          dyz
%%%     8   d(x2-y2)     d(x2-y2)     dxz          dxy
%%%     9   dz2          dxy          dyz          d(x2-y2)
%%% f  10   fy(3x2-y2)   fz3          fz3          fz3        (z(5z2-3r2))
%%%    11   fz(x2-y2)    fxz2         fxz2         fxz2       (x(5z2-3r2))  
%%%    12   fyz2         fyz2         fyz2         fyz2       (y(5z2-3r2))
%%%    13   fxz2         fz(x2-y2)    fz(x2-y2)    fz(x2-y2)  (z(x2-y2))
%%%    14   fxyz         fxyz         fxyz         fxyz       (2xyz)
%%%    15   fx(x2-3y2)   fx(x2-3y2)   fx(x2-3y2)   fx(x2-3y2) (x(x2-3y2))
%%%    16   fz3          fy(3x2-y2)   fy(3x2-y2)   fy(3x2-y2) (y(3x2-y2))
%%% --------------------------------------------------------------------- %%%
% clear all
function TBHmftn()

%% --- Input Arguments --- %%
% isEc     = 0;
% ecut     = 0.0005;
% isSO     = 2;
isSymSO  = 0;
isSymOn  = 0;
isSymCut = 0;
Rcut     = 15;
isSW     = 1;

wcal = ReadInput('input.txt');
isEc = wcal.isEC;
ecut = wcal.ecut;
isSO = wcal.isSO;
%%%--------------------------------------------------%%%

%% --- Choose your ftn58 resource --- %%
if 1==1
    [ftn58,~]=readwanhr(wcal.ref,0);   
%     save ftn58.mat
else 
    load ftn58.mat
end
%%%--------------------------------------------------%%%

%% --- Read unit vectors and atomic positions (from st_input.txt)
fid    = fopen('st_input.txt','r');
buf    = fgetl(fid);
system = buf;
buf    = fgetl(fid);
ABC    = str2num(buf);

for ii=1:3
    buf          = fgetl(fid);
    vector(ii,:) = str2num(buf);
end

buf = fgetl(fid);
atn = str2num(buf); % # of atom type

orbit = zeros(atn,12);
for ii=1:atn
    buf   = fgetl(fid);
    space = find(buf(:)==' ');
    norb  = length(space);
    atom{ii} = buf(1:(space(1)-1));
    orbid = str2num(buf((space(1)+1):end));
    orbit(ii,1:length(orbid)) = orbid;
end

buf = fgetl(fid);
nat = str2num(buf);

id2txt = {'s' 'px' 'py' 'pz' 'dz' 'dxz' 'dyz' 'dxy' 'dx2y2' ...
          'fz3' 'fxz2' 'fyz2' 'fz(x2-y2)' 'fxyz' 'fx(x2-3y2)' 'fy(3x2-y2)'};

id_at   = 1;
id_orb  = 1;
id_orbp = 1;
for ii=1:length(nat) %Type of atoms
    for jj=1:nat(ii) %# of atoms
        buf = fgetl(fid);
        ps  = str2num(buf);
        atominfo(id_at).Atom(1,:)       = atom{ii};
        atominfo(id_at).Position(1,1:3) = ps(1:3);
        atominfo(id_at).Norb            = ps(4);
        atominfo(id_at).OrbitIndex      = id_orb:(id_orb+ps(4)-1);
        
        temob    = [];
        orbit_n0 = find(orbit(ii,:)>0);
        for kk=1:length(orbit_n0)
            temob = strcat(temob,{' '},id2txt(1,orbit(ii,kk)));
        end
        atominfo(id_at).Orbit(1,:) = cellstr(temob);
        
        %%% Orbitps format:
        %%% #_of_orbital orbital_index(ftn58) type_of_atom frac_coord standard_orbital_id  
        for kk=1:ps(4)
            orbitps(id_orbp,1:7) = [id_orbp,id_orb+(kk-1),id_at,ps(1:3),orbit(ii,kk)];
            id_orbp = id_orbp + 1;
        end
        
        atominfo(id_at).OrbitID = orbit(ii,1:length(orbit_n0));
        id_at  = id_at + 1;
        id_orb = id_orb + ps(4);
    end 
end

%% --- Transfer orbital id (orbital switcher) --- %%
if isSW==1
    nbond   = ftn58(1,2);
    norb    = ftn58(1,1);
    %--- Rearrange the OrbitID to the standard order
    for iatom=1:length(atominfo)
        oid = sort(getfield(atominfo, {iatom}, 'OrbitID'));
        
        temob    = [];
        for kk=1:length(oid)
            temob = strcat(temob,{' '},id2txt(1,oid(kk)));
        end
        atominfo(iatom).Orbit(1,:) = cellstr(temob);
        atominfo(iatom).OrbitID = oid;       
    end
    
    
    if isSO==1
        orbitps = [orbitps;orbitps];%Update orbitals with SOC
        orbitps(id_orbp:end,1) = orbitps(id_orbp:end,1) + id_orbp-1;
        orbits  = sortrows(orbitps(1:norb/2,:),[3 7]);
        orbits  = [orbits;orbits];
        orbits(norb/2+1:end,1) = orbits(1:norb/2,1) + norb/2;
        TB2KKR  = [orbits(:,1) (1:norb)'];
        orbitps = [orbitps(:,1:6) orbits(:,7)]; 
    else
        orbits  = sortrows(orbitps,[3 7]);
        TB2KKR  = [orbits(:,2) (1:norb)'];
        orbitps = [orbitps(:,1:6) orbits(:,7)]; 
    end
    
    new_ij  = zeros(nbond,2);
    for ii=1:norb
        orb1_id = find(ftn58(2:end,2)==TB2KKR(ii,1)); 
        orb2_id = find(ftn58(2:end,3)==TB2KKR(ii,1));
        new_ij(orb1_id,1) = TB2KKR(ii,2);
        new_ij(orb2_id,2) = TB2KKR(ii,2);
    end

    ftn58(2:end,2:3) = new_ij(:,:);    
end
    
ftn58sparse.System = system;
ftn58sparse.Ainfo  = atominfo;
ftn58sparse.abc    = ABC;
ftn58sparse.BR     = vector;
ftn58sparse.Nat    = length(atominfo);
ftn58sparse.ver    = 'type1';
 
if isEc==1
nbond    = ftn58(1,2);
ibdelete = find(abs(ftn58(2:end,4))<ecut);
ftn58(ibdelete+1,:)=[];
nbond          = nbond-length(ibdelete);
ftn58(2:end,1) = (1:nbond)';
ftn58(1,2)     = nbond;
end

%% -- Actual Transfer Procedure (ftn58 to ftn58.sparse) --- %%

% Transfer to new format of ftn58.sparse (without SOC version)
if isSO==0    
    if isSymCut==1
        load ftn58sparse_raw.mat
        [distances,ib2dis,ftn58] = symcut(ftn58sparse);
        keepid        = find(ib2dis(:,2)<=Rcut);
        newftn58      = [ftn58(1,:);ftn58(keepid+1,:)];
        newftn58(1,2) = length(newftn58)-1;
        ftn58         = newftn58;
        figure;
        hist(ib2dis(:,2),length(distances));
    end
    
    if isSymOn==1
        ftn_0 = find(ftn58(2:end,5)==0&ftn58(2:end,6)==0&ftn58(2:end,7)==0);
        for ii=1:length(ftn_0)
            if (ftn58(ftn_0(ii)+1,2)~=ftn58(ftn_0(ii)+1,3)&& ...
                orbitps(ftn58(ftn_0(ii)+1,2),3)==orbitps(ftn58(ftn_0(ii)+1,3),3))
            ftn58(ftn_0(ii)+1,4) = 0.0;
            end
        end
%        ftn_d = find(abs(ftn58(2:end,6))>=6);
%        tmp_d = setdiff([1:length(ftn58)-1],ftn_d);
%        newftn58 = ftn58(tmp_d+1,:);
%        ftn58 = newftn58;        
    end
ftn58sparse.isSO=0;
ib1=find(ftn58(2:end,2) == ftn58(2:end,3))+1;
ib2=find(ftn58(2:end,2) ~= ftn58(2:end,3))+1;
ftn58sparse.norb    = ftn58(1,1);
ftn58sparse.Orbitps = orbitps;
ftn58sparse.ij      = [ftn58(ib1,2:3);ftn58(ib2,2:3);ftn58(ib2,[3 2])];
ftn58sparse.tt      = real([ftn58(ib1,4);ftn58(ib2,4);ftn58(ib2,4)]);
ftn58sparse.dd      = [ftn58(ib1,5:7);ftn58(ib2,5:7);-ftn58(ib2,5:7)];
end

%Transfer to new format of ftn58.sparse (DFT_SOC version)
if isSO==1
    if isSW~=1
        orbitps = [orbitps;orbitps];%Update orbitals with SOC
        orbitps(id_orbp:end,1) = orbitps(id_orbp:end,1) + id_orbp-1;
    end
    
    if isSymCut==1
        load ftn58sparse.mat
        [distances,ib2dis,ftn58] = symcut(ftn58sparse);
        keepid        = find(ib2dis(:,2)<=Rcut);
        newftn58      = [ftn58(1,:);ftn58(keepid+1,:)];
        newftn58(1,2) = length(newftn58)-1;
        ftn58         = newftn58;
        figure;
        hist(ib2dis(:,2),length(distances));
    end
%orbitps(id_orb:end,2) = orbitps(id_orb:end,2) + id_orb-1;

    if isSymOn==1
        ftn_0 = find(ftn58(2:end,5)==0&ftn58(2:end,6)==0&ftn58(2:end,7)==0);
        for ii=1:length(ftn_0)
            if (ftn58(ftn_0(ii)+1,2)~=ftn58(ftn_0(ii)+1,3)&& ...
                orbitps(ftn58(ftn_0(ii)+1,2),3)==orbitps(ftn58(ftn_0(ii)+1,3),3))
            ftn58(ftn_0(ii)+1,4) = 0.0;
            end
        end
    end
    
    if isSymSO==1
        ftn_SO = find(ftn58(2:end,5)~=0|ftn58(2:end,6)~=0|ftn58(2:end,7)~=0);
        ftn58(ftn_SO+1,4) = real(ftn58(ftn_SO+1,4));
    end

ftn58sparse.isSO=1;
ib1=find(ftn58(2:end,2) == ftn58(2:end,3))+1;
ib2=find(ftn58(2:end,2) ~= ftn58(2:end,3))+1;
ftn58sparse.norb        =  ftn58(1,1);

ftn58sparse.Orbitps = orbitps;
ftn58sparse.ij      = [ftn58(ib1,2:3);ftn58(ib2,2:3);ftn58(ib2,[3 2])];
ftn58sparse.tt      = [ftn58(ib1,4);ftn58(ib2,4);conj(ftn58(ib2,4))];
ftn58sparse.dd      = [ftn58(ib1,5:7);ftn58(ib2,5:7);-ftn58(ib2,5:7)];
end

% Transfer to new format of ftn58.sparse (manually SOC version)
% Add SOC %
if isSO==2
lamdaPbTe = [1.34 0.7361];    
lamdaSnTe = [0.4541 0.7361];
ipp       = [1,4];
[ftn58, ftn58SO]=SOftn58(1,ftn58,lamdaSnTe,ipp);

ftn58sparse.isSO=1;
ib1=find(ftn58(2:end,2)==ftn58(2:end,3))+1; % +1 for correct index 
ib2=find(ftn58(2:end,2)~=ftn58(2:end,3))+1;
ftn58sparse.norb=ftn58(1,1);

orbitps = [orbitps;orbitps];%Update orbitals with SOC
orbitps(id_orbp:end,1) = orbitps(id_orbp:end,1) + id_orbp-1;

ftn58sparse.Orbitps = orbitps;
ftn58sparse.ij=[ftn58SO(2:end,2:3);ftn58(ib1,2:3);ftn58(ib2,2:3);ftn58SO(2:end,[3 2]);ftn58(ib2,[3 2])];
ftn58sparse.tt=[ftn58SO(2:end,4);ftn58(ib1,4);ftn58(ib2,4);conj(ftn58SO(2:end,4));ftn58(ib2,4)];
ftn58sparse.dd=[ftn58SO(2:end,5:7);ftn58(ib1,5:7);ftn58(ib2,5:7);-ftn58SO(2:end,5:7);-ftn58(ib2,5:7)];
end

save ftn58sparse.mat ftn58sparse