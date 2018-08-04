%%%            Program for Constructing Super-unitcell               %%%
%%% ---------------------------------------------------------------- %%%
%%% Necessary Input files:                                           %%%
%%% ftn58sparse.mat                                                  %%%
%%% ---------------------------------------------------------------- %%%
%%% Input data description:                                          %%%
%%% U2STR  ==> The transfrom matrix from unitcell to supercell       %%%
%%% ---------------------------------------------------------------- %%%
%%% 19 Nov created by Hans                                           %%%
%%% ---------------------------------------------------------------- %%%
function SuperCell()

%% --- Input Argments --- %%
wcal = ReadInput('input.txt');
v1   = wcal.vec1(2:end-1);
sid  = find(isspace(v1));
vec1 = [str2double(v1(1:sid(1)-1)) str2double(v1(sid(1)+1:sid(2)-1)) ...
        str2double(v1(sid(2)+1:end))];
v2   = wcal.vec2(2:end-1);
sid  = find(isspace(v2));
vec2 = [str2double(v2(1:sid(1)-1)) str2double(v2(sid(1)+1:sid(2)-1)) ...
        str2double(v2(sid(2)+1:end))];
v3   = wcal.vec3(2:end-1);
sid  = find(isspace(v3));
vec3 = [str2double(v3(1:sid(1)-1)) str2double(v3(sid(1)+1:sid(2)-1)) ...
        str2double(v3(sid(2)+1:end))];

U2STR   = [vec1;vec2;vec3];
Dis     = [0 0 0.0;];
isVASP  = 1;

%% --- Actual Procedure --- %%%
instruct    = load(wcal.file);
if isfield(instruct,'ftn58sparse')
    ftn58sparse = instruct.ftn58sparse;
elseif isfield(instruct,'SPftn58sparse')
    ftn58sparse = instruct.SPftn58sparse;
else
    ftn58sparse = instruct.Sftn58sparse;
end

orbID = zeros(ftn58sparse.Nat,5);
for i=1:ftn58sparse.Nat
    temp  = getfield(ftn58sparse.Ainfo, {i}, 'Position');
    temno = getfield(ftn58sparse.Ainfo, {i}, 'Norb');
    temid = getfield(ftn58sparse.Ainfo, {i}, 'OrbitID');    
    ps(i,1:3) = temp;
    ps(i,4)   = temno;
    orbit_n0  = find(temid(:)>0);
    orbID(i,1:length(orbit_n0)) = temid;
end

BR  = ftn58sparse.BR;
so  = ftn58sparse.isSO;
abc = ftn58sparse.abc;
orb = ftn58sparse.Orbitps;
orb(:,2) = [1:size(orb,1)]';

%%%--- Change basis(from original to slab)
%%%--- Format (run over all atoms in unit cell) ---%%%
%%%--- 1:3 => original coord. 4 => # of orbitals 5:7 => coord. of new basis  
for ivec=1:3
    NBR(ivec,:) = BR(1,:)*U2STR(ivec,1) + BR(2,:)*U2STR(ivec,2) + BR(3,:)*U2STR(ivec,3);  
end

for ii=1:size(ps, 1)
    ps(ii, 5:7) = ps(ii,1:3)*BR/(NBR);
end
ps(:,8) = [1:size(ps,1)]';

%%%--- new vector grid
%%%--- ps1 is for a table for atoms in a 3D grid 
ps1 = [];
nbound = 4;
for ii=-nbound:nbound
    for jj=-nbound:nbound
        for kk=-nbound:nbound            
            for ss=1:size(ps,1)
                ps1 = [ps1; ps(ss, 1)+ii, ps(ss, 2)+jj, ps(ss, 3)+kk, ps(ss, 4), ps(ss,8)];
            end            
        end
    end
end

for ii=1:size(ps1, 1)
    ps1(ii,6:8)  = ps1(ii,1:3)*BR;
    ps1(ii,9:11) = ps1(ii,6:8)/(NBR);
end

% Create atoms in the slab unit cell
atomps = []; spatomps=[];
%%% column 1-3 in terms of lattice vector of bulk unit cell, 
%%% column 4 is the # of orbitals for the atom,
%%% column 5 labels the atom,
%%% column 6-8 is original atom position (BR basis)
%%% column 9-11 in terms of basis vector of supercell
%%% column 12-14 supercell position 
              
error = -1e-3;
ps1(:,9:11) = ps1(:,9:11) - repmat(Dis,size(ps1,1),1);
for ii=1:size(ps1, 1)
    if(((ps1(ii,9)>error)&&(ps1(ii,9)<=1.00+error))&&((ps1(ii,10)>error)&&(ps1(ii,10)<=1.00+error))&&((ps1(ii,11)>error)&&(ps1(ii,11)<=(1.0+error))))
        atomps = [atomps; ps1(ii,1:end)]; % result in a scheme of atom position within (1,1,1) 
    end
end

spid = setdiff([1:size(ps,1)],atomps(:,5));
for ii=1:length(spid)
    spatomps = [spatomps;[ps(spid(ii),1:4) spid(ii) ps(spid(ii),5:7) ps(spid(ii),5:7) ps(spid(ii),5:7)*NBR]];
end

% natom = unique(atomps(:,5));
% for iat=1:length(natom)
%     atid = find(atomps(:,5)==natom(iat));
%     atomps(atid,6:8) = repmat(atomps(atid(1),6:8),length(atid),1);
% end
for iat=1:size(atomps,1)
    atomps(iat,6:8)   = ps(atomps(iat,5),5:7);
    atomps(iat,12:14) = atomps(iat,9:11)*NBR;
end

% Aplot('Unit Cell',atomps(:,12:14)*[abc(1) 0 0;0 abc(2) 0;0 0 abc(3)]');   

%% --- Assign postion for each orbital --- %%%
%%%--- Format of orbitps 
%%%--- 1     => orbital ID 
%%%--- 2,3   => original orbital ID and atom ID in ftn58sparse
%%%--- 4:6   => orbital frac. coord. 
%%%--- 7:9   => atom frac. coord. 
%%%--- 10    => orbital standard ID (see table)
%%%--- 11:13 => xyz corrd. (BR base) 
orbitps     = [];
atomps = sortrows(atomps, 5);
for ii=1:size(atomps, 1)   %total # of super-atoms
    initorb = find(orb(:,3)==atomps(ii,5));
    for jj=1:atomps(ii, 4) %total # of orbitals belong to each super-atom
        orbitps = [orbitps; sum(atomps(1:ii-1,4)) + jj, orb(initorb(1),2)+jj-1,...
                   ii, atomps(ii,9:11), ps(atomps(ii,5),5:7),orbID(atomps(ii,5),jj),...
                   atomps(ii,9:11)*NBR];
    end
end
% orbitps      = sortrows(orbitps, 6);
orbitps(:,1) = [1:size(orbitps,1)];

if so==1
norbit_total = size(orbitps, 1); %% the number of orbitals in slab unit cell(spin counterpart not included)
orbitps = [orbitps; orbitps];    %% include the spin degree
orbitps(norbit_total+1:end, 1) = orbitps(1:norbit_total, 1) + norbit_total;
orbitps(norbit_total+1:end, 2) = orbitps(1:norbit_total, 2) + sum(ps(1:end, 4));
end

%% --- Prepare Superatoms Info (for inverted process) --- %%
if ~isempty(spatomps)
    sptable = [];
    nat = size(spatomps,1);
    for ii=1:size(atomps,1)
        component = repmat(atomps(ii,9:11),nat,1) - spatomps(:,9:11);
        component = round(component,2);
        parter    = find(mod(component(:,1),1)==0&mod(component(:,2),1)==0&mod(component(:,3),1)==0);
        ib1     = find(orb(:,3)==atomps(ii,5));
        sptable = [sptable; ib1 ib1];
        for ip=1:length(parter)
            ib2 = find(orb(:,3)==spatomps(parter(ip),5));
            sptable = [sptable;ib2 ib1];
            spatomps(parter(ip),6:8) = atomps(ii,6:8);
        end
    end
    
    % --- Creating the talbe for recording the orbit relations
    % --- 1st column -> orbit ID from the input model
    % --- 2nd column -> the corresponding orbit ID in small cells
    for ii=1:size(sptable,1)
        id = find(orbitps(:,2)==sptable(ii,2));
        sptable(ii,2) = id;
    end
    
    sporbitps     = [];
    spatomps = sortrows(spatomps, 5);
    for ii=1:size(spatomps, 1)   %total # of super-atoms
        initorb = find(orb(:,3)==spatomps(ii,5));
        for jj=1:spatomps(ii, 4) %total # of orbitals belong to each super-atom
            sporbitps = [sporbitps; sum(spatomps(1:ii-1,4)) + jj, orb(initorb(1),2)+jj-1,...
                         spatomps(ii,5), spatomps(ii,6:8), ps(spatomps(ii,5),5:7),orbID(spatomps(ii,5),jj),...
                         spatomps(ii,9:11)*NBR];
        end
    end
    sporbitps      = sortrows(sporbitps, 6);
    sporbitps(:,1) = [1:size(sporbitps,1)];
    
    if so==1
        spnorbit_total = size(sporbitps, 1); %% the number of orbitals in slab unit cell(spin counterpart not included)
        sporbitps = [sporbitps; sporbitps];  %% include the spin degree
        sporbitps(spnorbit_total+1:end, 1) = sporbitps(1:spnorbit_total, 1) + spnorbit_total;
        sporbitps(spnorbit_total+1:end, 2) = sporbitps(1:spnorbit_total, 2) + sum(ps(1:end, 4));
    end
end

%% --- Generating the vasp POSCAR of unit and supercell --- %%%
if (isVASP)
    Super2vasp(atomps,ftn58sparse,NBR)
end

if isempty(spatomps)
   sptable   = [[1:size(orbitps,1)]' [1:size(orbitps,1)]'];
   sporbitps = [];
end
save slab_info.mat orbitps sporbitps sptable atomps BR NBR U2STR
