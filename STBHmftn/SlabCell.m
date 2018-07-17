%%%            Program for Constructing Slab Super-unitcell          %%%
%%% ---------------------------------------------------------------- %%%
%%% Necessary Input files:                                           %%%
%%% ftn58sparse.mat                                                  %%%
%%% ---------------------------------------------------------------- %%%
%%% Input data description:                                          %%%
%%% nlayer    ==> number of layers in the slab supercell             %%%
%%% isCutSlab ==> if true, cutting a chosen surface                  %%%
%%%          ("uplimit" and "botlimit" are set for indicating your   %%%
%%%           up and bottom coordinates along hkl direction;         %%%
%%%           "xyz_dir" indicates the xyz cartesian axis)            %%%
%%% ---------------------------------------------------------------- %%%
%%% 09/23/16 Updated:                                                %%%
%%% Slab cutting is based on the frac_corrd of 3 unit vectors.       %%%
%%% Orbit positions are sort according to the vector_c.              %%%
%%% (xyz_dir = 3)                                                    %%%
%%% ---------------------------------------------------------------- %%%
function SlabCell()

%% --- Input Arguments --- %%
wcal = ReadInput('input.txt');
nlayer    = wcal.nL;
isCutSlab = wcal.isCut;
isVASP    = 1;
Rinv      = 1;

% nlayer    = 10;
% isCutSlab = 0;
% uplimit   = 8.1;
% botlimit  = 0.1;

LLB = wcal.LB(2:end-1);
sid = find(isspace(LLB));
botlimit = str2double(LLB(1:sid(1)-1)); 
uplimit  = str2double(LLB(sid(1)+1:end));

MI  = wcal.hkl(2:end-1);
sid = find(isspace(MI));
hkl = [str2double(MI(1:sid(1)-1)) str2double(MI(sid(1)+1:sid(2)-1)) ...
       str2double(MI(sid(2)+1:end))];


%%% Definition of slab unit vector (Based on conventional unit cell):
%%% 1) The first(S1) and second(S2) unit vectors of slab are orthogonal and
%%%    perpendicular to [hkl] or parallel to (hkl) surface.
%%% 2) Part of the slab's third(S3) unit vector parallels to [hkl].
plane_h = hkl(1);
plane_k = hkl(2);
plane_l = hkl(3);
%%%--------------------------------------------------%%%

%% --- Actual Procedure --- %%%
load(wcal.ref);

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

vector = ftn58sparse.BR;
so     = ftn58sparse.isSO;
abc    = ftn58sparse.abc;
xyz_dir   = 3;

% RT     = [1 -1 0;1 1 0;0 0 sqrt(2)]/sqrt(2);
% vector = vector*RT';
%for i=1:3
%    vector(:,i) = vector(:,i)*abc(1,i);
%end

vectora   = vector(1,:);
vectorb   = vector(2,:);
vectorc   = vector(3,:);
bulkbasis = [vectora; vectorb; vectorc];
xyz2abc   = inv(bulkbasis); %[xyz2abc][a;b;c] = I

%%% For the case that conventional unit cell is the same as its unit cell
if wcal.isConv
   xyz2abc = [1 0 0;0 1 0;0 0 1]; %[xyz2abc][a;b;c] = I 
end
%%% ----------------------------------------------------------------- %%%

abc2xyz   = inv(xyz2abc);%abc2xyz = bulkbasis
xyzvector = xyz2abc * bulkbasis; % to conventional cell (xyz2abc*abc2xyz) 

direction_index  = [plane_h, plane_k, plane_l];
direction_vector = plane_h * xyzvector(1,:) + plane_k * xyzvector(2,:) + plane_l * xyzvector(3,:);

% Construct a vector grid to search new unit cell (orthorhombic !?) 
nbound    = 2;
neighbor1 = [];
neighbor2 = [];
for ll = -nbound:nbound
    for mm = -nbound:nbound
        for nn = -nbound:nbound
            
            basis_index  = ll * abc2xyz(1,:) + mm * abc2xyz(2,:) + nn * abc2xyz(3,:); % ll*a+mm*b+nn*c 
            basis_vector = ll * vectora + mm * vectorb + nn * vectorc; % a,b,c in cart
            
            if((abs(dot(basis_index, direction_index))<=0.001)&&((abs(ll)>0)||(abs(mm)>0)||(abs(nn)>0)))
                neighbor1 = [neighbor1; ll, mm, nn, basis_vector, norm(basis_vector)];
            end
            if(abs(dot(basis_index, direction_index))>=0.0001)
                neighbor2 = [neighbor2; ll, mm, nn, basis_vector, norm(basis_vector)];
            end
            
        end
    end
end

[~, ii] = min(neighbor1(1:end, 7));
basis1  = neighbor1(ii, 4:6); %%% first basis vector

for ii=1:size(neighbor1, 1)
    theta = acos(basis1 * neighbor1(ii, 4:6)'/norm(basis1)/norm(neighbor1(ii, 4:6)));
    neighbor1(ii, 8) = theta;
end

candidate_basis2 = [];
for ii=1:size(neighbor1, 1)
    if((neighbor1(ii,8)>0.001)&&(neighbor1(ii,8)<=(pi/2+0.01)))
        candidate_basis2 = [candidate_basis2; neighbor1(ii,4:8)];
    end
end

[~, ii] = min(candidate_basis2(1:end, 4));
basis2  = candidate_basis2(ii, 1:3); %%% second basis vector

candidate_basis3 = [];
for ii=1:size(neighbor2, 1)
    theta1 = acos(basis1 * neighbor2(ii, 4:6)'/norm(basis1)/norm(neighbor2(ii, 4:6)));
    theta2 = acos(basis2 * neighbor2(ii, 4:6)'/norm(basis2)/norm(neighbor2(ii, 4:6)));
    neighbor2(ii, 8) = theta1;
    neighbor2(ii, 9) = theta2;
    
    if((theta1>0.001)&&(theta1<=pi/2)&&(theta2>0.001)&&(theta2<=pi/2))
        candidate_basis3 = [candidate_basis3; neighbor2(ii, 4:9)];
    end
end 

[temp,ii] = min(candidate_basis3(1:end, 4));
if (Rinv)
    basis3    = -candidate_basis3(ii,1:3); % final basis vector
else
    basis3    = candidate_basis3(ii,1:3); % final basis vector
end
BR  = [basis1;basis2;basis3];
hkl = [plane_h plane_k plane_l];

%%% Manually Adjust %%%
%basis3 = direction_vector;
%%% Manually Adjust %%%

transform_matrix = [basis1; basis2; basis3]';

%%%--- Change basis(from original to slab)
%%%--- Format (run over all atoms in unit cell) ---%%%
%%%--- 1:3 => original coord. 4 => # of orbitals 5:7 => coord. of new basis  
for ii=1:size(ps, 1)
    ps(ii, 5:7) = (ps(ii,1)*vectora + ps(ii,2)*vectorb + ps(ii,3)*vectorc) * transpose(inv(transform_matrix));%%% express the position in terms of basis vector
end

%%%--- new vector grid
%%%--- ps1 is for a table for atoms in a 3D grid 
ps1 = [];
nbound = 4;
for ii=-nbound:nbound
    for jj=-nbound:nbound
        for kk=-nbound:nbound            
            for ss=1:size(ps,1)
                ps1 = [ps1; ps(ss, 1)+ii, ps(ss, 2)+jj, ps(ss, 3)+kk, ps(ss, 4), ss];
            end            
        end
    end
end

for ii=1:size(ps1, 1)
    ps1(ii,6:8) = ps1(ii,1)*vectora + ps1(ii,2)*vectorb + ps1(ii,3)*vectorc;
    ps1(ii,9:11) = ps1(ii,6:8) * transpose(inv(transform_matrix));
end

% Create atoms in the slab unit cell
atomps = []; 
%%% cloumn 1-3 in terms of lattice vector of bulk unit cell, 
%%% cloumn 4 is the # of orbitals for the atom,
%%% cloumn 5 labels the atom,
%%% column 6-8 is the position in xyz coordinate (BR basis)
%%% cloumn 9-11 in terms of basis vector of slab unit cell
              
error = -1e-3;             
for ii=1:size(ps1, 1)
    if(((ps1(ii,9)>error)&&(ps1(ii,9)<=1.00+error))&&((ps1(ii,10)>error)&&(ps1(ii,10)<=1.00+error))&&((ps1(ii,11)>error)&&(ps1(ii,11)<=(1.0+error))))
        atomps = [atomps; ps1(ii,1:end)]; % result in a scheme of atom position within (1,1,1) 
    end
end

Aplot('Unit Cell',atomps(:,6:8)*[abc(1) 0 0;0 abc(2) 0;0 0 abc(3)]');     

%%%--- Cut slab structure
%%%--- superatoms has the same format as the atomps
superatomps = [];
for ii=0:nlayer-1
    for jj=1:size(atomps, 1)
        superatomps = [superatomps; atomps(jj, 1:10), atomps(jj, 11)+ii];
    end
end

superatomps   = sortrows(superatomps, 5);
super_cart    = [];
super_cart_rp = [];
for ii=1:size(superatomps)
    super_cart = [super_cart;superatomps(ii,9)*basis1 + superatomps(ii,10)*basis2 + ...
                     superatomps(ii,11)*basis3];
    for jj=-1:1
        super_cart_rp = [super_cart_rp;(superatomps(ii,9)+jj)*basis1 + superatomps(ii,10)*basis2 + ...
                         superatomps(ii,11)*basis3];
    end
end
super_cart(:,:)    = super_cart(:,:)*[abc(1) 0 0;0 abc(2) 0;0 0 abc(3)];
super_cart_rp(:,:) = super_cart_rp(:,:)*[abc(1) 0 0;0 abc(2) 0;0 0 abc(3)];

%%% For cutting specific surfaces on both sides %%%
if isCutSlab==1
superatomps   = superatomps(superatomps(1:end,xyz_dir+8)>botlimit & superatomps(1:end,xyz_dir+8)<uplimit, 1:end);
super_cart_rp =[];
for ii=1:size(superatomps,1)
    for jj=-1:1
        super_cart_rp = [super_cart_rp;(superatomps(ii,9)+jj)*basis1 + superatomps(ii,10)*basis2 + ...
                         superatomps(ii,11)*basis3];
    end
end
super_cart_rp(:,:) = super_cart_rp(:,:)*[abc(1) 0 0;0 abc(2) 0;0 0 abc(3)];
end
%%% ---------------------------------------------- %%%

% Aplot('SuperCell',super_cart_rp);

% Check structure (extending with 4 x 4 x 4)
if 1==1
Catomps = [];
for ii=0:3
    for jj=0:3
        for mm=0:3
            Catomps = [Catomps; atomps(:,1:8), atomps(:,9)+ii, atomps(:,10)+jj, atomps(:,11)+mm];
        end
    end
end

for ii=1:size(Catomps)
    CA_cart(ii,:) = Catomps(ii,9)*basis1 + Catomps(ii,10)*basis2 + Catomps(ii,11)*basis3;
end

%Aplot('Check',CA_cart*[abc(1) 0 0;0 abc(2) 0;0 0 abc(3)]');
end

%% --- Assign postion for each orbital --- %%%
%%%--- Format of orbitps 
%%%--- 1     => orbital ID 
%%%--- 2,3   => original orbital ID and atom ID in ftn58sparse
%%%--- 4:6   => orbital frac. coord. 
%%%--- 7:9   => atom frac. coord. 
%%%--- 10    => orbital standard ID (see table)
%%%--- 11:13 => xyz corrd. (BR base) 
orbitps     = [];
superatomps = sortrows(superatomps, xyz_dir+5);
for ii=1:size(superatomps, 1)   %total # of super-atoms
    for jj=1:superatomps(ii, 4) %total # of orbitals belong to each super-atom
        orbitps = [orbitps; sum(superatomps(1:ii-1,4)) + jj, sum(ps(1:superatomps(ii,5)-1,4))+jj,...
                   ii, superatomps(ii,9:11), ps(superatomps(ii,5),5:7),orbID(superatomps(ii,5),jj),...
                   superatomps(ii,9:11)*transform_matrix'];
    end
end
orbitps      = sortrows(orbitps, 6);
orbitps(:,1) = [1:size(orbitps,1)];

if so==1
norbit_total = size(orbitps, 1); %% the number of orbitals in slab unit cell(spin counterpart not included)
orbitps = [orbitps; orbitps];    %% include the spin degree
orbitps(norbit_total+1:end, 1) = orbitps(1:norbit_total, 1) + norbit_total;
orbitps(norbit_total+1:end, 2) = orbitps(1:norbit_total, 2) + sum(ps(1:end, 4));
end

%% --- Generating the vasp POSCAR of unit and supercell --- %%%
if (isVASP)
    ftn2vasp(atomps,superatomps,ftn58sparse,BR,nlayer)
end

save slab_info.mat superatomps orbitps transform_matrix atomps bulkbasis BR hkl
