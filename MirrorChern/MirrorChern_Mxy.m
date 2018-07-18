%%% Functon for calculating the Berry curvature along x-y                %%%
%%% -------------------------------------------------------------------- %%%
%%% Define Berry curvature Omega(i,j), where i and j stands for          %%%
%%% the direction of two current directions. (i,j=x,y,z)                 %%%
%%% -------------------------------------------------------------------- %%%
function [Berry,Ek] = MirrorChern_Mxy(kpoints,ftn58sparse)

%% --- Initial setup for QSHC calcultaion --- %%%
BR   = ftn58sparse.BR;
abc  = ftn58sparse.abc; 
Norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
tt   = ftn58sparse.tt;
dd   = ftn58sparse.dd;
gam  = 1e-5;
RAng = pi/4;

%--- Here the k vector is no longer dimensionless 
T  = [BR(:,1)*abc(1) BR(:,2)*abc(2) BR(:,3)*abc(3)]; 
DD = dd*T;
Rot = [cos(-RAng) -sin(-RAng);sin(-RAng) cos(-RAng)];
XX  = DD(:,1);
YY  = DD(:,2);
DX  = Rot(1,1)*XX + Rot(1,2)*YY;
DY  = Rot(2,1)*XX + Rot(2,2)*YY;
DZ = DD(:,3);

%--- Mxy operator
Mxy = zeros(6);
sx  = [0 1;1 0];
sy  = [0 -1i;1i 0];

pM = [0 1 0;1 0 0;0 0 1];

Mxy(1:3,1:3) = pM;
Mxy(4:6,4:6) = pM;

MMxy = kron((sx-sy),Mxy);

[Uxy,~] = eig(MMxy);

%% --- Actual Procedure --- %%%
Berry = zeros(Norb/2,6);
Ek    = zeros(Norb); 

kp(1:3) = kpoints(:);
Hsparse = sparse(ii,jj,exp(1i*dd*kp').*tt,Norb,Norb);
H0      = full(Hsparse);
H0      = (H0 + H0')/2;
[~,D]   = eig(H0);
Ek      = diag(D);   
            
vx = 1i*full(sparse(ii,jj,DX.*exp(1i*dd*kp').*tt,Norb,Norb));    % dH/dkx
vy = 1i*full(sparse(ii,jj,DY.*exp(1i*dd*kp').*tt,Norb,Norb));    % dH/dky
vz = 1i*full(sparse(ii,jj,DZ.*exp(1i*dd*kp').*tt,Norb,Norb));    % dH/dkz

vx = Uxy\vx*Uxy;
vy = Uxy\vy*Uxy;
vz = Uxy\vz*Uxy;

vx = (vx + vx')/2;
vy = (vy + vy')/2;
vz = (vz + vz')/2;

%%% Separate +/- Eigenstates %%%
Hxy  = Uxy\H0*Uxy;
Hu   = Hxy(1:Norb/2,1:Norb/2); Hu = (Hu + Hu')/2;
Hd   = Hxy(Norb/2+1:end,Norb/2+1:end); Hd = (Hd + Hd')/2;
[eu,Eu] = eig(Hu); Eu = diag(Eu);
[ed,~]  = eig(Hd);
U_o   = [eu;zeros(Norb/2)];
U_e   = [zeros(Norb/2);ed];
%%%            end            %%%

V_x_e = U_e'*vx*U_e;  % matrix elements for velocity for E field
V_x_o = U_o'*vx*U_o;  % matrix elements for velocity for E field
V_y_e = U_e'*vy*U_e;  % matrix elements for velocity for E field
V_y_o = U_o'*vy*U_o;  % matrix elements for velocity for E field
V_z_e = U_e'*vz*U_e;  % matrix elements for velocity for E field
V_z_o = U_o'*vz*U_o;  % matrix elements for velocity for E field
                        
ch_xy_e = zeros(Norb/2,Norb/2); 
ch_xy_o = zeros(Norb/2,Norb/2); 
ch_yz_e = zeros(Norb/2,Norb/2); 
ch_yz_o = zeros(Norb/2,Norb/2); 
ch_zx_e = zeros(Norb/2,Norb/2); 
ch_zx_o = zeros(Norb/2,Norb/2); 

for nb1 = 1:Norb/2
    for nb2 = nb1+1:Norb/2     % half-part of the matrix
        diff_energy = Eu(nb1) - Eu(nb2); 
        
        common_Vxp_e = V_x_e(nb2,nb1)/(diff_energy^2+1i*gam);
        common_Vyp_e = V_y_e(nb2,nb1)/(diff_energy^2+1i*gam);
        common_Vzp_e = V_z_e(nb2,nb1)/(diff_energy^2+1i*gam);
        common_Vxn_e = V_x_e(nb2,nb1)/(diff_energy^2-1i*gam);
        common_Vyn_e = V_y_e(nb2,nb1)/(diff_energy^2-1i*gam);
        common_Vzn_e = V_z_e(nb2,nb1)/(diff_energy^2-1i*gam);
        common_Vxp_o = V_x_o(nb2,nb1)/(diff_energy^2+1i*gam);
        common_Vyp_o = V_y_o(nb2,nb1)/(diff_energy^2+1i*gam);
        common_Vzp_o = V_z_o(nb2,nb1)/(diff_energy^2+1i*gam);
        common_Vxn_o = V_x_o(nb2,nb1)/(diff_energy^2-1i*gam);
        common_Vyn_o = V_y_o(nb2,nb1)/(diff_energy^2-1i*gam);
        common_Vzn_o = V_z_o(nb2,nb1)/(diff_energy^2-1i*gam);
        
        
        ch_xy_e(nb1,nb2) = imag(V_x_e(nb1,nb2)*common_Vyp_e + V_x_e(nb1,nb2)*common_Vyn_e);  
        ch_xy_o(nb1,nb2) = imag(V_x_o(nb1,nb2)*common_Vyp_o + V_x_o(nb1,nb2)*common_Vyn_o);  
        ch_yz_e(nb1,nb2) = imag(V_y_e(nb1,nb2)*common_Vzp_e + V_y_e(nb1,nb2)*common_Vzn_e);
        ch_yz_o(nb1,nb2) = imag(V_y_o(nb1,nb2)*common_Vzp_o + V_y_o(nb1,nb2)*common_Vzn_o);  
        ch_zx_e(nb1,nb2) = imag(V_z_e(nb1,nb2)*common_Vxp_e + V_z_e(nb1,nb2)*common_Vxn_e);
        ch_zx_o(nb1,nb2) = imag(V_z_o(nb1,nb2)*common_Vxp_o + V_z_o(nb1,nb2)*common_Vxn_o);  
    end                
end

ch_xy_e = ch_xy_e - ch_xy_e.'; 
ch_xy_o = ch_xy_o - ch_xy_o.'; 
ch_yz_e = ch_yz_e - ch_yz_e.'; 
ch_yz_o = ch_yz_o - ch_yz_o.'; 
ch_zx_e = ch_zx_e - ch_zx_e.'; 
ch_zx_o = ch_zx_o - ch_zx_o.';

Berry(:,1) = sum(ch_xy_e,2);
Berry(:,2) = sum(ch_xy_o,2);
Berry(:,3) = sum(ch_yz_e,2);
Berry(:,4) = sum(ch_yz_o,2);
Berry(:,5) = sum(ch_zx_e,2);
Berry(:,6) = sum(ch_zx_o,2);
