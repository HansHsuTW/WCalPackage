%%% Program to convert ftn58sparse.mat to case.vasp %%%
function ftn2cryl(isSftn,ftn58sparse)

if isSftn==0    
    %%% --- ftn58sparse --- %%
    BR    = ftn58sparse.BR;
else
    %%% --- Sftn58sparse --- %%
    BR    = ftn58sparse.BR2D;
end

Nat  = ftn58sparse.Nat;
SNat = repmat([1],Nat,1);

orbID = zeros(ftn58sparse.Nat,5);
for i=1:Nat
    temp       = getfield(ftn58sparse.Ainfo, {i}, 'Atom'); 
    temp_p     = getfield(ftn58sparse.Ainfo, {i}, 'Position');
    Atoms(i,:) = temp; 
    ps(i,:)    = temp_p;
end

sys    = ftn58sparse.System;
vector = BR.*repmat(ftn58sparse.abc,3,1);
ABC    = [1,1,1];

% RT = [cos(pi/2) sin(pi/2) 0;-sin(pi/2) cos(pi/2) 0;0 0 1];
% vector = vector*RT;
% ps = ps*RT;

fid = fopen('Unit.vasp','wt');
fprintf(fid,'%s\n',sys);
fprintf(fid,'%8.6f %8.6f %8.6f\n',ABC);
fprintf(fid,'%8.6f %8.6f %8.6f\n',vector');
for i=1:Nat
    fprintf(fid,'%s\t',Atoms(i,:));
end
fprintf(fid,'\n');
for i=1:Nat
    fprintf(fid,'%i\t',SNat(i));
end
fprintf(fid,'\n%s\n','Direct');
fprintf(fid,'%12.8f %12.8f %12.8f\n',ps');
fclose(fid);

end
    
%% --- Sftn58sparse --- %%
% Nat       = size(superatomps,1);
% [Spe,~,~] = unique(superatomps(:,5));
% Type      = length(Spe);
% Atom_unit = sortrows(superatomps,5);
% ps        = Atom_unit(:,9:11);
% vector    = [vector(1,:);vector(2,:);vector(3,:)*nlayer];
% ps        = [ps(:,1) ps(:,2) ps(:,3)/nlayer];
% 
% for i=1:Type
%     SNat(i) = length(find(superatomps(:,5)==Spe(i)));
% end
% 
% fid = fopen('Slab.vasp','wt');
% fprintf(fid,'%s\n','Slab');
% fprintf(fid,'%8.6f %8.6f %8.6f\n',ABC);
% fprintf(fid,'%8.6f %8.6f %8.6f\n',vector');
% for i=1:Type
%     fprintf(fid,'%s\t',Atoms(i,:));
% end
% fprintf(fid,'\n');
% for i=1:Type
%     fprintf(fid,'%i\t',SNat(i));
% end
% fprintf(fid,'\n%s\n','Direct');
% fprintf(fid,'%12.8f %12.8f %12.8f\n',ps');
% fclose(fid);
% 
% end