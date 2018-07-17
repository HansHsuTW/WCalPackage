%%% Program to convert ftn58sparse.mat to case.vasp %%%
function ftn2vasp(atomps,superatomps,ftn58sparse,BR,nlayer)

%% --- Unit cell --- %%
Nat       = size(atomps,1);
[Spe,~,~] = unique(atomps(:,5));
Type      = length(Spe);
Atom_unit = sortrows(atomps,5);
ps        = Atom_unit(:,9:11);

for i=1:Type
    SNat(i) = length(find(atomps(:,5)==Spe(i)));
end

orbID = zeros(ftn58sparse.Nat,5);
for i=1:ftn58sparse.Nat
    temp       = getfield(ftn58sparse.Ainfo, {i}, 'Atom');  
    Atoms(i,:) = temp; 
end

sys    = ftn58sparse.System;
vector = BR.*repmat(ftn58sparse.abc,3,1);
ABC    = [1,1,1];

fid = fopen('Unit.vasp','wt');
fprintf(fid,'%s\n',sys);
fprintf(fid,'%8.6f %8.6f %8.6f\n',ABC);
fprintf(fid,'%8.6f %8.6f %8.6f\n',vector');
for i=1:Type
    fprintf(fid,'%s\t',Atoms(i,:));
end
fprintf(fid,'\n');
for i=1:Type
    fprintf(fid,'%i\t',SNat(i));
end
fprintf(fid,'\n%s\n','Direct');
fprintf(fid,'%12.8f %12.8f %12.8f\n',ps');
fclose(fid);

%% --- Slab cell --- %%
Nat       = size(superatomps,1);
[Spe,~,~] = unique(superatomps(:,5));
Type      = length(Spe);
Atom_unit = sortrows(superatomps,5);
ps        = Atom_unit(:,9:11);
vector    = [vector(1,:);vector(2,:);vector(3,:)*nlayer];
ps        = [ps(:,1) ps(:,2) ps(:,3)/nlayer];

for i=1:Type
    SNat(i) = length(find(superatomps(:,5)==Spe(i)));
end

fid = fopen('Slab.vasp','wt');
fprintf(fid,'%s\n','Slab');
fprintf(fid,'%8.6f %8.6f %8.6f\n',ABC);
fprintf(fid,'%8.6f %8.6f %8.6f\n',vector');
for i=1:Type
    fprintf(fid,'%s\t',Atoms(i,:));
end
fprintf(fid,'\n');
for i=1:Type
    fprintf(fid,'%i\t',SNat(i));
end
fprintf(fid,'\n%s\n','Direct');
fprintf(fid,'%12.8f %12.8f %12.8f\n',ps');
fclose(fid);