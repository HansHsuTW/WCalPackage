%%% Program to convert ftn58sparse.mat to case.vasp %%%
function Super2vasp(atomps,ftn58sparse,BR)

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
for i=1:Type
    temp       = getfield(ftn58sparse.Ainfo, {Spe(i)}, 'Atom');  
    Atoms(i,:) = temp; 
end

sys    = ftn58sparse.System;
vector = BR.*repmat(ftn58sparse.abc,3,1);
ABC    = 1;

fid = fopen('Unit.vasp','wt');
fprintf(fid,'%s\n',sys);
fprintf(fid,'%8.6f \n',ABC);
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