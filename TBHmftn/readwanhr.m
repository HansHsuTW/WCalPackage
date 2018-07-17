function [ftn58 , ainfo]=readwanniercorrect(s,zero)
% see hamiltonian.F90
% s='g1wan_hr.dat';
% zero=0 ; %energy shift
fid = fopen(s,'r');

ainfo.title = fgetl(fid);
ainfo.norb  = str2num(fgetl(fid));
ainfo.nrpts = str2num(fgetl(fid));
ainfo.ndegen=fscanf(fid,'%i',[1 ainfo.nrpts]);
norb        = ainfo.norb;
nbondtemp   = ainfo.norb*ainfo.norb*ainfo.nrpts;

data     = fscanf(fid,'%f',[7 nbondtemp])';
ainfo.vec= data(1:norb*norb:end,1:3);
tfactor  = repmat(ainfo.ndegen,ainfo.norb*ainfo.norb,1);
data(:,6)= data(:,6)./tfactor(:);
data(:,7)= data(:,7)./tfactor(:);
ibondselect= find(data(:,4)<=data(:,5));
nbond      = length(ibondselect);
%            ftn58=[[norb nbond 0 0 0 0 0]; [(1:nbond)' data(:,4:5) data(:,6)+ i*data(:,7) data(:,1:3)]];
ftn58=[[norb nbond 0 0 0 0 0]; [(1:nbond)' data(ibondselect,4:5) data(ibondselect,6)+ ...
       1i*data(ibondselect,7) data(ibondselect,1:3)]];
           
return
            
