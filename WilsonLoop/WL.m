%%% Function to carry the Wilson Loop calculation %%%
function WL(ftn58sparse,ocnorb,nkx,nky)

kx   = linspace(0,2*pi,nkx+1);
kx   = kx(1:nkx);
ky   = linspace(0*pi,1*pi,nky);
kz   = pi;

%%-- Load data --%% 
Norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
tt   = ftn58sparse.tt;
dd   = ftn58sparse.dd;

%%-- Acutal Procedure --%%
D     = zeros(ocnorb,ocnorb);
theta = zeros(nky,ocnorb);

tic
%--- 1st run for [001]
for ikz=0:1:1
    kz = (ikz-1)*pi;
    filename = ['001_' num2str(ikz) '.mat'];
parfor iky=1:nky
    tic
    time_init=tic;
    
    %%% Prepare F(N-1,0) %%%
    Fn      = eye(ocnorb);
    kpoints = [kx(1) ky(iky) kz];
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2;
    [V,~]   = eig(HH);
    evec0   = V(1:end,1:ocnorb);
    eveco   = V(1:end,1:ocnorb);
    
    for ikx=2:nkx
        kpoints = [kx(ikx) ky(iky) kz];
        Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
        HH      = full(Hsparse);
        HH      = (HH+HH')/2;
        [V,~]   = eig(HH);
        evecn   = V(1:end,1:ocnorb);
        
        %%% Prepare Fmn %%%
        Fn    = Fn*(eveco'*evecn);
        eveco = evecn;
    end

    F0 = eveco'*evec0;
    Fn = Fn*F0;
    
    %%% Calculate theta %%%
    raw_theta(iky,:) = angle(eig(Fn))/pi;    
    theta(iky,:) = sortrows(angle(eig(Fn))/pi);
        
%     if mod(iky,10)==0
%         fprintf('%3i/%i: %.3fs [%.2fm]\n',iky,nky,toc,toc(time_init)/60);
%     end
end
    save(filename, 'ky', 'theta', 'ocnorb');
end
disp('[001]--finished');

%--- 2nd runt for [010]
for ikz=0:1:1
    kz = (ikz-1)*pi;
    filename = ['010_' num2str(ikz) '.mat'];
parfor iky=1:nky
    tic
    time_init=tic;
    
    %%% Prepare F(N-1,0) %%%
    Fn      = eye(ocnorb);
    kpoints = [kx(1) kz ky(iky)];
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2;
    [V,~]   = eig(HH);
    evec0   = V(1:end,1:ocnorb);
    eveco   = V(1:end,1:ocnorb);
    
    for ikx=2:nkx
        kpoints = [kx(ikx) kz ky(iky)];
        Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
        HH      = full(Hsparse);
        HH      = (HH+HH')/2;
        [V,~]   = eig(HH);
        evecn   = V(1:end,1:ocnorb);
        
        %%% Prepare Fmn %%%
        Fn    = Fn*(eveco'*evecn);
        eveco = evecn;
    end

    F0 = eveco'*evec0;
    Fn = Fn*F0;
    
    %%% Calculate theta %%%
    raw_theta(iky,:) = angle(eig(Fn))/pi;    
    theta(iky,:) = sortrows(angle(eig(Fn))/pi);
        
%     if mod(iky,10)==0
%         fprintf('%3i/%i: %.3fs [%.2fm]\n',iky,nky,toc,toc(time_init)/60);
%     end
end
    save(filename, 'ky', 'theta', 'ocnorb');
end
disp('[010]--finished');

%--- 3rd run for [100]
for ikz=0:1:1
    kz = (ikz-1)*pi;
    filename = ['100_' num2str(ikz) '.mat'];
parfor iky=1:nky
    tic
    time_init=tic;
    
    %%% Prepare F(N-1,0) %%%
    Fn      = eye(ocnorb);
    kpoints = [kz kx(1) ky(iky)];
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2;
    [V,~]   = eig(HH);
    evec0   = V(1:end,1:ocnorb);
    eveco   = V(1:end,1:ocnorb);
    
    for ikx=2:nkx
        kpoints = [kz kx(ikx) ky(iky)];
        Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,Norb,Norb);
        HH      = full(Hsparse);
        HH      = (HH+HH')/2;
        [V,~]   = eig(HH);
        evecn   = V(1:end,1:ocnorb);
        
        %%% Prepare Fmn %%%
        Fn    = Fn*(eveco'*evecn);
        eveco = evecn;
    end

    F0 = eveco'*evec0;
    Fn = Fn*F0;
    
    %%% Calculate theta %%%
    raw_theta(iky,:) = angle(eig(Fn))/pi;    
    theta(iky,:) = sortrows(angle(eig(Fn))/pi);
        
%     if mod(iky,10)==0
%         fprintf('%3i/%i: %.3fs [%.2fm]\n',iky,nky,toc,toc(time_init)/60);
%     end
end
    save(filename, 'ky', 'theta', 'ocnorb');
end
disp('[100]--finished');