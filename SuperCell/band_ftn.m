%%%          Band Structure Plot for ftn58sparse         %%%
%%% The spin polariztion is implemented in this function %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%
clear all

E_range = [-2 2];
Ef      = -2.0960;
isSP    = 0;

%% Initial info. %%
InFile = 'super_ftn58sparse.mat';

instruct    = load(InFile);
if isfield(instruct,'ftn58sparse')
    ftn58sparse = instruct.ftn58sparse;
elseif isfield(instruct,'SPftn58sparse')
    ftn58sparse = instruct.SPftn58sparse;
else
    ftn58sparse = instruct.Sftn58sparse;
end

norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR;
Sz   = [1 0;0 -1];
% Sz   = [0 -1i;1i 0];

%% Kpoints %%%
nk = 50;
p1 = [linspace(1/3,0,nk+1)' linspace(1/3,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p2 = [linspace(0.0,0.5,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p3 = [linspace(0.5,1/3,nk)' linspace(0.0,1/3,nk)' linspace(0.0,0.0,nk)'];
    
list    = [1 51 101 150];
kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:end,:)]*2*pi;

% p1 = [linspace(0.0,0.5,nk)' linspace(0.5,0.5,nk)' linspace(0.0,0.0,nk)'];
% p2 = [linspace(0.0,-0.5,nk)' linspace(-0.5,-0.5,nk)' linspace(0.0,0.0,nk)'];
%     
% list    = [1 151 301];
% kpoints = [p1(1:nk,:);p2(1:nk,:)]*2*pi;

% p1 = [linspace(0.5,0.0,nk)' linspace(0.0,0.0,nk)' linspace(0.0,0.0,nk)'];
% p2 = [linspace(0.0,0.0,nk)' linspace(0.0,0.5,nk)' linspace(0.0,0.0,nk)'];
% p3 = [linspace(0.0,0.5,nk)' linspace(0.5,0.5,nk)' linspace(0.0,0.0,nk)'];
%     
% list    = [1 10 20 30];
% kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:)]*2*pi;
  
%% Time-reversal operator %%
sy = [0 -1i;1i 0];
T  = 1i*kron(sy,eye(2));

%% Eigenvalue and Eigenvector calculations %%%
tic
eigvec = cell(size(kpoints,1),1);
nks    = length(kpoints);
Ek     = zeros(nks,norb);
SPz    = zeros(nks,norb);
for ik=1:nks
    time_init=tic;
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2 ;
    [vec, Etemp] = eig(HH);
     
    Ham{ik,1}    = HH;
    eigvec{ik,1} = vec;

%     A(ik*norb-norb+1:ik*norb,1) = ik;
%     A(ik*norb-norb+1:ik*norb,2) = diag(Etemp);
    
    Ek(ik,:)  = diag(Etemp);
    SSPz      = vec'*kron(Sz,eye(norb/2))*vec;
    SSPz      = (SSPz+SSPz')/2;
    SPz(ik,:) = diag(SSPz); 
    
    if mod(ik,10)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
    end
    
end
toc

[~, Emin] = min(abs(Ek(1,:)-Ef));
disp(Emin);

%% Plotting %%%
h = figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');

if isSP==1
    ncolor   = 1e4;
    spmap_p  = [linspace(1,1,ncolor+1)' linspace(0,1,ncolor+1)' linspace(0,1,ncolor+1)'];
    spmap_n  = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)'];
    spmap    = [spmap_p(1:ncolor,:);spmap_n];

    for ii=1:norb
        cplot(1:nks,Ek(:,ii)-Ef,SPz(:,ii),'-','LineWidth',1);
        colormap(spmap);
        hold on
    end
    colorbar('TickLabelInterpreter','LaTex','FontSize',18);
else
    for ii=1:norb
        plot(1:nks,Ek(:,ii)-Ef,'b-','LineWidth',1.5);
        hold on
    end
end

for il = 2:size(list,2)-1
line('XData', [list(il) list(il)], 'YData', [E_range(1) E_range(2)], 'LineStyle', '-', ...
    'LineWidth', 0.1, 'Color','k');
end
line('XData', [0 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');

axis([1 nks E_range(1) E_range(2)]);
ax = gca;
ax.TickDir    = 'out';
ax.FontSize   = 26;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = list(:);
ax.YTick      = [E_range(1) E_range(2)];
ax.XTickLabel = {'\bf{K}' '\bf{$\Gamma$}' '\bf{M}' '\bf{K}'};
% ax.XTickLabel = {};
ax.LineWidth  = 1.2;
ax.TickLabelInterpreter='latex';

% save eigvec.mat Ham eigvec