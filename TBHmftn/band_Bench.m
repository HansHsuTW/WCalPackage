%%%          Wannier TB-model Benchmark                  %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%
clear all

E_range = [-3 3];
Ef      = 4.3391;
load SnTesoc.mat

%% Initial info. %%
load ftn58sparse
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;

%% Kpoints %%%
load KP.mat
kpoints = kpoints(:,1:3);

%% Calculate band structrue %%%
tic
eigvec = cell(size(kpoints,1),1);
for ik=1:size(kpoints,1)
     kcolumnvec = kpoints(ik,:)'*2*pi;
     Hsparse    = sparse(ii,jj,exp(1i*dd*kcolumnvec).*tt,norb,norb);
     HH         = full(Hsparse);
     HH         = (HH+HH')/2;
     [eigvectemp, Etemp]=eig(HH);

    A(ik*norb-norb+1:ik*norb,1)=ik;
    A(ik*norb-norb+1:ik*norb,2)=diag(Etemp);
    
end
toc

%% Plotting %%%
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');

nks  = pW2k.nks;
nbnd = pW2k.nbnd;
Eng  = pW2k.Eofk;
 
for ib = 1:nbnd
    h1 = plot(1:nks,Eng(ib,:),'r-','LineWidth',1.5);
    hold on
end
      
h2 = plot(A(:,1),A(:,2)-Ef,'b.','MarkerSize',15);
hold on

line('XData', [0 A(end,1)], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');

%%% Plotting Details %%%
axis([A(1,1) A(end,1) E_range(1) E_range(2)]);
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');
ls = legend([h1 h2],'\bf{DFT}','\bf{W90}','Location','best');
set(ls,'interpreter','LaTex','FontSize',22);

ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.TickLabelInterpreter='latex';
ax.XTickLabel = {};
ax.LineWidth = 0.5;