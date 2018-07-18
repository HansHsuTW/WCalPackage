%%% Program for calculating the mirror Chern number (Mxy) in SnTe        %%%
%%% -------------------------------------------------------------------- %%%
clear all

%% Initial setup for i, j, and k  %%%       
Nk    = [101,101];    
isRT  = 0;
ival  = 3;
np    = 6;

%% K-mesh setup %%%
load ftn58sparse_SnTe.mat

N1 = Nk(1); N2 = Nk(2);
p1 = linspace(0,2*pi,N1+1); 
p2 = linspace(0,2*pi,N2+1); 
p1 = p1(1:N1); 
p2 = p2(1:N2); 

kk = 1;
for ii=1:N1
    for jj=1:N2
            kpts(kk,:) = [p1(ii) p1(ii) p1(ii)+p2(jj)] ;
            kk = kk + 1;
    end
end
nks = length(kpts);

%% --- Berry phase calculation --- %%%
norb = ftn58sparse.norb;

if isRT ==0
    c = parcluster('local');
    c.NumWorkers = np;
    parpool(c, c.NumWorkers);
    
    Berry = zeros(nks,norb/2,6);
    Ek    = zeros(nks,norb);
    tic
    parfor ik=1:nks
            kpoints          = kpts(ik,:);
            [Btmp, Etmp]     = MirrorChern_Mxy(kpoints,ftn58sparse);
            Berry(ik,:,:) = [Btmp];
            Ek(ik,:)      = Etmp;               
    end   
    save MBerry.mat Berry Ek
else
    load MBerry_xy.mat
end

%% --- QSHE at various energies --- %%%
cquantum = 1; % the unit becomes (e^2/hbar)
kB       = 8.6173324E-5;
BR       = ftn58sparse.BR;
abc      = ftn58sparse.abc;
Emu      = [-1,1];
nE       = 2;
Sk       = [5 6];  

T        = [BR(:,1)*abc(1) BR(:,2)*abc(2) BR(:,3)*abc(3)];
G        = inv(T);
Ainv     = norm(cross(G(:,1)+G(:,2)+G(:,3),G(:,3)));
Dens     = 2*pi*Ainv/N1/N2;             % density, return in unit of cm      
COF      = Dens*cquantum;
mu       = linspace(Emu(1),Emu(2),nE);  % chemical potnetial
cond_ss  = zeros(size(mu,2),length(Sk));         

for imu = 1:length(mu)
    u  = mu(imu); 
    for ii=1:size(cond_ss,2)
        Berry_new = Berry(:,1:ival,:);
        cond_ss(imu,ii) = COF*sum(sum(Berry_new(:,:,Sk(ii))));
    end
end

%% --- Plot Figure --- %%%
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');
clf
for ii=1:size(cond_ss,2)
    plot(mu,cond_ss(:,ii),'Linewidth',1.0);
    hold on
end
line('XData',[mu(1) mu(end)], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
ls = legend('\bf{$\sigma^z_{yx}$}','\bf{$\sigma^y_{zx}$}');
set(ls,'interpreter','LaTex','FontSize',20);
set(gca,'Ticklength',[0.05 0]);
set(gca,'Fontsize',15);
xlabel('\bf{$\mu$ (eV)}','interpreter','LaTex');
ylabel('\bf{$\sigma$ ($\frac{e}{\hbar} \Omega^{-1}cm^{-1}$)}','interpreter','LaTex');
ax = gca;
ax.TickLabelInterpreter='latex';
ax.FontSize = 20;