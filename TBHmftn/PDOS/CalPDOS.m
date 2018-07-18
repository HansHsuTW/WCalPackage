%%%          Project orbit weight on band structure      %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%
clear all 

load eigvec.mat

%% Input
Ef      = 5.1534;
iorb    = [4:6 10:12]; 
E_range = [-2 2];

%% Actual Procedure
norb = size(Ek,2);
nks  = size(Ek,1);
kB   = 8.6173324E-5;

orbweig = zeros(nks,norb);
for ik=1:nks
    wave{ik} = round(conj(eigvec{ik}).*eigvec{ik},1);
    orbweig(ik,:) = sum(wave{ik}(iorb,:),1);
end

%% Plot
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
           'papertype','usletter','numbertitle','off','name','kx_kz',...
           'PaperPositionMode','manual','paperorientation','portrait',...
           'color','w');
       
ncolor   = 1e4;
map = [linspace(1,1,ncolor+1)' linspace(1,0,ncolor+1)' linspace(1,0,ncolor+1)'];
       
 for ib=1:norb
     cplot(1:nks,Ek(:,ib)-Ef,orbweig(:,ib),'-','LineWidth',3);
     colormap(map);
     hold on
 end

for ii=1:norb
    plot(1:nks,Ek(:,ii)-Ef,'k--','LineWidth',0.5);
    hold on
end 

line('XData', [0 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');

for il = 2:size(list,2)-1
line('XData', [list(il) list(il)], 'YData', [E_range(1) E_range(2)], 'LineStyle', '-', ...
    'LineWidth', 0.1, 'Color','k');
end
axis([1 nks E_range(1) E_range(2)]);
ax = gca;
ax.FontSize   = 24;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = list(:);
ax.XTickLabel = {'\bf{K}' '\bf{L}' '\bf{$\Gamma$}' '\bf{X}' '\bf{W}' '\bf{$\Gamma$}' '\bf{K}'};
ax.LineWidth  = 0.5;
ax.TickLabelInterpreter='latex';