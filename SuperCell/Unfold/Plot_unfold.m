%%% plot unfolding band structure %%%
%%% by Hans 4 Aug 2018            %%%
clear all

%% inupt %%
Ef      = -2.0960;
E_range = [-3 3]; 

klabel = {'\bf{K}' '\bf{$\Gamma$}' '\bf{M}' '\bf{K}'};
kid    = [1 51 101 150];

%% Plot %%
load unfold.mat
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','a4','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','landscape',...
       'color','w');

for ii=1:norbss
    plotweight = Weight(1:nks,ii)*20+1e-5;
    h1 = plot(1:nks,Ek(1:nks,ii)-Ef,'LineStyle','-','Color',[0.6 0.6 0.6],'LineWidth',1);
    hold on
    h2 = scatter(1:nks,Ek(1:nks,ii)-Ef,plotweight, 'filled','MarkerFaceColor','b');
    hold on
end

box on

%%% Plotting Details %%%
line('XData', [1 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');

for il = 1:size(kid,2)
line('XData', [kid(il) kid(il)], 'YData', [E_range(1) E_range(2)], 'LineStyle', '-', ...
    'LineWidth', 0.1, 'Color','k');
end

axis([1 nks E_range(1) E_range(2)]);
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');
ls = legend([h1 h2],'supercell','unfolded supercell','Location','best');
% legend('boxon');
set(ls,'interpreter','LaTex','FontSize',20);
set(ls,'Color','w');

ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = kid;
ax.TickLabelInterpreter='latex';
ax.XTickLabel = klabel;
ax.LineWidth = 0.5;