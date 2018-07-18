%% Plotting %%
Ef = 5.3384;

ncolor = 1e4;
map    = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)'];
 
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');
pcolor(KK(2:end-1,2:end),EE(2:end-1,2:end)-Ef,...
       ((As(2:end-1,2:end)*1.0+Ab(2:end-1,2:end)*10))/pi),shading interp;
colormap(map);
caxis([0 400]);
box on

ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');

ax = gca;
ax.FontSize   = 20;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTickLabel = { };
ax.LineWidth  = 2;
ax.TickLabelInterpreter='latex';