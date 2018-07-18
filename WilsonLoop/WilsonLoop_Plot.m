%%% Script for plotting the Wilson Loop result %%%
%%% ------------------------------------------ %%%
clear all;

load 100_0.mat
%%-- Plotting Detail --%%
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off','name','W1',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');

for iorb=1:ocnorb
    plot(ky(1:end)/pi,theta(1:end,iorb),'r.');
    hold on
end

axis([ky(1)/pi ky(end)/pi -1 1]);
ylabel('\bf{$\theta$($\pi$)}','FontSize',26,'Interpreter','Latex');
xlabel('\bf{$k_{y}$($\pi$)}','FontSize',26,'Interpreter','Latex');

ax = gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.TickLabelInterpreter='latex';