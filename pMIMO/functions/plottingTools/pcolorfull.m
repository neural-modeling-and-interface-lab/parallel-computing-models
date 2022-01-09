function pcolorfull(C)
% Show the full matrix with pcolor

[a,b] = size(C);
C2 = vertcat(C,zeros(1,b));
C2 = horzcat(C2,zeros((a+1),1));
h = pcolor(C2);
shading flat; 
% shading faceted ;
set(gca, 'XTick', [], 'YTick', []);
% set(gca, 'TickDir', 'out', 'XTick', [1:b], 'YTick', [1:a]);
% set(gca, 'TickDir', 'out', 'XTick', [], 'YTick', [1:a]);