% Yuan Chen
% PlotMeshNetwork.m

clearvars;
close all;

load('PlotNetwork.mat');

figure();
hold on;

for i = 0:sectorMeshSize
    line([0, sectorMeshSize], -[i, i], 'LineStyle', '-', 'LineWidth', .1, 'color', .65*[1, 1, 1]);
    line([i, i], -[0, sectorMeshSize], 'LineStyle', '-', 'LineWidth', .1, 'color', .65*[1, 1, 1]);
end

for i = 1:N
    currentLocation = LocationArray{i};
    
    for j = i+1:N
        if (Adjacency(i, j) == 1)
            neighborLocation = LocationArray{j};
            line([currentLocation(1), neighborLocation(1)], -[currentLocation(2), neighborLocation(2)], 'LineWidth', 0.25);
        end
    end
    endc


for i = 1:N
    if (adversaryIndicator(i) == 0)
        currentLocation = LocationArray{i};
        plot(currentLocation(1), -currentLocation(2), 'ok', 'MarkerSize', 7, 'MarkerFaceColor', [.8, .8, .8]); 
    end
end

for i = 1:N
    if (adversaryIndicator(i) == 1)
        currentLocation = LocationArray{i};
        plot(currentLocation(1), -currentLocation(2), 'dk', 'MarkerSize', 6, 'MarkerFaceColor', [.8, .1, .05]);
    end
end

axisFactor = GridDimension;
title(sprintf('Network of %d Agents', N), 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'XTick', [], 'YTick', [], 'XColor', [1 1 1], 'YColor', [1 1 1]);
axis(axisFactor*[-.1 1.1 -1.1 .1]);