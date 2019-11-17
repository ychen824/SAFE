% Yuan Chen
% PlotMeshSAFE.m

clear all;
close all;

load('WarehouseNetwork.mat');
load('MeshSAFEAverage.mat');

currentFig = figure();
set(currentFig, 'Position', [200 100 550 250]);
subplot(1, 2, 1);
plot(iterationVector, averageErrorNorm(1, :), 'b', 'LineWidth', 1);
hold on;

for i = 2:N
    plot(iterationVector, averageErrorNorm(i, :), 'b', 'LineWidth', 1);
end

title('SAFE Performance', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Iterations (t)', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$||x_n(t) - \theta_{{\cal{I}}_n}^*||_2$', 'interpreter', 'latex', 'FontSize', 12);
xlim([0, numIterations]);

subplot(1, 2, 2);
plot(iterationVector, averageErrorNormR(1, :), 'b', 'LineWidth', 1);
hold on;

for i = 2:N
    plot(iterationVector, averageErrorNormR(i, :), 'b', 'LineWidth', 1);
end

title('CIRFE Performance', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Iterations (t)', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$||x_n(t) - \theta_{{\cal{I}}_n}^*||_2$', 'interpreter', 'latex', 'FontSize', 12);
xlim([0, numIterations]);

%%
grid = [sqrt(m), sqrt(m)];

SAFEWorstEstimate = reshape(SAFEWorstEstimate, grid);
CIWorstEstimate = reshape(CIWorstEstimate, grid);

allError = figure();
set(allError, 'Position', [500, 300, 750, 225]);
subplot(1, 3, 1);
imshow(recoveredImage/255);
title('Ground Truth', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
subplot(1, 3, 2);
imshow(SAFEWorstEstimate/255);
title('SAFE Max Abs. Error', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
subplot(1, 3, 3);
imshow(CIWorstEstimate/255);
title('CIRFE Max Abs. Error', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

%%
allError = figure();
set(allError, 'Position', [800, 100, 550, 300]);
subplot(1, 2, 1);
image(finalSAFEAbsError);
title('SAFE Agent Abs. Error', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Agent Number')
ylabel('Position Abs. Error')

subplot(1, 2, 2);
image(finalCIAbsError);
title('CIRFE Agent Abs. Error', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Agent Number')
ylabel('Position Abs. Error')