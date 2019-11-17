% Yuan Chen
% MeshSAFEMultiTrial.m

clearvars;
close all;

set(0,'defaulttextinterpreter','none')

%%
load('WarehouseNetwork.mat'); %Load predetermined network

thetaImage = imread('Baboon.jpg');
thetaImage = rgb2gray(thetaImage);
thetaImage = imresize(thetaImage, [sqrt(m), sqrt(m)]);

thetaStar = double(reshape(thetaImage, m, 1));

recoveredImage = reshape(thetaStar, [sqrt(m), sqrt(m)]);

figure();
imshow(recoveredImage/255);

[P, ~] = size(DH);

%%
totalMeasureDim = sum(numVisible);
normAddMatrix = zeros(N, totalMeasureDim);

currentPos = 1;
for i = 1:N
    currentDim = numVisible(i);
    normAddMatrix(i, currentPos:currentPos + currentDim - 1) = 1;
    currentPos = currentPos + currentDim;
end

normAddMatrix = sparse(normAddMatrix);

%% Parameter Selection for SIU algorithm
%securityIndex = 1.6/N;
%numberAttacked = max(round(securityIndex*N)-1, 0);
evalues = sort(eig(L));

%b = 2/(evalues(2)+evalues(N));
b = 1/evalues(N);
a = 1;
%a = 1/(1 - 2*securityIndex);
tau1 = .26;
tau2 = .001;
tauGamma = 0.25;

capGamma = 40;

kappa1 = 1 + sqrt(N);

L2 = evalues(2);
%a = (0.5*b*L2)/(kappa1);
%b = b/10;
%a = .1*b;
numErrorPoints = 126;
strideLength = 1; % The number of iterations between error points

numIterations = 1 + (numErrorPoints-1)*strideLength;
fprintf('Total Number of Iterations: %d\n', numIterations);

% Evaluate threshold evolution
gammaHistory = zeros(1, numErrorPoints);
gamma1History = zeros(1, numErrorPoints);
gamma2History = zeros(1, numErrorPoints);

gamma = capGamma;
gammaHistory(1) = gamma;
% Set Threshold Gamma

for t = 2:numIterations
    gammaHistory(t) = capGamma/((t-1)^(tauGamma));
end

%% Simulate Distributed Algorithm

attackVector = sparse(255 - thetaStar);
errorNormAddMatrix = sparse(kron(speye(N), ones(1, m)));
trueParam = sparse(BigP*kron(ones(N, 1), thetaStar));
calH = cat(1, HArray{:});
trueMeasurement = calH*thetaStar;
calNH = cat(1, normalHArray{:});
calAH = calH - calNH;
measurementAttack = 2*calAH*attackVector;
sparseIM = speye(m);
sparseIN = speye(N);
currentAttackIndex = 1;


numTrials = 100;

averageErrorNorm = zeros(N, numErrorPoints);
averageErrorNormR = zeros(N, numErrorPoints);

averageEstimate = sparse(N*m, 1);
averageEstimateR = sparse(N*m, 1);

varTheta = var(thetaStar);

noiseSigma = 50;
for currentTrial = 1:numTrials
    taMeasurement = zeros(totalMeasureDim, 1);

    EstimateErrorNormR = zeros(N, numErrorPoints);
    EstimateErrorNorm = zeros(N, numErrorPoints);

    EstimateHistory = sparse(N*m, 1);
    EstimateHistoryR = sparse(N*m, 1);

    stackedError = EstimateHistory - trueParam;
    stackedErrorR = EstimateHistoryR - trueParam;

    currentErrorNormSquared = errorNormAddMatrix*(stackedError.^2);
    currentErrorNormSquaredR = errorNormAddMatrix*(stackedErrorR.^2);

    EstimateErrorNorm(:, 1) = sqrt(currentErrorNormSquared);
    EstimateErrorNormR(:, 1) = sqrt(currentErrorNormSquaredR);

    for t = 2:numIterations
        %currentAttackSet = [currentAttackIndex];

        alpha = a/((t-1)^(tau1));
        beta = b/((t-1)^tau2);
        %dynamicsMatrixFirstPart = speye(N) - beta*L;
        dynamicsMatrixFirstPart = speye(N*m) - beta*censoredL;
        gamma = capGamma/((t-1)^(tauGamma));

        currentNoise = noiseSigma*randn(totalMeasureDim, 1);

        %randMeasurementAttack = randAttackVector*attackVector;
        %randMeasurement = trueMeasurement+randMeasurementAttack + currentNoise;
        randMeasurement = trueMeasurement + measurementAttack + currentNoise;

        taMeasurement = ((t-2)/(t-1))*taMeasurement + 1/(t-1) * randMeasurement;

        fixedInnovation = DH*EstimateHistory - taMeasurement;
        randInnovation = DH*EstimateHistoryR - randMeasurement;

        %fixedZ = fixedInnovation.^2;
        %randZ = randInnovation.^2;

        %fixedZ = sqrt(normAddMatrix*fixedZ);
        %randZ = sqrt(normAddMatrix*randZ);

        fixedZ = abs(fixedInnovation);

        fixedK = min(gamma./fixedZ, ones(P, 1));
        %randK = min(gamma./randZ, ones(N, 1));

        %fixedKMatrix = sparse(diag(fixedK));
        %randKMatrix = sparse(diag(randK));

        fixedKMatrix = spdiags(fixedK(:), 0, P, P);

        %fixedDynamics = sparse(kron(dynamicsMatrixFirstPart - alpha*fixedKMatrix, sparseIM));
        %randDynamics = sparse(kron(dynamicsMatrixFirstPart - alpha*sparseIN, sparseIM));

        %fixedDynamics = sparse(kron(dynamicsMatrixFirstPart, sparseIM) - alpha*DH'*fixedKMatrix*DH);
        %randDynamics = sparse(kron(dynamicsMatrixFirstPart, sparseIM) - alpha*DH'*DH);

        fixedDynamics = dynamicsMatrixFirstPart - sparse(alpha*DH'*fixedKMatrix*DH);
        randDynamics = dynamicsMatrixFirstPart - sparse(alpha*DH'*DH);

        EstimateHistory = fixedDynamics*EstimateHistory + alpha*sparse(DH'*fixedKMatrix*taMeasurement);
        EstimateHistoryR = randDynamics*EstimateHistoryR + alpha*sparse(DH'*taMeasurement);

        if (mod(t-1, strideLength) == 0)
            stackedError = EstimateHistory - trueParam;
            stackedErrorR = EstimateHistoryR - trueParam;

            currentErrorNormSquared = errorNormAddMatrix*(stackedError.^2);
            currentErrorNormSquaredR = errorNormAddMatrix*(stackedErrorR.^2);

            EstimateErrorNorm(:, floor((t-1)/strideLength)+1) = sqrt(currentErrorNormSquared);
            EstimateErrorNormR(:, floor((t-1)/strideLength)+1) = sqrt(currentErrorNormSquaredR);
        end
        %gamma = r*gamma;


        if (mod(t, 10) == 0)
            fprintf('Iteration %d complete \n', t);
        end
        %currentAttackIndex = mod(currentAttackIndex, N)+1;
        %fprintf('\n');
    end
    averageErrorNorm = averageErrorNorm + EstimateErrorNorm;
    averageErrorNormR = averageErrorNormR + EstimateErrorNormR;
    averageEstimate = averageEstimate + EstimateHistory;
    averageEstimateR = averageEstimateR + EstimateHistoryR;
    fprintf('=== Trial %d complete === \n', currentTrial);
end

averageErrorNorm = averageErrorNorm/numTrials;
averageErrorNormR = averageErrorNormR/numTrials;
averageEstimate = averageEstimate/numTrials;
averageEstimateR = averageEstimateR/numTrials;
%% Plot results
tempIterationVector = 1:numErrorPoints;
iterationVector = (tempIterationVector-1)*strideLength;
currentFig = figure();
set(currentFig, 'Position', [200 100 550 250]);
subplot(1, 2, 1);
%line1 = plot(iterationVector, errorBound, 'r--', 'LineWidth', 2.5);
%hold on;
plot(iterationVector, averageErrorNorm(1, :), 'b', 'LineWidth', 1);
hold on;

for i = 2:N
    plot(iterationVector, averageErrorNorm(i, :), 'b', 'LineWidth', 1);
end

title('SAFE Performance', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Iterations (t)', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$||x_n(t) - \theta_{{\cal{I}}_n}^*||_2$', 'interpreter', 'latex', 'FontSize', 12);
%legend([line1, line2], {'Error Bound', 'Agent Errors'}, 'Location', 'Southwest');
xlim([0, numIterations]);
%ylim([0, 10^2]);
%export_fig FixedAttackHigh.png -transparent -m2 -painters

subplot(1, 2, 2);
%line1 = plot(iterationVector, errorBound, 'r--', 'LineWidth', 2.5);
%hold on;
plot(iterationVector, averageErrorNormR(1, :), 'b', 'LineWidth', 1);
hold on;

for i = 2:N
    plot(iterationVector, averageErrorNormR(i, :), 'b', 'LineWidth', 1);
end

title('CIRFE Performance', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Iterations (t)', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$||x_n(t) - \theta_{{\cal{I}}_n}^*||_2$', 'interpreter', 'latex', 'FontSize', 12);
%legend([line1, line2], {'Error Bound', 'Agent Errors'}, 'Location', 'Southwest');
xlim([0, numIterations]);
%ylim([0, 10^2]);

%export_fig SAFEWarehouse.png -transparent -m2 -painters

save('MeshSAFEAverage.mat', 'numIterations', 'iterationVector', 'averageErrorNormR', 'averageErrorNorm', 'N');

%%
[totalDim, ~] = size(EstimateHistory);

finalSAFEEstimates = averageEstimate;
finalCIEstimates = averageEstimateR;

finalSAFEEstimatesReshape = reshape(finalSAFEEstimates, [m, N]);
finalCIEstimatesReshape = reshape(finalCIEstimates, [m, N]);

%flatP = cat(2, PArray{:});

thetaStarReshape = reshape(trueParam, [m, N]);

%thetaStarReshape = kron(ones(1, N), thetaStar);
%thetaStarReshape = flatP*kron(ones(N, 1), thetaStar);

finalSAFEAbsError = abs(finalSAFEEstimatesReshape - thetaStarReshape);
finalCIAbsError = abs(finalCIEstimatesReshape - thetaStarReshape);

grid = [sqrt(m), sqrt(m)];

[SAFEMax, SAFEMaxArgs] = max(finalSAFEAbsError, [], 2);
[CIMax, CIMaxArgs] = max(finalCIAbsError, [], 2);

SAFEAgentMax = reshape(max(SAFEMax, [], 2), grid);
CIAgentMax = reshape(max(CIMax, [], 2), grid);

SAFEAgentMean = reshape(mean(finalSAFEAbsError, 2), grid);
CIAgentMean = reshape(mean(finalCIAbsError, 2), grid);

SAFEAgentMeanEstimate = full(reshape(mean(finalSAFEEstimatesReshape, 2), grid));
CIAgentMeanEstimate = full(reshape(mean(finalCIEstimatesReshape, 2), grid));

SAFEWorstEstimate = zeros(m, 1);
CIWorstEstimate = zeros(m, 1);

for i = 1:m
    SAFEWorstEstimate(i) = finalSAFEEstimatesReshape(i, SAFEMaxArgs(i));
    CIWorstEstimate(i) = finalCIEstimatesReshape(i, CIMaxArgs(i));
end
%%
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

%export_fig SAFEWarehouseComparison.png -transparent -m2 -painters

save('MeshSAFEAverage.mat', 'recoveredImage', 'SAFEWorstEstimate', 'CIWorstEstimate', 'grid', '-append');

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

%export_fig SAFEWarehouseAgent.png -transparent -m2 -painters

save('MeshSAFEAverage.mat', 'finalSAFEAbsError', 'finalCIAbsError', 'varTheta', '-append');
