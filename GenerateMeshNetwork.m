% Yuan Chen

clear all;
close all;

set(0,'defaulttextinterpreter','none')

%% Set up Agents

cyberMeshSize = 25; 

N = cyberMeshSize^2; %Number of Agents
sectorMeshSize = 230;
GridDimension = sectorMeshSize;
unitDimension = GridDimension/(cyberMeshSize-1);
sensingRange = round(1.5*unitDimension);
maxVisible = min(sectorMeshSize^2, (2*sensingRange + 1)^2);

interestRange = round(2*sensingRange);
maxInterest = min(sectorMeshSize^2, (2*interestRange + 1)^2);

m = sectorMeshSize^2;

if (N < 4) %At least 4 agents
    return;
end

HArray = cell(N, 1);
LocationArray = cell(N, 1); %x, y positions of agents
PArray = cell(N, 1);
BigPArray = cell(N, 1);

sectorLocationArray = cell(N, 1); % Assign agent to a particular "sector"

sectorIndexArray = cell(m, 1);
sectorCount = zeros(m, 1);

sectorMembers = cell(m, 1);

for i = 1:m
    sectorMembers{i} = cell(0, 1);
end



for i = 1:N
    %tempH = zeros(1, m); %Initially, no agent has any measurement
    rowIndex = floor((i-1)/cyberMeshSize)+1;
    columnIndex = mod((i-1), cyberMeshSize)+1;
    
    currentLocation = [rowIndex*unitDimension; columnIndex*unitDimension] - 0.99*[unitDimension; unitDimension];
    
    if (rowIndex == cyberMeshSize)
        currentLocation(1) = 0.99*GridDimension;
    end
    if (columnIndex == cyberMeshSize)
        currentLocation(2) = 0.99*GridDimension;
    end
    
    %LocationArray{i} = GridDimension*rand(2, 1);
    LocationArray{i} = currentLocation;
    sectorLocation = ceil(LocationArray{i});
    sectorLocationArray{i} = sectorLocation;
    
    visibleLocationArray = cell(maxVisible, 1);
    interestedLocationArray = cell(maxInterest, 1);
    %visibleLocationArray{1} = sectorLocation;
    
    totalVisibleLocations = 0;
    totalInterestedLocations = 0;
    
%     for offset = -sensingRange:sensingRange
%         for offset2 = -sensingRange:sensingRange
%             if (((sectorLocation(1) + offset) > 0) && ((sectorLocation(1) + offset) <= sectorMeshSize))
%                 if (((sectorLocation(2) + offset2) > 0) && ((sectorLocation(2) + offset2) <= sectorMeshSize))
%                     totalVisibleLocations = totalVisibleLocations + 1;
%                     tempVisibleLocation = [sectorLocation(1) + offset; sectorLocation(2) + offset2];
% 
%                     visibleLocationArray{totalVisibleLocations} = tempVisibleLocation;
%                 end
%             end
%         end
%     end

    for offset = -interestRange:interestRange
        for offset2 = -interestRange:interestRange
            if (((sectorLocation(1) + offset) > 0) && ((sectorLocation(1) + offset) <= sectorMeshSize))
                if (((sectorLocation(2) + offset2) > 0) && ((sectorLocation(2) + offset2) <= sectorMeshSize))
                    totalInterestedLocations = totalInterestedLocations + 1;
                    tempInterestedLocation = [sectorLocation(1) + offset; sectorLocation(2) + offset2];
                    
                    interestedLocationArray{totalInterestedLocations} = tempInterestedLocation;
                    
                    if ((abs(offset) <= sensingRange) && (abs(offset2) <= sensingRange))
                        totalVisibleLocations = totalVisibleLocations + 1;
                        visibleLocationArray{totalVisibleLocations} = tempInterestedLocation;
                    end
                end
            end
        end
        
    end
    
    tempH = zeros(totalVisibleLocations, m);
    
    for currentVisLoc = 1:totalVisibleLocations
        currentLoc = visibleLocationArray{currentVisLoc};
        
        currentSectorIndex = (currentLoc(1)-1)*sectorMeshSize + currentLoc(2);
        tempH(currentVisLoc, currentSectorIndex) = 1;
    end
    
    sectorIndex = (sectorLocation(1)-1)*sectorMeshSize + sectorLocation(2);
    
    sectorIndexArray{i} = sectorIndex;
    %tempH(sectorIndex) = 1;
    
    
    sectorCount(sectorIndex) = sectorCount(sectorIndex) + 1;
    
    sectorMembers{sectorIndex}{sectorCount(sectorIndex)} = i;
    
    HArray{i} = sparse(tempH);
    
    clear tempH;
    
    %tempP = spalloc(m, 1, totalInterestedLocations);
    interestedIndices = [];
    for currentIntLoc = 1:totalInterestedLocations
        currentLoc = interestedLocationArray{currentIntLoc};
        
        currentSectorIndex = (currentLoc(1)-1)*sectorMeshSize + currentLoc(2);
        %tempP(currentSectorIndex) = 1;
        interestedIndices = cat(1, interestedIndices, currentSectorIndex);
    end
    tempP = sparse(interestedIndices, ones(totalInterestedLocations, 1), ones(totalInterestedLocations, 1), m, 1);
    PArray{i} = spdiags(tempP, 0, m, m);
    BigPArray{i} = tempP;
    
    clear tempP;
    
    if (mod(i, 50) == 0)
        fprintf('--- Agent %d Interests Assigned ---\n', i);
    end
end
fprintf('=== All Interests Assigned ===\n');

BigP = spdiags(cat(1, BigPArray{:}), 0, N*m, N*m);

% Compute Distance Matrix Between Agents
DistanceMatrix = zeros(N);

for i = 1:N
    for j = i:N
        distance = norm(LocationArray{i} - LocationArray{j});
        DistanceMatrix(i, j) = distance;
        DistanceMatrix(j, i) = distance;
    end
end
edgeThreshold = 1.5*unitDimension;

l2 = 0;

% while (l2 < 1e-9)
%     Adjacency = sparse((DistanceMatrix < edgeThreshold*ones(N)) - eye(N));
%     Diagonal = diag(sum(Adjacency, 2));
%     L = sparse(Diagonal - Adjacency);
% 
%     evalues = sort(eig(L));
%     l2 = evalues(2);
%     edgeThreshold = 1.05*edgeThreshold;
% end

Adjacency = sparse((DistanceMatrix < edgeThreshold*ones(N)) - eye(N));
Diagonal = diag(sum(Adjacency, 2));
L = sparse(Diagonal - Adjacency);

numSelected = 160;
numRemaining = N-numSelected;
remainingIndices = zeros(numRemaining, 1);
selectedAgents = zeros(numSelected, 1);

strongAdversarySector = 5;

strongAdversaries = cell2mat(sectorMembers{strongAdversarySector});
numAdversaries = sectorCount(strongAdversarySector);

subL2 = 0;
totalSelections = 0;
selectionLimit = 2;
allIndices = 1:N;

%% Generate Censored Laplacian
%numNonzeroBlocks = nnz(L);
%censoredL = spalloc(N*m, N*m, numNonzeroBlocks*m);
rowIndices = [];
columnIndices = [];
values = [];
for i = 1:N
    currentDiag = sparse(m, m);
    for j = 1:N
        if (i == j)
            continue;
        end
        if (L(i, j) ~= 0)
            currentBlock = PArray{i}*PArray{j};
            currentDiag = currentDiag + currentBlock;
            [currentRIndices, currentCIndices, currentValues] = find(currentBlock);
            rowIndices = cat(1, rowIndices, currentRIndices + (i-1)*m);
            columnIndices = cat(1, columnIndices, currentCIndices + (j-1)*m);
            values = cat(1, values, -currentValues);
            %censoredL((i-1)*m + 1:i*m, (j-1)*m + 1:j*m) = -PArray{i}*PArray{j};
        end
    end
    [currentRIndices, currentCIndices, currentValues] = find(currentDiag);
    rowIndices = cat(1, rowIndices, currentRIndices + (i-1)*m);
    columnIndices = cat(1, columnIndices, currentCIndices + (i-1)*m);
    values = cat(1, values, currentValues);
    %censoredL((i-1)*m + 1:i*m, (i-1)*m + 1:i*m) = currentDiag;
    
    if (mod(i, 50) == 0)
        fprintf('--- Agent %d Censored Laplacian Computed ---\n', i);
    end
end

censoredL = sparse(rowIndices, columnIndices, values, N*m, N*m);
fprintf('=== Censored Laplacian Computed ===\n');

%% Choose adversaries
numSelected = 70;
adversaries = randsample(1:N, numSelected);

adversaryIndicator = zeros(N, 1);
adversaryIndicator(adversaries) = 1;

%% Plot Graph of Agents
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
end


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

save('PlotNetwork.mat', 'axisFactor', 'N', 'LocationArray', 'adversaryIndicator', 'sectorMeshSize', 'GridDimension', 'Adjacency')

%% Compute parameters
DH = sparse(blkdiag(HArray{:}));

numVisible = zeros(N, 1);
for i = 1:N
    temp = size(HArray{i});
    numVisible(i) = temp(1);
end

calH = cat(1, HArray{:});
observeGrammian = calH'*calH;

normalHArray = HArray;

for i = adversaries
    currentHSize = size(HArray{i});
    newH = sparse(currentHSize(1), currentHSize(2));
    normalHArray{i} = newH;
end

normalCalH = cat(1, normalHArray{:});
normalObserveGrammian = normalCalH'*normalCalH;
%% Generate Incidence Matrix
%[destination, source, ~] = find(Adjacency); %Note, adjacency matrix of undirected graph is symmetric
% 
% incidence = sparse(zeros(length(source), N));
% flagCollector = sparse(zeros(N, length(source)));
% 
% for i = 1:length(source)
%     flagCollector(source(i), i) = 1;
%     incidence(i, source(i)) = 1;
%     incidence(i, destination(i)) = -1;
% end
% 
% messageDifferenceMatrix = sparse(kron(incidence, speye(m)));
% normSumMatrix = sparse(kron(speye(length(source)), ones(1, m)));
% fprintf('--- Incidence Computed ---\n');

%export_fig WarehouseNetwork.png -transparent -m2 -painters

%strongAdversarySector = 5;

%strongAdversaries = cell2mat(sectorMembers(strongAdversaryIndex));

save('WarehouseNetwork.mat', 'L', 'N', 'm', 'DH', 'HArray', 'normalHArray', 'PArray', 'censoredL', 'BigP', 'numVisible', 'sectorMembers', 'sectorCount', 'adversaries');
%save('EXNetwork.mat', 'Adjacency', 'L', 'l2', 'N', 'numSelected', 'remainingIndices', 'selectedAgents');
