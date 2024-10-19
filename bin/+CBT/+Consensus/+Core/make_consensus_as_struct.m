function [consensusStruct, cache] = make_consensus_as_struct(barcodes, barcodeBitmasks, barcodeAliases, otherBarcodeData, clusterThresholdScoreNormalized, cache, rawBarcodes, rawBgs)

    % Validate/standardize barcodes cell array
    numBarcodes = size(barcodes, 1);
    if (numBarcodes < 2) 
         if numBarcodes == 1 % in case only one barcode, return it instead of an error, todo: make this nicer 
            consensusStruct = struct;
            consensusStruct.finalConsensusKey = '1';
            consensusStruct.rawBarcode = rawBarcodes{1};
            consensusStruct.rescaledBarcode = barcodes{1};
            consensusStruct.bitmask = barcodeBitmasks{1};
            consensusStruct.name = barcodeAliases{1};
            consensusStruct.rawBgs = rawBgs;
            consensusStruct.cache = cache;
            consensusStruct.timestamp = datetime('now');
            consensusStruct.formatVersion = '0.1.0';            
            return;
         else    
            error('Consensus requires multiple barcodes');
         end
    end
    if not(iscell(barcodes))
        barcodes = mat2cell(barcodes, ones(numBarcodes, 1), size(barcodes, 2));
    else
        validateattributes(barcodes, {'cell'}, {'column'}, 1);
    end

    % Validate/standardize bitmasks cell array
    numBitmasks = size(barcodeBitmasks, 1);
    if (numBitmasks ~= numBarcodes)
        error('There must be the same number of bitmasks as barcodes');
    end
    if not(iscell(barcodeBitmasks))
        barcodeBitmasks = mat2cell(barcodeBitmasks, ones(numBitmasks, 1), size(barcodeBitmasks, 2));
    else
        validateattributes(barcodeBitmasks, {'cell'}, {'column'}, 2);
    end

    % Validate individual barcodes & bitmasks (length & shape)
    commonLength = size(barcodes{1}, 2);
    if not(all(cellfun(@(rowVect) isrow(rowVect) & (length(rowVect) == commonLength), [barcodes; barcodeBitmasks])))
        error('All barcodes and bitmasks must be row vectors of equal length');
    end

    % Validate/standardize barcode aliases
    if (nargin < 3) || isempty(barcodeAliases)
        barcodeAliases = cell(numBarcodes, 1);
        for barcodeNum=1:numBarcodes
            barcodeAliases{barcodeNum} = '';
        end
    else
        % Validate individual aliases as strings
        %  (empty strings or char row vectors)
        validateattributes(barcodeAliases, {'cell'}, {'column', 'numel', numBarcodes}, 3);
        if not(iscellstr(barcodeAliases)) || not(all(cellfun(@(barcodeAlias) (isempty(barcodeAlias) || isrow(barcodeAlias)), barcodeAliases)))
            error('Barcode aliases must be a column cell vector of strings');
        end
    end

    % Validate/standardize other data
    if (nargin < 4) || isempty(otherBarcodeData)
        otherBarcodeData = arrayfun(@(x) {},  zeros(numBarcodes, 1), 'UniformOutput',false);
    end

    % Validate/standardize cluster thresholds
    defaultClusterThresholdScoreNormalized = 0.75;
    bestPossibleScore = sqrt(commonLength); % potential overestimate since this length implies that bitmasks contain all 1s
    if (nargin < 5)
        clusterThresholdScoreNormalized = defaultClusterThresholdScoreNormalized;
    else
        validateattributes(clusterThresholdScoreNormalized, {'numeric'}, {'scalar', 'nonnegative', '<=', 1}, 5);
    end
    clusterThresholdScore = clusterThresholdScoreNormalized*bestPossibleScore;

    % Validate/standardize cache
    if (nargin < 6)
        cache = containers.Map();
    else
        validateattributes(cache, {'containers.Map'}, {}, 6);
    end


    inputs = struct;

    inputs.barcodes = barcodes;
    inputs.barcodeBitmasks = barcodeBitmasks;
    inputs.barcodeAliases = barcodeAliases;
    inputs.otherBarcodeData = otherBarcodeData;
    inputs.clusterThresholdScore = clusterThresholdScore;
    inputs.rawBarcodes = rawBarcodes;
    inputs.rawBgs = rawBgs;

    barcodeStructsMap = containers.Map();
    consensusKeyPool = cell(numBarcodes, 1);
    for keyNum=1:numBarcodes
        barcodeKey = num2str(keyNum);
        consensusKeyPool{keyNum} = barcodeKey;
        barcodeStruct = struct;
        barcodeStruct.maxWeight = 1;
        barcodeStruct.indexWeights = barcodeBitmasks{keyNum};
        barcodeStruct.alias = barcodeAliases{keyNum};
        barcodeStruct.barcode = barcodes{keyNum};
        barcodeStruct.parents = {};
        barcodeStruct.bestScore = NaN;
        barcodeStruct.xcorrAtBest = NaN;
        barcodeStructsMap(barcodeKey) = barcodeStruct;
    end

    keyList = cell(2*numBarcodes - 1, 1);
    keyList(1:numBarcodes) = consensusKeyPool;
    for keyNum=(numBarcodes+1):(2*numBarcodes - 1)
        keyList{keyNum} = '';
    end
    consensusMergingTree = zeros((numBarcodes - 1), 3);
    relevantCacheKeys = cell(numBarcodes - 1, 1);
    import CBT.Consensus.CoreCacheBridge.get_best_pair_indices;
    import CBT.Consensus.Core.merge_consensus_pair;
    for subtreeIndex = 1:(numBarcodes - 1)
        [keyA, keyB, bestScore, xcorrAtBest, flipTFAtBest, circShiftAtBest, cache, subtreeCacheKeys] = get_best_pair_indices(consensusKeyPool, barcodeStructsMap, cache);
        relevantCacheKeys{subtreeIndex} = subtreeCacheKeys;...
        bestScoreNormalized = bestScore/bestPossibleScore;
        keyAB = ['[', keyA, ',', keyB ']'];
        keyList{numBarcodes + subtreeIndex} = keyAB;

        structA = barcodeStructsMap(keyA);
        structB = barcodeStructsMap(keyB);

        [barcodeAB, indexWeightsAB] = merge_consensus_pair(structA.barcode, structB.barcode, structA.indexWeights, structB.indexWeights, flipTFAtBest, circShiftAtBest);

        barcodeStructAB = struct;
        barcodeStructAB.maxWeight = structA.maxWeight + structB.maxWeight;
        barcodeStructAB.indexWeights = indexWeightsAB;
        barcodeStructAB.alias = '';
        barcodeStructAB.barcode = barcodeAB;
        barcodeStructAB.parents = { ...
            {keyA, false, 0}, ... % A wasn't be reoriented into sync with B (prior to merging A & B)...
            {keyB, flipTFAtBest, circShiftAtBest} ... % instead B was reoriented (flipped/cyclically-permutated as necessary) into sync with A (prior to merging A & B)
        };
        barcodeStructAB.bestScore = bestScoreNormalized;
        barcodeStructAB.xcorrAtBest = xcorrAtBest;
        barcodeStructsMap(keyAB) = barcodeStructAB;

        % setxor effectively replaces keyA and keyB with keyAB
        % since the former should be in the pool and the latter
        % shouldn't be
        consensusKeyPool = setxor(consensusKeyPool, {keyA, keyB, keyAB});

        idxKeyA = find(strcmp(keyList, keyA), 1);
        idxKeyB = find(strcmp(keyList, keyB), 1);
        consensusMergingTree(subtreeIndex, :) = [idxKeyA, idxKeyB, (1 - bestScoreNormalized)];
    end
    relevantCacheKeys = vertcat(relevantCacheKeys{:}); % unify into single string cell column
    relevantCacheKeys = unique(relevantCacheKeys); % remove all the duplicates

    subcache = [containers.Map(); cache]; % create a full copy of the cache
    remove(subcache, setdiff(keys(cache), relevantCacheKeys)); % remove anything from the copy that isn't related to this consensus


    getKeyComponentNums = @(k) str2num(strrep(strrep(strrep(k, ']', ''), '[', ''), ',', '; ')); %#ok<ST2NM>
    formatSortedDash = @(componentNums) strjoin(arrayfun(@(x) num2str(x), num2str(sort(componentNums)), 'UniformOutput', false), '-');
    keyListSimplified = cellfun(@(k) formatSortedDash(getKeyComponentNums(k)), keyList, 'UniformOutput', false);

    clusterAssignmentsMatrix = cluster(consensusMergingTree, 'cutoff', 1 - clusterThresholdScoreNormalized, 'criterion', 'distance');
    numClusters = max(clusterAssignmentsMatrix);
    clusterKeys = cell(numClusters, 1);
    clusterResultStructs = cell(numClusters, 1);
    clusterSizes = zeros(numClusters, 1);
    import CBT.Consensus.Helper.extract_consensus_aligned_barcodes;
    for clusterNum=1:numClusters
        componentNums = find(clusterAssignmentsMatrix == clusterNum);
        clusterKey = keyList{find(strcmp(keyListSimplified, formatSortedDash(componentNums)), 1)};
        clusterSizes(clusterNum) = length(componentNums);
        clusterKeys{clusterNum} = clusterKey;
        clusterResultStruct = struct;
        [ ...
            clusterResultStruct.barcodeKeys, ...
            clusterResultStruct.alignedBarcodes, ...
            clusterResultStruct.alignedBarcodeBitmasks, ...
            clusterResultStruct.barcodes, ...
            clusterResultStruct.barcodeBitmasks, ...
            clusterResultStruct.flipTFs, ...
            clusterResultStruct.circShifts ...
            ] = extract_consensus_aligned_barcodes(clusterKey, barcodeStructsMap);
        clusterResultStructs{clusterNum} = clusterResultStruct;
    end

    consensusStruct = struct;
    consensusStruct.inputs = inputs;
    consensusStruct.keyList = keyList;
    consensusStruct.keyListSimplified = keyListSimplified;
    consensusStruct.finalConsensusKey = consensusKeyPool{1};
    consensusStruct.barcodeStructsMap = barcodeStructsMap;
    consensusStruct.consensusMergingTree = consensusMergingTree;
    consensusStruct.clusterKeys = clusterKeys;
    consensusStruct.clusterAssignmentsMatrix = clusterAssignmentsMatrix;
    consensusStruct.clusterResultStructs = clusterResultStructs;
    consensusStruct.cache = subcache;
    consensusStruct.timestamp = datetime('now');
    consensusStruct.formatVersion = '0.1.0';
end