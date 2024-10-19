function [probBinding, competitors] = calculate_competitive_binding_probs(ntSeq, competitors, skipProbCalcCompetitors, isBitsmartEncodedTF)
    % CALCULATE_COMPETITIVE_BINDING_PROBS - Calculate Competitive binding
    %   probabilities for multiple types of ligands
    %
    % Calculates the equilibrium binding probability for each basepair
    % along a DNA strand for multiple types of competing ligands
    %
    % The transfer matrix method is from
    %  Teif, Nucl. Ac. Res. 35, e80 (2007)
    %   http://nar.oxfordjournals.org/content/35/11/e80.full
    %   http://nar.oxfordjournals.org/content/42/15/e118.full
    %
    % Inputs:
    %   ntSeq
    %     DNA sequence with standard character/matlab uint8 encoding
    %       (unless isBitsmartEncodedTF is set to true, in which case
    %        it should be encoded using NtSeq.Core.get_bitsmart_ACGT)
    %   competitors
    %     struct containing binding constants and bulk concentrations for
    %     each of the competing ligands
    %   skipProbCalcCompetitors
    %     list of competitors (named as fields in competitors struct) for
    %     which to skip calculating probabilities
    %   isBitsmartEncodedTF
    %     set to true if and only if the ntSeq is bitsmart encoded
    %
    % Outputs:
    %   probBinding
    %     probability of each competitor binding
    %   competitors
    %     same as input but with any newly computed binding constant values
    %     included
    %
    % Authors:
    %   Saair Quaderi (current version)
    %   Erik Lagerstedt (earlier version)
    %   Tobias Ambj�rnsson (earlier version)


    % TODO Note: The ordering of states represented in the transfer
    %  matrices and associated vectors appear to be a bit more
    %  complicated than necessary (since they are that way in the paper)
    %  They can probably be simplified and the matrix multiplications
    %  by which these are used to calculate probabilities can probably
    %  also be expanded explicitly to multiplications and additions of
    %  scalar non-zero values in a manner that runs significantly faster
    %  than the current approach with non-sparse matrices and mtimes

    if nargin < 3
        skipProbCalcCompetitors = {};
    end
    
    if nargin < 4
        isBitsmartEncodedTF = false;
    end
    
    if isBitsmartEncodedTF
        validateattributes(ntSeq, {'uint8'}, {'row'});
        ntBitsmartSeq = ntSeq;
    else
        if isempty(ntSeq)
            ntIntSeq = zeros(0,1);
        elseif ischar(ntSeq)
            validateattributes(ntSeq, {'char'}, {'row'});
            ntIntSeq = nt2int(ntSeq);
        elseif isa(ntSeq, 'uint8')
            validateattributes(ntSeq, {'uint8'}, {'row'});
            ntIntSeq = ntSeq;
        end
        clear ntSeq;
        
        import NtSeq.Core.get_bitsmart_ACGT;
        ntBitsmartSeq = get_bitsmart_ACGT(ntIntSeq);
        clear ntIntSeq;
    end
    
    if any(ntBitsmartSeq > (2^5 - 1)) || any(ntBitsmartSeq == 0)
        error('Sequences with gaps or nonstandard nucleotide labels are not supported for competitive binding theory curve generation');
    end

    subsystemNames = fieldnames(competitors);
    if any(strcmp(subsystemNames, 'Main'))
        error('Main is reserved');
    end
    if any(strcmp(subsystemNames, 'Free'))
        error('Free basepair subsystem is reserved');
    end
    subsystemNames(strcmp(subsystemNames, 'cooperativity')) = [];
    subsystemNames = [{'Free'}; subsystemNames(:)];
    numSubsystems = length(subsystemNames);
    
    numBpStates = 0;
    import CBT.Core.gen_binding_constants_mat;
    for subsystemNum = 1:numSubsystems
        subsystemName = subsystemNames{subsystemNum};
        if strcmp(subsystemName, 'Free')
            bpsOccupied = 1;
        else
            if not(isfield(competitors.(subsystemName), 'bulkConc'))
                    error(['Bulk concentration was not specified for ''', subsystemName , '''']);
            end
            if not(isfield(competitors.(subsystemName), 'bindingConstantsACGT'))
                if not(isfield(competitors.(subsystemName), 'bindingConstantRules'))
                    error(['Neither explicit binding constants nor rules for the binding constants were specified for ''', subsystemName , '''']);
                else
                    competitors.(subsystemName).bindingConstantsACGT = gen_binding_constants_mat(competitors.(subsystemName).bindingConstantRules);
                end
            end
            bindingConstantsACGT = competitors.(subsystemName).bindingConstantsACGT;
            bpsOccupied = ndims(bindingConstantsACGT);
            if bpsOccupied < 2
                error(['Binding for ''', subsystemName ,''' does not occupy two or more basepairs']);
            end
            if any(size(bindingConstantsACGT) ~= 4)
                error(['Binding constants for ''', subsystemName ,''' must be specified for four types of nucleotides for each occupied basepair']);
            end
        end
        stateIdx.(subsystemName).Start = numBpStates + 1;
        stateIdx.(subsystemName).End = stateIdx.(subsystemName).Start + (bpsOccupied - 1);
        numBpStates = stateIdx.(subsystemName).End;
    end

    transferMatrix.Main = zeros(numBpStates, numBpStates);
    constraints.allowedState.beforeFirstBp = zeros(1, numBpStates);
    constraints.allowedState.afterLastBp = zeros(numBpStates, 1);

    seqBindingConstants = struct;

    import CBT.Core.get_sequence_binding_constants;
    for subsystemNum = 1:numSubsystems
        subsystemName = subsystemNames{subsystemNum};
        switch subsystemName
            case 'Free'
                constraints.allowedState.afterLastBp(stateIdx.(subsystemName).End) = 1;
            otherwise
                bindingIdx.(subsystemName) = sub2ind([numBpStates, numBpStates], stateIdx.(subsystemName).End, stateIdx.(subsystemName).End - 1);
                seqBindingConstants.(subsystemName) = competitors.(subsystemName).bulkConc * get_sequence_binding_constants(...
                    competitors.(subsystemName).bindingConstantsACGT,...
                    ntBitsmartSeq);
        end

        constraints.allowedState.beforeFirstBp(stateIdx.(subsystemName).End) = 1;
        constraints.allowedState.beforeFirstBp(stateIdx.(subsystemName).End) = 1;
        constraints.allowedState.beforeFirstBp(stateIdx.(subsystemName).End) = 1;
        
        transferMatrix.(subsystemName) = zeros(numBpStates, numBpStates);
        for subsystemNumB = 1:numSubsystems
            subsystemNameB = subsystemNames{subsystemNumB};
            if isfield(competitors, 'cooperativity') && isfield(competitors.cooperativity, subsystemName) && isfield(competitors.cooperativity.(subsystemName), subsystemNameB)
                cooperativityAB = competitors.cooperativity.(subsystemName).(subsystemNameB);
            else
                cooperativityAB = 1;
            end
            transferMatrix.(subsystemName)(stateIdx.(subsystemName).Start, stateIdx.(subsystemNameB).End) = cooperativityAB;
        end
        for idx = (stateIdx.(subsystemName).Start + 1):(stateIdx.(subsystemName).End - 1)
            transferMatrix.(subsystemName)(idx, idx - 1) = 1;
        end
        transferMatrix.Main = transferMatrix.Main + transferMatrix.(subsystemName);
    end

    numBasepairs = length(ntBitsmartSeq);
    constraints.allowedState.beforeFirstBp = constraints.allowedState.beforeFirstBp / norm(constraints.allowedState.beforeFirstBp);

    import CBT.Core.right_left_mtx;
    direction = {'Rightwards'; 'Leftwards'};
    numStartSides = length(direction);
    for startSideNum = 1:numStartSides
        startSideName = direction{startSideNum};
        mtimesReverse = false;
        % tmpConstraintMat = [];
        % tmpNormConstants = [];
        switch startSideName
            case 'Rightwards'
                % tmpConstraintMat = zeros(1, numBpStates, numBasepairs);
                % tmpNormConstants = zeros(1, numBasepairs);
                % tmpVectorPrev = constraints.allowedState.beforeFirstBp;
                tmpVectorPrev2 = constraints.allowedState.beforeFirstBp;
                seqBpIdxStart = 1;
                seqBpIdxInterval = 1;
                seqBpIdxEnd = numBasepairs;
            case 'Leftwards'
                mtimesReverse = true;
                % tmpConstraintMat = zeros(numBpStates, 1, numBasepairs);
                % tmpNormConstants = zeros(numBasepairs, 1);
                % tmpVectorPrev = constraints.allowedState.afterLastBp;
                tmpVectorPrev2 = constraints.allowedState.afterLastBp;
                seqBpIdxStart = numBasepairs;
                seqBpIdxInterval = -1;
                seqBpIdxEnd = 1;
        end
        tmpTransferMatrix = transferMatrix.Main;
        tmpBindingConstants = zeros((numSubsystems - 1), numBasepairs);
        tmpBindingIdxs = zeros(numSubsystems - 1, 1);
        numInserts = 0;
        for subsystemNum = 1:numSubsystems
            subsystemName = subsystemNames{subsystemNum};
            switch subsystemName
                case 'Free'
                otherwise
                    numInserts = numInserts + 1;
                    tmpBindingIdxs(numInserts) = bindingIdx.(subsystemName);
                    tmpBindingConstants(numInserts, :) = seqBindingConstants.(subsystemName);
            end
        end

        [normConstants.(startSideName), probConstraintMat.(startSideName)] = ...
            right_left_mtx(seqBpIdxStart, seqBpIdxInterval, seqBpIdxEnd, ...
                                    tmpBindingConstants, tmpBindingIdxs, ...
                                    tmpVectorPrev2, tmpTransferMatrix, ...
                                    mtimesReverse);
    end

    tmpNumerators = [1, cumprod(normConstants.Rightwards(1:(numBasepairs - 1)) ./ normConstants.Leftwards(1:(numBasepairs - 1)))];
    tmpDenominators = normConstants.Leftwards .* mtimes(constraints.allowedState.beforeFirstBp, probConstraintMat.Leftwards(:, 1, 1));
    probNormalizationVector = tmpNumerators ./ tmpDenominators;

    import CBT.Core.binding_probs_mtx;
    skipProbCalcCompetitors = [skipProbCalcCompetitors(:); {'Free'}];
    for subsystemNum = 1:numSubsystems
        subsystemName = subsystemNames{subsystemNum};
        if any(strcmp(skipProbCalcCompetitors, subsystemName))
            continue;
        end
        
        tmpTransferMatrix = transferMatrix.(subsystemName);

        tmpBindingProb = binding_probs_mtx( ...
                            numBasepairs, ...
                            seqBindingConstants.(subsystemName), ...
                            bindingIdx.(subsystemName), ...
                            tmpTransferMatrix, ...
                            constraints.allowedState.beforeFirstBp, ...
                            constraints.allowedState.afterLastBp, ...
                            probConstraintMat.Leftwards, ...
                            probConstraintMat.Rightwards );
        probBinding.(subsystemName) = tmpBindingProb .* probNormalizationVector; 
    end
end
