function [bindingConstantsMat] = gen_binding_constants_mat(bindingConstantRules)
    % GEN_BINDING_CONSTANTS_MAT - Expands out the binding constant rules
    %  specified to create a 4x4x...x4 matrix with k dimensions where k is
    %  the number of nucleotides in a bound sequence and the four
    %  indices (1, 2, 3, and 4) in the jth dimension are associated with
    %  an A, C, G, and T respectively in the jth character in the sequence
    %  to be bound to
    %
    % Inputs:
    %   bindingConstantRules
    %     Nx2 cell array with N > 1 rules
    %     - where first column represents a specification of nucleotide
    %     sequence characters with length k
    %     - where second column represents the binding constant for the
    %     specified sequence
    %     - where higher-indexed rows take precedence over lower-indexed
    %     rows
    %
    % Outputs:
    %   bindingConstantsMat
    %     4x4x...x4 matrix with k dimensions where k is the number of 
    %     nucleotides in a bound sequence and the four indices
    %     (1, 2, 3, and 4) in the jth dimension are associated with
    %     an A, C, G, and T respectively in the jth character in the
    %     sequence to be bound to
    %
    % Authors:
    %   Saair Quaderi
    
    import NtSeq.Core.get_bitsmart_ACGT;
    import NtSeq.Core.uint8_in_binary;
    
    numRules = size(bindingConstantRules, 1);
    if numRules < 1
        error('At least one rule must be specified');
    end
    if not(iscell(bindingConstantRules)) || not(isequal(size(bindingConstantRules), [numRules, 2]))
        error('Rules must be specified in an Nx2 cell array where the first column contains IUPAC nucleotide codes and the second column contains the associated numeric binding constants.');
    end
    seqSpecs = bindingConstantRules(:, 1);
    seqSpecLens = cellfun(@numel, seqSpecs);
    seqSpecLen = seqSpecLens(1);
    if any(seqSpecLens ~= seqSpecLen)
        error('Varying sequence specification lengths are not supported');
    end
    bindingConstantsMatSize = repmat(4, [1, seqSpecLen]);
    bindingConstantsMat = NaN(bindingConstantsMatSize);
    [bitsmartN] = get_bitsmart_ACGT('N');
    isValid = false;
    try
        isValid = all(arrayfun(@(bitsmartVal) bitor(bitsmartVal, bitsmartN) == bitsmartN, get_bitsmart_ACGT([seqSpecs{:}])));
    catch
    end
    if not(isValid)
        error('Only ''non-gap'' IUPAC nucleotide codes are supported.')
    end
    
    isValid = false;
    try
        bindingConstantVals = [bindingConstantRules{:,2}]';
        isValid = isnumeric(bindingConstantVals) && (numel(bindingConstantVals) == numRules);
    catch
    end
    if not(isValid)
        error('The second rules column must contain a numeric binding constant in each row');
    end
    for ruleNum=1:numRules
        mat_logical = uint8_in_binary(get_bitsmart_ACGT(seqSpecs{ruleNum}), 4);
        idxs = mat2cell(mat_logical, ones([1, size(mat_logical, 1)]), 4);
        bindingConstantsMat(idxs{:}) = bindingConstantVals(ruleNum);
    end
end