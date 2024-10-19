function [probBinding] = cb_netropsin_vs_yoyo1_plasmid(ntSeq, concNetropsin_molar, concYOYO1_molar, untrustedMargin, onlyYoyo1Prob, roundingPrecision)
    % CB_NETROPSIN_VS_YOYO1 = Calculate Competitive binding
    %   probabilities for Netropsin vs YOYO1 on a plasmid
    %
    % Calculates the equilibrium binding probability for each basepair
    % along a DNA strand for the two types of competing ligands
    %
    % The transfer matrix method is from
    %  Teif, Nucl. Ac. Res. 35, e80 (2007)
    %   http://nar.oxfordjournals.org/content/35/11/e80.full
    %   http://nar.oxfordjournals.org/content/42/15/e118.full
    %
    % Inputs:
    %   ntSeq
    %     DNA sequence with standard character/matlab uint8 encoding
    %   concNetropsin_molar
    %     Netropsin concentration (units molar)
    %   concYOYO1_molar
    %     YOYO1 concentration (units molar)
    %   untrustedMargin (optional, defaults to 1000)
    %     since calculate_competitive_binding_probs gives questionable
    %     results influenced by the choice of the initial and final state,
    %     at the start and end of the sequence, this parameter specifies
    %     how much temporary wrap-around padding (with padding values 
    %     taken from other end of the provided sequence as if it was
    %     cyclical) must be prepended to the ends of the sequence prior to
    %     computations such that the probabilities are at a reasonable
    %     equilibrium in the states surrounding the original sequence area
    %     which is extracted after computations are complete
    %   onlyYoyo1Prob (optional, default false)
    %     only outputs the probability of YOYO1 bind (array not struct of
    %     arrays)
    %   roundingPrecision (optional, default 8)
    %     rounds answers to roundingPrecision digits to the right of the
    %     decimal point (the results are not that accurate anyway so by
    %     rounding to some fixed precision it should make it easier to check if
    %     results are essentially the same when negligible rounding errors
    %     may be the result of negligible asymmetries in the
    %     implementation)
    %
    % Outputs:
    %   probBinding
    %      probability of binding for YOYO1 and Netropsin along DNA
    %      as a struct if onlyYoyo1Prob is not true, or as an array with
    %      only the probability of YOYO1 binding if onlyYoyo1Prob is true
    %
    % Authors:
    %   Saair Quaderi
    
    if (nargin < 4) || isempty(untrustedMargin)
        untrustedMargin = 1000;
    end

    if (nargin < 5) || isempty(onlyYoyo1Prob)
        onlyYoyo1Prob = false;
    end
    
    if (nargin < 6) || isempty(roundingPrecision)
        roundingPrecision = 8;
    end
    
    import CBT.get_binding_constant_rules;
    competitors.Netropsin.bindingConstantRules = get_binding_constant_rules('Netropsin');
    competitors.Netropsin.bulkConc = concNetropsin_molar;

    competitors.Yoyo1.bindingConstantRules = get_binding_constant_rules('Yoyo1');
    competitors.Yoyo1.bulkConc = concYOYO1_molar;

    if onlyYoyo1Prob
        skipProbCalcCompetitors = {'Netropsin'};
    else
        skipProbCalcCompetitors = {};
    end
    seqLen = length(ntSeq);
    if untrustedMargin > 0
        repeatedTimes = ceil(untrustedMargin/seqLen);
        paddedNtSeq = repmat(ntSeq(:), [repeatedTimes, 1]);
        paddedNtSeq = [paddedNtSeq((end - untrustedMargin) + (1:untrustedMargin)); paddedNtSeq; paddedNtSeq(1:untrustedMargin)]';
    else
        untrustedMargin = 0;
        paddedNtSeq = ntSeq;
    end
    unpaddedIdxs = untrustedMargin + (1:seqLen);
    

    if isempty(paddedNtSeq)
        paddedNtIntSeq = zeros(0,1);
    elseif ischar(paddedNtSeq)
        validateattributes(ntSeq, {'char'}, {'row'});
        paddedNtIntSeq = nt2int(paddedNtSeq);
    elseif isa(paddedNtSeq, 'uint8')
        validateattributes(paddedNtSeq, {'uint8'}, {'row'});
        paddedNtIntSeq = paddedNtSeq;
    end

    import NtSeq.Core.get_bitsmart_ACGT;
    paddedNtBitsmartSeq = get_bitsmart_ACGT(paddedNtIntSeq);
    isBitsmartEncodedTF = true;
    
    import CBT.Core.calculate_competitive_binding_probs;
    probBinding = calculate_competitive_binding_probs(paddedNtBitsmartSeq, competitors, skipProbCalcCompetitors, isBitsmartEncodedTF);
    
    % ntBitsmartSeq = paddedNtBitsmartSeq(unpaddedIdxs);
    structFieldnames = fieldnames(probBinding);
    numStructFieldnames = length(structFieldnames);
    for structFieldnameNum = 1:numStructFieldnames
        structFieldname = structFieldnames{structFieldnameNum};
        probBinding.(structFieldname) = round(probBinding.(structFieldname)(unpaddedIdxs), roundingPrecision);
    end
    if onlyYoyo1Prob
        probBinding = probBinding.Yoyo1;
    end
end