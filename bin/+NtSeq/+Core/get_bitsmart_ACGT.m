function [bitsmartACGT_uint8] = get_bitsmart_ACGT(ntSeq)
    % GET_BITSMART_ACGT - converts a nucleotide sequence specified as
    %  nucleotide integers or as sequence chars in accordance with Matlab's
    %  variation of IUPAC nucleotide codes into a custom encoding
    %
    % Inputs:
    %   ntSeq
    %     The nucleotie sequence specified as uint8 nucleotide integers or
    %     as sequence chars.
    %     The encoding provided must be in accordance with Matlab's 
    %     variation of IUPAC nucleotide codes 
    %
    % Outputs:
    %   bitsmartACGT_uint8
    %     The nucleotide sequence specified as uint8 integers in our custom
    %     bitsmart encoding which differs substantially from Matlab's
    %     nucleotide integer encoding. See section below for details,
    %     noting that dnaNtInt refers to matlab's uint8 nucleotide integers
    %
    % BITSMART ENCODING DETAILS:
    %
    %  The custom encoding represents whether a Adenosine, Cytidine,
    %  Guanine, or Thymidine would fit the description encoded, each with
    %  an individual bit, such that faster bitwise operators can be used
    %  when processing the data
    %  http://mathworks.com/help/bioinfo/ref/nt2int.html#brdesnv-1
    %
    % The custom encoding can be summarized as follows:
    %
    % bitsmartACGT_uint8      bitsmartACGT_logicalArr       dnaNtInt    dnaDef           dnaNtCode    dnaNtDesc
    %                         [0 0 * -  A C G T]
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % 32                      [0 0 1 0  0 0 0 0]            0           [0]               *           Unknown (any character not in table)
    % 8                       [0 0 0 0  1 0 0 0]            1           [1]               A           Adenosine
    % 4                       [0 0 0 0  0 1 0 0]            2           [2]               C           Cytidine
    % 2                       [0 0 0 0  0 0 1 0]            3           [3]               G           Guanine
    % 1                       [0 0 0 0  0 0 0 1]            4           [4]               T           Thymidine
    % 10                      [0 0 0 0  1 0 1 0]            5           [1|3]             R           Purine (A or G)
    % 5                       [0 0 0 0  0 1 0 1]            6           [2|4]             Y           Pyrimidine (T or C)
    % 3                       [0 0 0 0  0 0 1 1]            7           [3|4]             K           Keto (G or T)
    % 12                      [0 0 0 0  1 1 0 0]            8           [1|2]             M           Amino (A or C)
    % 6                       [0 0 0 0  0 1 1 0]            9           [2|3]             S           Strong interaction (3 H Bonds) (G  or C)
    % 9                       [0 0 0 0  1 0 0 1]            10          [1|4]             W           Weak interaction (2 H Bonds) (A or T)
    % 7                       [0 0 0 0  0 1 1 1]            11          [2|3|4]           B           Not A (C or G or T)
    % 11                      [0 0 0 0  1 0 1 1]            12          [1|3|4]           D           Not C (A or G or T)
    % 13                      [0 0 0 0  1 1 0 1]            13          [1|2|4]           H           Not G (A or C or T)
    % 14                      [0 0 0 0  1 1 1 0]            14          [1|2|3]           H           Not T (A or C or G)
    % 15                      [0 0 0 0  1 1 1 1]            15          [1|2|3|4]         N           Any nucleotide (A or C or G or T)
    % 16                      [0 0 0 1  0 0 0 0]            16          [1|2|3|4]*        -           Gap of indeterminate length
    %
    % Authors:
    %   Saair Quaderi
    
    if isempty(ntSeq)
        ntIntSeq = zeros(0,1);
    elseif ischar(ntSeq)
        validateattributes(ntSeq, {'char'}, {'row'});
        ntIntSeq = nt2int(ntSeq);
    elseif isa(ntSeq, 'uint8')
        validateattributes(ntSeq, {'uint8'}, {'row'});
        ntIntSeq = ntSeq;
    else
        error('Nucleotide sequence must be uint8 or char sequence');
    end
    
    bitsmartTranslationArr = uint8([32 8 4 2 1 10 5 3 12 6 9 7 11 13 14 15 16]);
    bitsmartACGT_uint8 = bitsmartTranslationArr(ntIntSeq + 1);
end