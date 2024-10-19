function [ mat_logical ] = uint8_in_binary( vect_uint8, k)
    % UINT8_IN_BINARY - Converts a vector of N uint8 values into a matrix
    %   of Nxk logical containing the bit values of the binary
    %   representations of each of the uint8 values
    %
    %   If k is specified, only th last k bits of the uint8s are provided.
    %   If k is not specified, k defaults to 8
    %
    % Inputs:
    %   vect_uint8
    %     a vector of N uint8 values
    %   k (optional, defaults to 8)
    %      the number of bits from the right end of the binary
    %      representation to express in the Nxk output matrix
    %
    % Outputs:
    %   mat_logical
    %     an Nxk matrix of logical values where each jth row contains the
    %     last k bits of the binary representation of the jth element in
    %     the uint8 input vector
    %
    % Authors:
    %   Saair Quaderi
    
    if nargin < 2
        k = 8;
    end
    mat_logical = logical(rem(floor(double(vect_uint8(:))*pow2(1 - k:0)),2));
end

