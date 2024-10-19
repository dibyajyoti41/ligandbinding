function [varargout] = data_hash( varargin )
    % DATA_HASH - Generates a hash of the matlab data provided.
    %
    %  This just wraps the third party functionality available in
    %  ThirdParty.DataHash.DataHash, so see its documentation for details
    import ThirdParty.DataHash.DataHash;
    [varargout{1:nargout}] = feval(@DataHash, varargin{:});
end