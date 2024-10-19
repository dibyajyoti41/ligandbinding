function [varargout] = feval_with_structified_args(fn, argStruct)
    % FEVAL_WITH_STRUCTIFIED_ARGS - evaluates a function with the 
    %  fieldnames and their values in the struct as interspersed arguments
    %
    % Inputs:
    %   fn
    %     the function to run
    %   argStruct
    %     the struct containing the arguments for the function as fields
    %     and values
    %
    % Outputs: (variable)
    %   the outputs of the evaluated function with the provided inputs
    %   using the same nargout as this function
    %
    % e.g.
    %   s = struct('A', 1, 'B', 2);
    %   feval_with_structified_args(fn, s);
    %   % is equivalent to fn('A', 1, 'B', 2)l
    %
    % Authors:
    %   Saair Quaderi
    
    [varargout{1:nargout}] = feval(@(s, fn2) feval(@(x) fn2(x{:}),...
        permute([...
            fieldnames(s),...
            cellfun(...
                @(fname) s.(fname),...
                fieldnames(s),...
                'UniformOutput', false)],...
            [2 1])),...
        argStruct,...
        fn);
end