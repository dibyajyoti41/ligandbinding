function wrapperFn = arg_swap(fn, argIndices)
    % ARG_SWAP - returns a function wrapper that manipulates arguments
    %   This is a higher order function. The function it returns returns
    %    the result of fn for arguments rearranged in accordance to
    %    argIndices
    %   The rearrangement may change the order of arguments, remove
    %    arguments, and/or duplicate arguments. If fn takes k
    %    arguments, wrapperFn with argIndices with value (1:k) would behave
    %    just like fn.
    %
    %   (Special case: if argIndices contains an index number greater than
    %    the number of inputs provided to wrapperFn, the tail of argIndices
    %     is ignored to exclude the index number and all subsequent values)
    %
    % Inputs:
    %   fn
    %     function handle for function being wrapped
    %   argIndices (optional)
    %     a row vector containing ordere indices of arguments to be passed
    %       on to fn from wrapperFn
    %     defaults to [] (e.g. no arguments are passed on)
    %
    %  Outputs:
    %    wrapperFn
    %      function that returns the result of fn for arguments rearranged
    %      in accordance to argIndices
    %
    %  Examples:
    %    For instance, if you have a function:
    %      fn1 = @(a, b, c) c*b - a;
    %    the output produced by arg_swap(fn1, [2, 1, 3])
    %      would be equivalent to @(a, b, c) c*a - b;
    %    the output produced by arg_swap(fn1, [3, 1, 1])
    %      would be equivalent to @(a, ~, c) a*a - c;
    %
    %   if you have a function that takes no inputs and want to call
    %     it as a callback function which will receive some inputs, 
    %     all of which should be ignored, you can just wrap it with arg_swap:
    %     fnNoInputs = @() disp('Example');
    %     fnIgnoreInputs = arg_swap(fnNoInputs);
    %
    % Authors:
    %   Saair Quaderi

    validateattributes(fn, {'function_handle'}, {'scalar'});
    if nargin < 2
        argIndices = [];
    elseif not(isempty(argIndices))
        validateattributes(argIndices, {'numeric'}, {'positive', 'integer', 'row'});
    end
    
    function varargout = wrapper(varargin)
        unprovidedArgIdx = find([argIndices, nargin + 1] > nargin, 1);
        argIndices = argIndices(1:(unprovidedArgIdx - 1));
        argsSwapped = varargin(argIndices);
        [varargout{1:nargout}] = feval(fn, argsSwapped{:});
    end
    wrapperFn = @wrapper;
end

