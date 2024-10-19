function [hUIElems] = generate_ui_elems(fn_create_control, overrideArgStructs, defaultArgStruct)
    if (nargin < 1) || isempty(fn_create_control)
        fn_create_control = @uicontrol;
    end
    if (nargin < 2) || isempty(overrideArgStructs)
        overrideArgStructs = cell(0, 1);
    end
    if (nargin < 3) || isempty(defaultArgStruct)
        defaultArgStruct = struct();
    end

    function [childObj] = generate_control(overrideArgStruct)
        import Fancy.Utils.feval_with_structified_args;
        import Fancy.Utils.merge_structs;
        
        argStruct = merge_structs(defaultArgStruct, overrideArgStruct);
        childObj = feval_with_structified_args(fn_create_control, argStruct);
    end
    hUIElems = cellfun(@(overrideArgStruct) generate_control(overrideArgStruct), overrideArgStructs, 'UniformOutput', false);
    hUIElems = reshape([hUIElems{:}], size(hUIElems));
end