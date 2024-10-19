function [] = reposition_nrm(hGOs, cellVectPosArrNrm)
    import Fancy.UI.FancyPositioning.set_position;
    isValidNonPlaceholder = arrayfun(@(gObj) isvalid(gObj) & not(isa(gObj, 'matlab.graphics.GraphicsPlaceholder')), hGOs);
    hasValidNonemptyPosition = cellfun(@(pos) not(isempty(pos)) & isequal(size(pos), [1, 4]) & isnumeric(pos), cellVectPosArrNrm);
    tryReposition = isValidNonPlaceholder(:) & hasValidNonemptyPosition(:);
    arrayfun(@(idx) set_position(hGOs(idx), cellVectPosArrNrm{idx}, 'Normalized'), find(tryReposition));
end