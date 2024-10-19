function [indices] = reorient_circ_indices(indices, flipTF, circShift)
    indices = circshift(indices, circShift, 2);
    if flipTF
        indices = fliplr(indices);
    end
end