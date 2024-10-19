function [reorientedBarcode, reorientedBarcodeBitmask] = reorient_barcode_and_bitmask(initBarcode, initBarcodeBitmask, flipTF, circShift)
    import Barcoding.Reorienting.reorient_circ_indices;

    [indices] = reorient_circ_indices(1:length(initBarcode), flipTF, circShift);
    reorientedBarcode = initBarcode(indices);
    reorientedBarcodeBitmask = initBarcodeBitmask(indices);
end