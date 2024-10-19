function [barcodeKeys, alignedBarcodes, alignedBarcodeBitmasks, barcodes, barcodeBitmasks, flipTFs, circShifts] = extract_consensus_aligned_barcodes(consensusKey, barcodeStructsMap)
    import CBT.Consensus.Helper.extract_consensus_components;
    import Barcoding.Reorienting.reorient_barcode_and_bitmask;

    [barcodeKeys, barcodes, barcodeBitmasks, flipTFs, circShifts] = extract_consensus_components(consensusKey, barcodeStructsMap);
    numBarcodes = size(barcodeKeys, 1);
    alignedBarcodes = cell(numBarcodes, 1);
    alignedBarcodeBitmasks = cell(numBarcodes, 1);
    for barcodeNum=1:numBarcodes
        [alignedBarcodes{barcodeNum}, alignedBarcodeBitmasks{barcodeNum}] = reorient_barcode_and_bitmask(barcodes{barcodeNum}, barcodeBitmasks{barcodeNum}, flipTFs(barcodeNum), circShifts(barcodeNum));
    end
end