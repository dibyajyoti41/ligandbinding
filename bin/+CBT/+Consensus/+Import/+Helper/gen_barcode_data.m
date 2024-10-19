function [barcodeGenData] = gen_barcode_data(alignedKymo)
    skipDoubleTanhAdjustment = true;
    import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
    edgeDetectionSettings = get_default_edge_detection_settings(skipDoubleTanhAdjustment);

    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    [leftEdgeIdxs, rightEdgeIdxs, alignedKymoMoleculeMask] = approx_main_kymo_molecule_edges(alignedKymo, edgeDetectionSettings);
    
    leftEdgeIdx = round(nanmean(leftEdgeIdxs));
    rightEdgeIdx = round(nanmean(rightEdgeIdxs));
    
    if rightEdgeIdx <= leftEdgeIdx
        disp('ah');
    end
    
    % Determine indices for rotated barcode with background cropped out
    barcodeIdxs = leftEdgeIdx:rightEdgeIdx;
    rawBarcode = nanmean(alignedKymo, 1);
    
    bgIndices = true(1, size(alignedKymo, 1));
    bgIndices(barcodeIdxs) = 0;

    nonBarcodeVals = alignedKymo(:, bgIndices);
    nonBarcodeVals = nonBarcodeVals(:);
    nonBarcodeDistFit = fitdist(nonBarcodeVals, 'Normal');
    barcodeGenData.bgMeanApprox =  nonBarcodeDistFit.mu;
    barcodeGenData.bgStdApprox = nonBarcodeDistFit.sigma;
    
    rawBg = nanmean(rawBarcode(bgIndices));
  
    barcodeGenData.rawBarcode = rawBarcode(barcodeIdxs);
    
    barcodeGenData.rawBarcodeLeftEdgeIndex = leftEdgeIdx;
    barcodeGenData.rawBarcodeRightEdgeIndex = rightEdgeIdx;
    barcodeGenData.rawBg = rawBg;
end
