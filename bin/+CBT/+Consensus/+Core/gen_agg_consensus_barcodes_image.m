function [consensusImage] = gen_agg_consensus_barcodes_image(individualBarcodes, individualBitmasks, consensusBarcode, consensusBitmask)
    import Barcoding.Visualizing.colorize_bitmasked_barcode;

    hueGradients.individualBarcodes.inclusions = [0 0 0; 0 1 1]; % black to cyan
    hueGradients.individualBarcodes.exclusions = [1 0 0; 1 1 1]; % red to white

    hueGradients.consensusBarcodes.inclusions = [0 0 0; 0 1 0]; % black to green
    hueGradients.consensusBarcodes.exclusions = [1 0 0; 1 1 1]; % red to white


    numBarcodes = size(individualBarcodes, 1);
    images = cell(numBarcodes, 1);
    for bacodeNum=1:numBarcodes
        alignedBarcode = individualBarcodes{bacodeNum};
        alignedBitmask = individualBitmasks{bacodeNum};
        images{bacodeNum} = colorize_bitmasked_barcode( ...
            alignedBarcode, ...
            alignedBitmask, ...
            hueGradients.individualBarcodes.inclusions, ...
            hueGradients.individualBarcodes.exclusions);
    end

    consensusImage = colorize_bitmasked_barcode(...
        consensusBarcode, ...
        consensusBitmask, ...
        hueGradients.consensusBarcodes.inclusions, ...
        hueGradients.consensusBarcodes.exclusions);
    consensusImage = cell2mat([images; {consensusImage}]);
end