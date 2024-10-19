function [clusterThreshold, quitConsensus] = get_cluster_threshold(consensusSettings)
    % get_cluster_threshold - Retrieve the threshold for clustering barcodes
    clusterThreshold = consensusSettings.barcodeClusterLimit;
    quitConsensus = false;
    if (consensusSettings.promptForBarcodeClusterLimit == 0)
       return 
    end
    % Set a threshold for grouping barcodes.
    answer = inputdlg(...
        {'Limit for when barcodes are considered part of same cluster:'},... % prompt
        'Limit for similarity',... % dialog title
        1,... % number of lines
        {num2str(clusterThreshold)}... % default value
    );
    if isempty(answer)
        quitConsensus = true;
    else 
        clusterThreshold = str2double(answer{1});
    end
end