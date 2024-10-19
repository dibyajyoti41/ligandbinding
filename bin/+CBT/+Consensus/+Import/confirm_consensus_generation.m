function [quitConsensus] = confirm_consensus_generation(numBarcodes)
    continueGenPrompt = sprintf('Generate consensus for %d barcodes?', numBarcodes);
    continueGenChoice = questdlg(continueGenPrompt, 'Consensus Generation Confirmation', 'Continue', 'Continue');
    quitConsensus = not(strcmp(continueGenChoice, 'Continue'));
end