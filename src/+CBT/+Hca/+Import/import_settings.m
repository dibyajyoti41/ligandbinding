function [sets] = import_settings(setsName)
    % import settings from file

    % settings
    sets = ini2struct( setsName);
    % model settings. Exactly the one used in Contig Assembly paper. But should
    % add other models too.
    
    import CBT.Hca.Core.Theory.choose_cb_model;
    [sets.model ] = choose_cb_model(sets.theoryGen.method,sets.model.pattern);

    
    
  

end

