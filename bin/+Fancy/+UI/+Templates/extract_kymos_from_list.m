function [ hcaSessionStruct ] = extract_kymos_from_list( lm, timeFramesNr )
    % extract_kymos_from_list
    
    [selectedItems, ~] = get_selected_list_items(lm);

    % put the selected items in a structure
    kymoStructs = cell(1,size(selectedItems,1));
    kymoNames =  cell(1,size(selectedItems,1));
    for it=1:size(selectedItems,1)   
        kymoStructs{end-it+1} = selectedItems{it,2};
        kymoNames{end-it+1} =  selectedItems{it,1};
    end
    
    % define session structure. Is this the best place to define
    % it?
    hcaSessionStruct = struct();

     % put kymos in the session structure
    import CBT.Hca.UI.edit_kymographs_fun;
    hcaSessionStruct = edit_kymographs_fun(hcaSessionStruct,kymoStructs,kymoNames,timeFramesNr);
            
            


end

