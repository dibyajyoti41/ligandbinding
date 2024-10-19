classdef (Sealed) AppDataPoolMgr < handle
    % APPDATAPOOLMGR - Application data pool manager helps manage data
    %   pools shared by the various components of a running app
    %
    % Note: An instance of APPDATAPOOLMGR should be retrieved
    %   using get_instance since it will be defined as a persistent
    %   variable using the singleton design pattern.
    %
    % http://mathworks.com/help/matlab/ref/persistent.html
    % http://mathworks.com/help/matlab/matlab_oop/controlling-the-number-of-instances.html
    %
    % Authors:
    %   Saair Quaderi
    
    properties (Constant)
        Version = [0, 0, 1];
    end
    properties (Constant, Access = private)
        DefaultVerbose = false;
    end
    properties (SetAccess = private)
        DataMap %containers.Map mapping strings (data pool ids) to containers.Map which maps data item ids to data item values
        UUID % unique ID for debugging purposes
        Verbose = Fancy.AppMgr.AppDataPoolMgr.DefaultVerbose
    end
    methods
        function [appDataPoolMgr] = AppDataPoolMgr(uuid, verboseTF)
            % APPDATAPOOLMGR - Constructor for app data pool manager
            %   objects
            %
            % Inputs:
            %  uuid (optional, defaults to randomly generated UUID)
            %    a string containing a randomly generated universally
            %    unique identifier
            %    (see wikipedia page for Universally_unique_identifier)
            %  verboseTF (optional, defaults to DefaultVerbose property)
            %    if true, more information may be printed to the console
            %    about what is going on than would otherwise be printed
            %    (potentially useful for debugging)
            %
            % Outputs:
            %   appDataPoolMgr
            %     class instance of app data pool manager
            %
            % Authors:
            %   Saair Quaderi
            
            if (nargin < 1)
                uuid = [];
            end
            if (nargin < 2)
                verboseTF = [];
            end
            if isempty(uuid)
                uuid = char(java.util.UUID.randomUUID);
            end
            if isempty(verboseTF)
                verboseTF = Fancy.AppMgr.AppDataPoolMgr.DefaultVerbose;
            end
            appDataPoolMgr.DataMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
            appDataPoolMgr.UUID =  uuid;
            appDataPoolMgr.Verbose = verboseTF;

            if appDataPoolMgr.Verbose
                fprintf('Generated new app data pool manager: %s\n', uuid);
            end
        end
    end
    methods (Static)
        function [dataPoolIDs, dataPoolsItemCounts] = get_pool_ids()
            % GET_POOL_IDS - returns data pool IDs and how man items each
            %   pool contains
            %
            % Outputs:
            %   dataPoolIDs
            %     cell column vector containing pool ID strings
            %   dataPoolsItemCounts
            %     column vector containing number of items in each
            %     of the pools
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            dataPoolIDs = keys(appDataMap);
            dataPoolIDs = dataPoolIDs(:);
            if nargout > 1
                dataPoolsItemCounts = cellfun(@(dataPoolID) size(appDataMap(dataPoolID), 1), dataPoolIDs);
            end
        end
        function [dataItemIDs] = get_pool_data_item_ids(dataPoolID)
            % GET_POOL_DATA_ITEM_IDS - returns data item IDs (keys for
            %   mapped values) for a data pool where the data pool is
            %   for which to get keys is specified by its data pool ID
            %
            % Inputs:
            %   dataPoolID
            %     a string containing the ID for a data pool (managed by
            %     the app data pool manager instance returned
            %     by AppDataPoolMgr.get_instance)
            %
            %
            % Outputs:
            %   dataItemIDs
            %     cell column vector containing data item IDs for the items
            %     in the specified data pool (empty if data pool is not
            %     found)
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            if appDataMap.isKey(dataPoolID)
                datapoolMap = appDataMap(dataPoolID);
                dataItemIDs = keys(datapoolMap);
                dataItemIDs = dataItemIDs(:);
            else
                dataItemIDs = cell(0,1);
            end
        end
        
        function [dataItems, dataItemIDs, dataItemsFoundMask] = get_data_items(dataPoolID, dataItemIDs)
            % GET_POOL_DATA_ITEM_IDS - returns data item information
            %   for a data pool where the data pool is specified by its
            %   data pool ID (allowing request of only specific data items
            %   in the pool by passing in associated dataItemIDs)
            %
            % Inputs:
            %   dataPoolID
            %     a string containing the ID for a data pool (managed by
            %     the app data pool manager instance returned
            %     by AppDataPoolMgr.get_instance)
            %   dataItemIDs (optional; defaults to all data item IDs for
            %       data pool)
            %     cell column vector containing data item IDs for the items
            %     in the specified data pool (empty if data pool is not
            %     found)
            %
            % Outputs:
            %   dataItems
            %     cell column vector containing the data items associated
            %     with the dataItemIDs (cell is empty for cases where data
            %     item ID was not found as a key in the data pool that
            %     was specified)
            %   dataItemIDs
            %     cell column vector containing data item IDs for the items
            %     in the data pool (same as the optional input if it was
            %     provided)
            %   dataItemsFoundMask
            %     logical column vector with true for data item IDs which
            %     were found and false for ones which weren't found (useful
            %     if an empty array as specified for missing values in
            %     dataItems is an actual non-missing data item value which
            %     needs to be stored and distinguished from missing data
            %     items
            %
            % Authors:
            %   Saair Quaderi
            
            if nargin < 2
                dataItemIDs = AppDataPoolMgr.get_pool_data_item_ids(dataPoolID);
            end
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            numItems = length(dataItemIDs);
            dataItemsFoundMask = false(numItems, 1);
            dataItems = cell(numItems, 1);
            dataPoolFoundTF = appDataMap.isKey(dataPoolID);
            if not(dataPoolFoundTF)
                return;
            end
            if dataPoolFoundTF
                datapoolMap = appDataMap(dataPoolID);
                for itemNum = 1:numItems
                    dataItemID = dataItemIDs{itemNum};
                    dataItemsFoundMask(itemNum) = datapoolMap.isKey(dataItemID);
                    if dataItemsFoundMask(itemNum)
                        dataItems{itemNum} = datapoolMap(dataItemID);
                    end
                end
            end
        end
        function [dataItem, dataItemID, dataItemFoundTF] = get_data_item(dataPoolID, dataItemID)
            % GET_POOL_DATA_ITEM_IDS - returns data item information
            %   for a data item with a specific ID in a specific data pool
            %
            % Inputs:
            %   dataPoolID
            %     a string containing the ID for a data pool (managed by
            %     the app data pool manager instance returned
            %     by AppDataPoolMgr.get_instance)
            %   dataItemID 
            %     the ID (key) for the data item to retreive from the data
            %     pool
            %
            % Outputs:
            %   dataItem
            %     the value of the data item associated with the data item
            %     ID (empty if data item was not found)
            %   dataItemID
            %     same as the input value (only an output to be consistent
            %     with get_data_items functionality)
            %   dataItemFoundTF
            %     logical value which is true if a data item with the
            %     specified ID was found and false otherwise
            %
            % Authors:
            %   Saair Quaderi
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            [dataItems, ~, dataItemsFoundMask] = appDataPoolMgr.get_data_items(dataPoolID, {dataItemID});
            dataItem = dataItems{1};
            dataItemFoundTF = dataItemsFoundMask(1);
        end
        function [dataItemsRemoved, dataItemsRemovedMask] = remove_data_items(dataPoolID, dataItemIDsToRemove)
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            dataTypeFound = appDataMap.isKey(dataPoolID);
            if not(dataTypeFound)
                error('Unknown data pool');
            end
            
            [dataItemsRemoved, ~, dataItemsRemovedMask] = AppDataPoolMgr.get_data_items(dataPoolID, dataItemIDsToRemove);
            
            datapoolMap = appDataMap(dataPoolID);
            datapoolMap.remove(dataItemIDsToRemove(dataItemsRemovedMask));
        end
        
        function [dataItem, dataItemsRemovedTF] = remove_data_item(dataPoolID, dataItemID)
            import Fancy.AppMgr.AppDataPoolMgr;
            [dataItems, dataItemsRemovedTF] = AppDataPoolMgr.remove_data_items(dataPoolID, {dataItemID});
            dataItem = dataItems{1};
        end
        
        function [] = update_data_items(dataPoolID, dataItemIDs, dataItems)
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            numItems = numel(dataItemIDs);
            if length(dataItems) < numItems
                error('Mismatch is ID and item count');
            end
            appDataMap = appDataPoolMgr.DataMap;
            dataPoolFound = appDataMap.isKey(dataPoolID);
            if not(dataPoolFound)
                error('Unknown data pool');
            end
            datapoolMap = appDataMap(dataPoolID);
            for itemNum = 1:numItems
                datapoolMap(dataItemIDs{itemNum}) = dataItems{itemNum};
            end
            appDataMap(dataPoolID) = datapoolMap;
            appDataPoolMgr.DataMap = appDataMap;
        end
        function [] = update_data_item(dataPoolID, dataItemID, dataItem)
            import Fancy.AppMgr.AppDataPoolMgr;
            AppDataPoolMgr.update_data_items(dataPoolID, {dataItemID}, {dataItem})
        end
        function [hasDataPool] = has_data_pool(dataPoolID)
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            hasDataPool = appDataMap.isKey(dataPoolID);
        end
        
        function [] = create_new_data_pool(dataPoolID, keyType, valueType)
            if nargin < 2
                keyType = 'char';
            end
            if nargin < 3
                valueType = 'any';
            end
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            if appDataMap.isKey(dataPoolID)
                error('Data pool already exists');
            end
            appDataMap(dataPoolID) = containers.Map('KeyType', keyType, 'ValueType', valueType);
            appDataPoolMgr.DataMap = appDataMap;
        end
        function [datapoolMap] = remove_data_pool(dataPoolID)
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            if appDataMap.isKey(dataPoolID)
                datapoolMap = appDataMap(dataPoolID);
            end
            appDataMap.remove(dataPoolID);
        end
        function [] = clear_data_pool(dataPoolID)
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            appDataMap = appDataPoolMgr.DataMap;
            if appDataMap.isKey(dataPoolID)
                datapoolMap = appDataMap(dataPoolID);
                dataItemIDs = keys(datapoolMap);
                dataItemIDs = dataItemIDs(:);
                datapoolMap.remove(dataItemIDs);
            end
        end
        function [] = total_reset()
            import Fancy.AppMgr.AppDataPoolMgr;
            appDataPoolMgr = AppDataPoolMgr.get_instance();
            delete(appDataPoolMgr);
            AppDataPoolMgr.get_instance();
        end
        function [appDataPoolMgr] = get_instance(varargin)
            persistent persistentAppDataPoolMgr;
            import Fancy.AppMgr.AppDataPoolMgr;
            if isempty(persistentAppDataPoolMgr) || not(isvalid(persistentAppDataPoolMgr))
                persistentAppDataPoolMgr = AppDataPoolMgr();
                appDataPoolMgr = persistentAppDataPoolMgr;
            else
                appDataPoolMgr = persistentAppDataPoolMgr;
            end
        end
    end
end