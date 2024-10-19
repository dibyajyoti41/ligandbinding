classdef DataExporter < handle
    % DATAEXPORTER - Functions for exporting data from DBM
    
    properties (Access = private)
        dbmODW
        dbmOSW
    end
    methods
        function [dbmDE] = DataExporter(dbmODW, dbmOSW)
            import OldDBM.General.DataWrapper;
            import OldDBM.General.SettingsWrapper;
            validateattributes(dbmODW, {'OldDBM.General.DataWrapper'}, {'scalar'}, 1);
            validateattributes(dbmOSW, {'OldDBM.General.SettingsWrapper'}, {'scalar'}, 1);
            dbmDE.dbmODW = dbmODW;
            dbmDE.dbmOSW = dbmOSW;
        end

        function [] = export_dbm_session_struct_mat(dbmDE, defaultOutputDirpath)
            import OldDBM.General.Export.export_dbm_session_struct_mat;
            if nargin < 2
                defaultOutputDirpath = dbmDE.dbmOSW.get_default_export_dirpath('session');
            end
            
            export_dbm_session_struct_mat(dbmDE.dbmODW, dbmDE.dbmOSW, defaultOutputDirpath);
        end
        
        function [] = export_raw_kymos(dbmDE, defaultOutputDirpath)
            import OldDBM.General.Export.export_raw_kymos;
            if nargin < 2
                defaultOutputDirpath = dbmDE.dbmOSW.get_default_export_dirpath('raw_kymo');
            end
            
            export_raw_kymos(dbmDE.dbmODW, defaultOutputDirpath);
        end

        function [] = export_aligned_kymos(dbmDE, defaultOutputDirpath)
            import OldDBM.General.Export.export_aligned_kymos;
            if nargin < 3
                defaultOutputDirpath = dbmDE.dbmOSW.get_default_export_dirpath('aligned_kymo');
            end
            
            export_aligned_kymos(dbmDE.dbmODW, defaultOutputDirpath);
        end

        function [] = export_aligned_kymo_time_avs(dbmDE, defaultOutputDirpath)
            import OldDBM.General.Export.export_aligned_kymo_time_avs;
            if nargin < 3
                defaultOutputDirpath = dbmDE.dbmOSW.get_default_export_dirpath('aligned_kymo_time_avg');
            end
            
            export_aligned_kymo_time_avs(dbmDE.dbmODW, defaultOutputDirpath);
        end

        function [kymoStructs] = extract_kymo_structs(dbmDE)
            import OldDBM.General.Export.extract_kymo_structs;
            
            kymoStructs = extract_kymo_structs(dbmDE.dbmODW);
        end
        
        function [experimentCurveStructs] = extract_experimental_structs(dbmDE)
            import OldDBM.General.Export.extract_experimental_structs;
            
            experimentCurveStructs = extract_experimental_structs(dbmDE.dbmODW);
        end
    end
end