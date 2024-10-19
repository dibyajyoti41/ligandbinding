classdef FancyListMgrBtnSet < handle
    % FANCYLISTMGRBTNSET - Fancy List Manager Button SET
    properties
        BUTTON_HEIGHT_PX = 30
        BUTTON_WIDTH_PX = 0
        
        NUM_BUTTON_COLS = 1
        
        Buttons = cell(0, 1);
    end
    
    methods
        function [flmbs] = FancyListMgrBtnSet()
        end
        function [] = add_button(flmbs, button)
            flmbs.Buttons = [flmbs.Buttons; {button}];
        end
    end
end