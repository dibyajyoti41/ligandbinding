classdef FancyListMgrBtn < handle
    % FANCYLISTMGRBTN - Fancy List Manager Button
    properties
        InitButtonText = ''
        InitCallbackFn = @Fancy.UI.FancyList.FancyListMgrBtn.nop
        InitEnable = 'on'
        Handle = gobjects(0)
    end
    
    methods
        function [flmb] = FancyListMgrBtn(initButtonText, initCallbackFn, initEnable)
            if nargin >= 1
                flmb.InitButtonText = initButtonText;
            end
            if nargin >= 2
                flmb.InitCallbackFn = initCallbackFn;
            end
            if nargin >= 3
                flmb.InitEnable = initEnable;
            end
        end
        
        function [] = instantiate(flmb, flm)
            hButton = flmb.Handle;
            if isa(hButton, 'matlab.graphics.GraphicsPlaceholder')
                hButton = uicontrol('Parent', flm.ParentHandle, 'Style', 'pushbutton', 'String', flmb.InitButtonText, 'Enable', flmb.InitEnable);
                iptaddcallback(hButton, 'Callback', {flmb.InitCallbackFn, flm});
                flmb.Handle = hButton;
            end
        end
    end
    
    methods(Static)
        function [] = nop(~, ~, varargin)
        end
    end
end