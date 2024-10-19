function [ lm ] = create_button_set( lm, fun )
% creates a button for a given function fun

import Fancy.UI.FancyList.FancyListMgrBtnSet;
flmbs2 = FancyListMgrBtnSet();
flmbs2.NUM_BUTTON_COLS = 1;
flmbs2.add_button(fun());

lm.add_button_sets(flmbs2);

end

