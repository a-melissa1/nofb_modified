% Draws a fixation cross
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function DrawFixationCross_basic(sci, pp, centerX, centerY)
    Screen('DrawLine', sci.wnd, pp.fc_lc, centerX, centerY - pp.fc_size, centerX, centerY + pp.fc_size, pp.fc_lw);% vertical line
    Screen('DrawLine', sci.wnd, pp.fc_lc, centerX - pp.fc_size, centerY, centerX + pp.fc_size, centerY, pp.fc_lw);% horizontal line
end