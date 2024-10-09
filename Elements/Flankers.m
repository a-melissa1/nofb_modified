% Draws the long flankers around the vernier.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function Flankers(sci, pp, centerX, centerY, fl_offset, n_flankers, flen, lw, vert)
    fl_pos = fl_offset;
    for i=1:n_flankers
        if vert == 1
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-fl_pos, centerY - flen, centerX-fl_pos, centerY + flen, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX+fl_pos, centerY - flen, centerX+fl_pos, centerY + flen, lw);% right line
        else
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-flen, centerY-fl_pos, centerX+flen, centerY-fl_pos, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-flen, centerY+fl_pos, centerX+flen, centerY+fl_pos, lw);% right line
        end
        fl_pos = fl_pos + fl_offset;
    end
end