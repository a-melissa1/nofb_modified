% Draws the gapped flankers around the vernier.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function Flankers_gaps(sci, pp, centerX, centerY, fl_offset, n_flankers, lw, vert)
    flen = pp.VernierLength/2;
    fl_pos = fl_offset;
    outer_fl_len = pp.flanker_len - flen - pp.pac_mouth;
    y_outer_fl_start = flen + pp.pac_mouth;
    for i=1:n_flankers
        if vert == 1
            % Top Row
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-fl_pos, centerY+y_outer_fl_start, centerX-fl_pos, centerY+y_outer_fl_start+outer_fl_len, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX+fl_pos, centerY+y_outer_fl_start, centerX+fl_pos, centerY+y_outer_fl_start+outer_fl_len, lw);% right line
            % Mid Row
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-fl_pos, centerY - flen, centerX-fl_pos, centerY + flen, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX+fl_pos, centerY - flen, centerX+fl_pos, centerY + flen, lw);% right line
            % Bottom Row
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-fl_pos, centerY - y_outer_fl_start, centerX-fl_pos, centerY - y_outer_fl_start - outer_fl_len, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX+fl_pos, centerY - y_outer_fl_start, centerX+fl_pos, centerY - y_outer_fl_start - outer_fl_len, lw);% right line
        else
            % Right Row
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX+y_outer_fl_start, centerY-fl_pos, centerX+y_outer_fl_start+outer_fl_len, centerY-fl_pos, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX+y_outer_fl_start, centerY+fl_pos, centerX+y_outer_fl_start+outer_fl_len, centerY+fl_pos, lw);% right line
            % Mid Row
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-flen, centerY-fl_pos, centerX+flen, centerY-fl_pos, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-flen, centerY+fl_pos, centerX+flen, centerY+fl_pos, lw);% right line
            % Left Row
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-y_outer_fl_start, centerY-fl_pos, centerX-y_outer_fl_start-outer_fl_len, centerY-fl_pos, lw);% left line
            Screen('DrawLine', sci.wnd, pp.LineColor, centerX-y_outer_fl_start, centerY+fl_pos, centerX-y_outer_fl_start-outer_fl_len, centerY+fl_pos, lw);% right line


        end
        fl_pos = fl_pos + fl_offset;
    end
end