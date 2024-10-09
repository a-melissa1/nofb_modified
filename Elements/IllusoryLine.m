% Draws the bar that connects the two pacmen.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function IllusoryLine(sci, pp, centerX,centerY, pac_side, vert)
    y_pac_offset = (pp.VernierLength)/2 + (pp.pac_mouth)/2;
    if vert
        lrect = [0,0,pp.x_pac_offset*2,pp.pac_mouth];
        if ~pac_side
            lrect = [0,0,(pp.x_pac_offset*2)-(pp.pac_radius*2),pp.pac_mouth];
        end
        l1 = CenterRectOnPoint(lrect,centerX,centerY-y_pac_offset);
        l2 = CenterRectOnPoint(lrect,centerX,centerY+y_pac_offset);
        Screen('FillRect', sci.wnd, pp.BgLum , l1);
        Screen('FillRect', sci.wnd, pp.BgLum , l2);
    else
        lrect = [0,0,pp.pac_mouth,pp.x_pac_offset*2];
        if ~pac_side
            lrect = [0,0,pp.pac_mouth, (pp.x_pac_offset*2)-(pp.pac_radius*2)];
        end
        l1 = CenterRectOnPoint(lrect, centerX-y_pac_offset, centerY);
        l2 = CenterRectOnPoint(lrect, centerX+y_pac_offset, centerY);
        Screen('FillRect', sci.wnd, pp.BgLum , l1);
        Screen('FillRect', sci.wnd, pp.BgLum , l2);
    end
end