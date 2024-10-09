% Draws 4 pacmen centered around the specified position.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function Pacman(sci, pp, centerX, centerY, pac_side, vert)
    prect = [0,0,pp.pac_radius*2,pp.pac_radius*2];

    y_pac_offset = (pp.VernierLength)/2 + (pp.pac_mouth)/2;

    if vert==1
        p1 = CenterRectOnPoint(prect,centerX-pp.x_pac_offset,centerY-y_pac_offset);
        p2 = CenterRectOnPoint(prect,centerX-pp.x_pac_offset,centerY+y_pac_offset);
        p3 = CenterRectOnPoint(prect,centerX+pp.x_pac_offset,centerY-y_pac_offset);
        p4 = CenterRectOnPoint(prect,centerX+pp.x_pac_offset,centerY+y_pac_offset);
    else
        p1 = CenterRectOnPoint(prect,centerX-y_pac_offset,centerY-pp.x_pac_offset);
        p2 = CenterRectOnPoint(prect,centerX-y_pac_offset,centerY+pp.x_pac_offset);
        p3 = CenterRectOnPoint(prect,centerX+y_pac_offset,centerY-pp.x_pac_offset);
        p4 = CenterRectOnPoint(prect,centerX+y_pac_offset,centerY+pp.x_pac_offset);
    end

    Screen('FillOval', sci.wnd, pp.LineColor, p1);
    Screen('FillOval', sci.wnd, pp.LineColor, p2);
    Screen('FillOval', sci.wnd, pp.LineColor, p3);
    Screen('FillOval', sci.wnd, pp.LineColor, p4);

    % pacman mouth
    if vert==1
        mrect = [0,0,pp.pac_radius,pp.pac_mouth];

        if pac_side % inside
            m1 = CenterRectOnPoint(mrect,centerX-pp.x_pac_offset+(pp.pac_radius/2),centerY-y_pac_offset);
            m2 = CenterRectOnPoint(mrect,centerX-pp.x_pac_offset+(pp.pac_radius/2),centerY+y_pac_offset);
            m3 = CenterRectOnPoint(mrect,centerX+pp.x_pac_offset-(pp.pac_radius/2),centerY-y_pac_offset);
            m4 = CenterRectOnPoint(mrect,centerX+pp.x_pac_offset-(pp.pac_radius/2),centerY+y_pac_offset);
        else % outside
            m1 = CenterRectOnPoint(mrect,centerX-pp.x_pac_offset-(pp.pac_radius/2),centerY-y_pac_offset);
            m2 = CenterRectOnPoint(mrect,centerX-pp.x_pac_offset-(pp.pac_radius/2),centerY+y_pac_offset);
            m3 = CenterRectOnPoint(mrect,centerX+pp.x_pac_offset+(pp.pac_radius/2),centerY-y_pac_offset);
            m4 = CenterRectOnPoint(mrect,centerX+pp.x_pac_offset+(pp.pac_radius/2),centerY+y_pac_offset);
        end
    else
        mrect = [0,0, pp.pac_mouth, pp.pac_radius];
        
        if pac_side % inside
            m1 = CenterRectOnPoint(mrect, centerX-y_pac_offset, centerY-pp.x_pac_offset+(pp.pac_radius/2));
            m2 = CenterRectOnPoint(mrect, centerX-y_pac_offset, centerY+pp.x_pac_offset-(pp.pac_radius/2));
            m3 = CenterRectOnPoint(mrect, centerX+y_pac_offset, centerY-pp.x_pac_offset+(pp.pac_radius/2));
            m4 = CenterRectOnPoint(mrect, centerX+y_pac_offset, centerY+pp.x_pac_offset-(pp.pac_radius/2));
        else % outside
            m1 = CenterRectOnPoint(mrect, centerX-y_pac_offset, centerY-pp.x_pac_offset-(pp.pac_radius/2));
            m2 = CenterRectOnPoint(mrect, centerX-y_pac_offset, centerY+pp.x_pac_offset+(pp.pac_radius/2));
            m3 = CenterRectOnPoint(mrect, centerX+y_pac_offset, centerY-pp.x_pac_offset-(pp.pac_radius/2));
            m4 = CenterRectOnPoint(mrect, centerX+y_pac_offset, centerY+pp.x_pac_offset+(pp.pac_radius/2));
        end
    end
    
    Screen('FillRect', sci.wnd, pp.BgLum, m1);
    Screen('FillRect', sci.wnd, pp.BgLum , m2);
    Screen('FillRect', sci.wnd, pp.BgLum , m3);
    Screen('FillRect', sci.wnd, pp.BgLum , m4);
end
