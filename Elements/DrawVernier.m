% Draws the vernier.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function DrawVernier(sci, pp, x_pix, y_pix, signed_ofs_pix, oriDeg)

    %--- Coordinates for upper vertical segment - the lower segment is point-symmetrical WRT (0,0).
    x =  signed_ofs_pix*[1,1];
    y =  [pp.VernierGap/2, pp.VernierLength/2];
    xyRot = [cosd(oriDeg), sind(oriDeg); -sind(oriDeg), cosd(oriDeg)];
    Screen('DrawLines', sci.wnd, xyRot*[x,-x; y,-y], pp.LineWidth, pp.LineColor, [x_pix, y_pix], 2);
    return
end