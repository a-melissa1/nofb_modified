% Generates a horizontal color gradient. The gradient is gaussian in shape, with
%a given center and sd and color for the middle and edges.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function combined = colorgrad_horz(w, h, h_center, windowsz, sigma, top_col, bottom_col)
    x = [1:1:h];
    left_edge = h_center-floor(windowsz/2);
    y = normpdf(x,left_edge,sigma);
    oneline = [y(1:left_edge), ones(1, windowsz)*max(y), y(left_edge+1:end)];
    oneline = oneline(1:h);
    oneline = oneline/ max(oneline);
    oneline = oneline*(top_col-bottom_col);
    oneline = (oneline+bottom_col);
    oneline = transpose(oneline);
    combined = repmat(oneline,[1,w,3]);
    combined = cast(combined,"uint8");
end