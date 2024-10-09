% Draws all the distractors/flankers of the stimulus.
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.

function DrawFlankers(sci, pp, centerX, centerY, cond, vert)
    
    flen = pp.flanker_len;
    lw = pp.LineWidth;

    if cond ==1
        % 1 just vernier
        nflankers = 0;
        % no flankers
    elseif cond==2
        % 2 short flankers
        nflankers = 1;
    elseif  cond==3
        nflankers=8;

    end

    if nflankers ~= 0
        Flankers(sci, pp, centerX, centerY, pp.fl_offset, nflankers, flen, lw, vert);
    end



end