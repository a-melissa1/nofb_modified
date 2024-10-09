function DrawFixationCross_phat(sci, pp, centerX, centerY)
    Screen('DrawLine', sci.wnd, [100,50,50], centerX, centerY - pp.fc_size, centerX, centerY + pp.fc_size, pp.fc_lw);% vertical line
    Screen('DrawLine', sci.wnd, [100,50,50], centerX - pp.fc_size, centerY, centerX + pp.fc_size, centerY, pp.fc_lw);% horizontal line
end