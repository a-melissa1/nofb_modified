% FOR THE EVENMOORE Experiment.
% Writes important Data to the DV file
%
%
%################# H I S T O R Y #####################
% 26.04.2024 (LS = Lisa Schwetlick, LPSY@EPFL):
%   * First version.


function writeppp2DV(p, pp, adaptiveUsed, myPEST, myQUEST, pest, practiceTestLevels, dontRoundJitterX)
    lpsy.writeDvPmtr('Screen distance', p.ScreenDist, 'cm', 0);
    if isscalar(p.PracticeTrials)
        lpsy.writeDvPmtr('Practice trials',p.PracticeTrials,[]);
        % if ~isempty(practiceTestLevels) && ~isnan(practiceTestLevels(1))
        %     lpsy.writeDvPmtr('Practice test levels', sprintf('[%.1f,...,%.1f] %s', pest.pp2user(practiceTestLevels([1,end])),pest.userUnit));
        % end
    else
        lpsy.writeDvPmtr('Practice trials',length(p.PracticeTrials),[]);
        %lpsy.writeDvPmtr('Practice test levels', p.PracticeTrials, pest.userUnit);
    end
    lpsy.writeDvPmtr('Start trial repetitions (min/max)',p.StartTrialReps,'');
    lpsy.writeDvPmtr('Trials',p.Trials,[]);
    lpsy.writeDvPmtr('Volume of negative feedback beep', p.NegFeedbackVol, '%');
    lpsy.writeDvPmtr('Volume of warning beep (e.g. early response)', p.WarningsVol, '%');
    lpsy.writeDvPmtr('Vernier duration', p.VernierDuration, 'ms', 0, lpsy.frame2msec(pp.VernierDuration));
    lpsy.writeDvPmtr('Vernier length (total, including gap)', p.VernierLength, 'arcsec', 0, lpsy.pix2arcsec(pp.VernierLength));
    lpsy.writeDvPmtr('Vernier vertical gap size', p.VernierGap, 'arcsec', 0, lpsy.pix2arcsec(pp.VernierGap));
    lpsy.writeDvPmtr('Vernier offset (possibly PEST start)', p.VernierOffset, 'arcsec', 0, lpsy.pix2arcsec(pp.VernierOffset));
    lpsy.writeDvPmtr('Vernier orientation', p.VernierOri, 'degree');


    if (adaptiveUsed)
        lpsy.writeDvPmtr('PEST/QUEST: Target parameter',p.AdaptivePmtr);
        pr = p.TestValueRange;
        if isempty(pr)
            %--- Report "-1", because we cannot report "empty".
            pr = [-1,-1];
        end
        if ~isempty(myPEST)
            lpsy.writeDvPmtr('PEST: Start step width (WRT start value)',p.PestStartStep,'%');
            lpsy.writeDvPmtr('PEST: Stop step width (WRT current value)',p.PestStopStep,'%');
            lpsy.writeDvPmtr('PEST: WALD constant (adaptiveness)',p.PestWALD,[]);
            lpsy.writeDvPmtr('PEST: Smallest test value', pr(1), pest.userUnit, 0, pest.pp2user(pp.TestValueRange(1)));
            lpsy.writeDvPmtr('PEST: Biggest test value', pr(2), pest.userUnit, 0, pest.pp2user(pp.TestValueRange(2)));
        elseif ~isempty(myQUEST)
            lpsy.writeDvPmtr('QUEST: Threshold SD (sigma of prior p.d.f.)',p.QuestThreshPriorSD,'dB');
            lpsy.writeDvPmtr('QUEST: Slope (CumGauss sigma)',p.QuestSigma,'dB');
            lpsy.writeDvPmtr('QUEST: Lapse rate', 100*p.QuestLapseRate, '%',1);        
            lpsy.writeDvPmtr('QUEST: Smallest test value', pr(1), pest.userUnit, 0, pest.pp2user(pp.TestValueRange(1)));
            lpsy.writeDvPmtr('QUEST: Biggest test value', pr(2), pest.userUnit, 0, pest.pp2user(pp.TestValueRange(2)));
        end  
    end
    lpsy.writeDvPmtr('Response timeout (WRT vernier onset, 0:infinite)', p.ResponseTimeout, 's');    
    if numel(p.ITI)>1
        lpsy.writeDvPmtr('ITI (time from response to next vernier onset)',sprintf('[%.2f, %.2f] s', p.ITI));    
    else
        lpsy.writeDvPmtr('ITI (time from response to next vernier onset)',p.ITI,'s');    
    end
    
    lpsy.writeDvPmtr('Line rendering mode', 'aaL');    
    

    lpsy.writeDvPmtr('Line width', p.LineWidth, 'arcsec', 0, lpsy.pix2arcsec(pp.LineWidth));    
    lpsy.writeDvPmtr('--> Line width', lpsy.arcsec2pix(p.LineWidth), 'pix', 1, pp.LineWidth);       
    
    lpsy.writeDvPmtr('Line Color', mat2str(p.LineColor));
    
    lpsy.writeDvPmtr('Fixation Cross Size', p.fc_size, 'arcsec', 0, lpsy.pix2arcsec(pp.fc_size));
    lpsy.writeDvPmtr('Fixation Cross Line Width', p.fc_lw, 'arcsec', 0, lpsy.pix2arcsec(pp.fc_lw));
    lpsy.writeDvPmtr('Fixation Cross Duration', p.fc_dur, 's');
    lpsy.writeDvPmtr('Experiment Condition', p.cond, []);
    
    lpsy.writeDvPmtr('Pacman X Offset', p.x_pac_offset, 'arcsec', 0, lpsy.pix2arcsec(pp.x_pac_offset));
    lpsy.writeDvPmtr('Pacman Radius', p.pac_radius, 'arcsec', 0, lpsy.pix2arcsec(pp.pac_radius));
    lpsy.writeDvPmtr('Pacman Mouth', p.pac_mouth, 'arcsec', 0, lpsy.pix2arcsec(pp.pac_mouth));
    
    lpsy.writeDvPmtr('Flanker Offset', p.fl_offset, 'arcsec', 0, lpsy.pix2arcsec(pp.fl_offset));
    lpsy.writeDvPmtr('Number of Flankers', p.n_flankers, [], 1, pp.n_flankers);
    lpsy.writeDvPmtr('Flanker Length', p.flanker_len, 'arcsec', 0, lpsy.pix2arcsec(pp.flanker_len));
    
    lpsy.writeDvPmtr('Vernier Orientation', p.VernierOri, 'degrees');
    
    lpsy.writeDvPmtr('Peripheral Presentation', p.periph, []);
    lpsy.writeDvPmtr('Peripheral Offset', p.periph_offs, 'arcdeg', 0, lpsy.pix2arcdeg(pp.periph_offs));
    lpsy.writeDvPmtr('Peripheral Direction', p.periphDir, []);
    
    lpsy.writeDvPmtr('Background Luminance', p.BgLum, '%');


    %--- Reporting exact relative luminance values is tricky, because we use
    %    a non-identity hardware gamma table, which causes another rounding 
    %    operation. So we don't report an exact value at all rather than reporting 
    %    a value that is potentially as inaccurate as the specified value.
    lpsy.writeDvPmtr('Background luminance', p.BgLum, '%', 1, p.BgLum); 

    
    
    
    if dontRoundJitterX
        lpsy.writeDvPmtr('X position jitter (±, continuous)', p.PosJitter(1), 'pix', 1, pp.PosJitter(1));
    else
        lpsy.writeDvPmtr('X position jitter (±, whole pixels)', p.PosJitter(1), 'pix', 1, pp.PosJitter(1));
    end
    pj = p.PosJitter .* [1,1];
    lpsy.writeDvPmtr('Y position jitter (±, whole pixels)', pj(2), 'pix', 1, pp.PosJitter(2));
    lpsy.writeDvPmtr('Bug mode (0: no bug)', p.BugMode, []);
    lpsy.writeDvPmtr('Screen #', p.Screen, []);


end