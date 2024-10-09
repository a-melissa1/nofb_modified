function NoFB_noQuest_run

subjnr = input('Subject Code? ','s');
 
p.TestMode = 0;
 
%% Basic Setup Parameters
p.ScreenDist = 74;
 
%% Experiment Details
p.PracticeTrials = 10;%8
p.LineColor = [0,0,0];
 
% Vernier
p.VernierLength = 4000;%1260; % arcsec
p.VernierGap = 0;
p.LineWidth = 120; % arcsec
p.VernierDuration = 50;%200; %75; % 100 ms
p.VernierOri = 90;
 
%Flankers
p.periph = 1;
p.periph_offs = 5; % deg
 
 
p.x_pac_offset= 16000; % arcsecs
p.pac_radius= 1800; % arcsecs
p.pac_mouth= 2300; % arcsecs
p.fl_offset= 950; % arcsecs
p.n_flankers= 10; % number
p.flanker_len= p.VernierLength/4;
 
% Fixation Cross
p.fc_lw = 100;
p.fc_size=500;
 
%% Conditions
% 1 just vernier
% 2 2 flankers
% 3 8 flankers
% with mask m or without. With mask should stop feedback and stop
% processing.
 
%% Chaining the Conditions
conditions_cond = [1 2 3];
conditions_start_value = [100, 200, 200];
%ix = randperm(5);
%conditions_cond = conditions_cond(ix);
%conditions_start_value = conditions_start_value(ix);
 
conditions_ntrials = [60, 60, 60];
for m = [0, 1]
%% Call stimulus program.
for cc = [1:length(conditions_cond)]
cond = conditions_cond(cc);
p.Trials = conditions_ntrials(cc);
p.VernierOffset = conditions_start_value(cc);
p.mask = m;
p.cond = cond;
p.DvFilename = convertStringsToChars(fullfile('results',strcat('nofb_', subjnr, '_cond',num2str(m), '_',num2str(cond),'.dv')));
p_pvcell = struct2pvcell(p);
NoFB_noQuest(p_pvcell{:});
disp('done')
end
end
 
end
