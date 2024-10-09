function NoFB(varargin)
 
 
global sci p pp data
 
addpath("Elements");
 
%% --- Parse function parameters
p = inputParser;
p.addParamValue('TestMode', 0, @isscalar); % 0: Production.
% 1: Generate responses automatically w/o speeding up trial timing.
% 2: Interactively vary x stimulus position and vernier offset (vernier rendering test).
% 9: Run the experiment with the fewest shortcuts possible, while rushing through fast and quietly.
p.addParamValue('BugMode', 0, @isscalar); % 0: no bug. 1: 'BgLum' bug (see history for 25.01.2022).
p.addParamValue('Screen', [], @isscalar); % Only needed if using other than the 1st screen or other than specified in system.setup.
p.addParamValue('DvFilename',[], @(x) isempty(x) || ischar(x)); % Filename (possibly with path) for results. []=open file selection box; ''=don't save (for testing).
p.addParamValue('ScreenDist', 200, @(x) isscalar(x) && x>0); % [cm] Distance eye-to-monitor.
p.addParamValue('InfoScreenFct', @(wnd,txt,isTest) lpsy.showInfoScreen(wnd,txt,isTest));
% Handle of a info screen presentation function, pointing to: [isEsc] = ShowInfo(wnd, txt, isTest);
% See lpsy.showInfoScreen() for a prototype and more info.
p.addParamValue('NegFeedbackVol', 20, @isscalar); % [%] Volume for acoustical feedback for wrong response.
p.addParamValue('WarningsVol', 20, @isscalar); % [%] Volume of low one for warning beeps (too early response, no response).
 
p.addParamValue('TestValueRange', [], @(x) isempty(x) || numel(x)==2); % [same units as controlled parameter] Absolute test level limits. This range might be further limited when needed.
 
p.addParamValue('PracticeTrials', 8, @isvector) % Number of practice trials, or a vector of test level values, with one value per practice trial.
% If not a vector and if a paramter will be PEST-controlled, the first half of the practice trials will be presented at
% 1.2*startLevel, the second half at 0.8*startLevel.
% A given left/right stimulus alternative is never shown more than twice in a row. Trials with early or no response
% do not count and will be repeated.
p.addParamValue('StartTrialReps', [3,6], @isvector); % Number of repetitions of the first regular trial (min/max). At least "min" trials the start test have to be valid and answered
% correctly in a row before proceeding with PEST/QUEST, unless at most "max" valid trials have been presented. Only the last trial,
% whether correct or incorrect, is taken as a regular trial; the other trials will be reported as practice trials.
% For disabling the feature, pass [0,0], [0,1], or [1,1].
p.addParamValue('Trials', 80, @isscalar) % How many trials (excluding practice trials)?
p.addParamValue('VernierOffset', 75, @isscalar); % [arcsec] X-delta between upper and lower vernier segment. This is possibly just the PEST starting value, if PEST contrals the vernier offset.
p.addParamValue('VernierDuration', 20, @isscalar); % [ms]
p.addParamValue('ResponseTimeout', 30, @isscalar); % [s] Starting to count from stimulus offset (0:infinite).
p.addParamValue('ITI', [0.8 1.2], @(x) numel(x)>0 && numel(x)<=2); % [s] Basically the time between the response and the next stimulus onset.
% This can be an interval, in which case the ITI is randomly drawn from this interval after each trial).
% removed option Line Mode. Only aaL is used in this experiment
 
%% DISPLAY PARAMS
 
p.addParamValue('LineWidth', 43, @(x) isscalar(x) && x>=1); % [arcsec] In case of anti-aliasing, this refers to the virtual width of the line (which has to be 1 pixel minimum).
% In the protocol, this parameter is also reported in units of pixels. When using anti-aliasing, aim for a
% width of ??.5 pixels (i.e., 1.5, 2.5, 3.5,...).
p.addParamValue('VernierLength', 1260, @isscalar); % [arcsec] Total length (including gap).
p.addParamValue('VernierGap', 60, @isscalar); % [arcsec] Vertical gap. Only if this value is 0, there will be no gap. Otherwise, there will be at least a small gap.
p.addParamValue('PosJitter', 1.9, @(x) numel(x)>0 && numel(x)<=2); % [pix,�] �-Jitter of the whole stimulus for alleviating systematic imperfections of
% the screen surface or anti-aliasing artifacts. Normally, the jitter is subject to round-off in that the stimulus
% will be displaced an according integer multiple of pixels only, if at all. This only exception are displacements
% along X if one of the anti-aliasing modes is used and the jitter value is not whole-numbered.
% Pass a 2-value vector for specifying different values for X and Y.
p.addParamValue('BgLum', 50, @isscalar); % [%] Background Luminance
% To Do: Line color as Luminance?
% p.addParamValue('LineLum', 100, @isscalar); % [%]
p.addParamValue('LineColor', [200,100,100], @isvector); % [rgb] Line Color
p.addParamValue('col', [200,200,200], @isvector); % [rgb] Color of ???
 
p.addParamValue('cond', 1, @isscalar); % Condition Number
 
p.addParamValue('fc_size', 1000, @isscalar); % [arcsec] fixation cross line legth
p.addParamValue('fc_LineColor', [0,0,0], @isvector); % [rgb] fixation cross line color
p.addParamValue('fc_lw', 150, @isscalar); % [arcsec] fixation cross line width
p.addParamValue('fc_dur', 0.9, @isscalar); % [secs] fixation cross will be +- 0.5
 
p.addParamValue('x_pac_offset', 10000, @isscalar); % [arcsec] pacman x offset
p.addParamValue('pac_radius', 2000, @isscalar); % [arcsec] pacman radius
p.addParamValue('pac_mouth', 1000, @isscalar); % [arcsec] height of pacman mouth
p.addParamValue('fl_offset', 1000, @isscalar); % [arcsec] flanker offset
p.addParamValue('n_flankers', 10, @isscalar); % [#] number of flankers
p.addParamValue('flanker_len', 5000, @isscalar); % [arcsec] length of flankers
 
p.addParamValue('VernierOri', 0, @isscalar); % [degrees] Can only be 0 or 90.
 
p.addParamValue('periph', 0, @isscalar); % [0, 1] target in the center (0) or the periphery (1)
p.addParamValue('periph_offs', 8, @isscalar); % [arcdeg] taget offset
p.addParamValue('periphDir', 2, @isscalar); % [0, 1, 2] direction of the stimulus, alternating randomly (0), left (1) or right (2)
p.addParamValue('mask', 1, @isscalar); % [0, 1] mask on or off
 
%--- Assign all input parameters to the structure 'p'.
p.parse(varargin{:});
p = p.Results;
if p.TestMode==9
%--- Speed up things in TestMode 9.
p.NegFeedbackVol = 0;
p.WarningsVol = 0;
p.ITI = 0.05;
p.VernierDuration = 30;
end
assert(any(strcmpi(p.AdaptivePmtr,{'none','Vofs','Vdur'})), 'Invalid string for ''AdaptivePmtr'' parameter.');
assert(any(strcmpi(p.AdaptiveProc,{'pest','quest'})), 'Invalid string value for ''AdaptiveProc'' parameter.');
assert(p.VernierOri==0 || p.VernierOri==90, 'Invalid option for ''VernierOri'' parameter.');
%--- 'pp' might have been used as a global in another programs and
% will not be defined here than by adding fields, so clear any
% left-overs without removing it as a global variable.
pp = [];
 
%%
try
%% --- Open the screen.
 
lpsy.init();
Screen('Preference', 'SkipSyncTests', 1);
if isempty(p.Screen)
p.Screen = lpsy.getSysPmtr('Several','DefaultScreenNo',0);
end
pp.BugMode = p.BugMode;
pp.BgLum = min(255, round(p.BgLum/100*255));
wnd = lpsy.openWindow(p.Screen, pp.BgLum);
%--- For accurate anti-aliasing we have to rely on a good monitor
% calibration (i.e. gamma correction), and we need blending to be enabled.
% This is also used for text rendering.
% Both (gamma and blending) should not have an effect on the stimulus
% if anti-aliasing is actually not used, provided we don't make mistakes
% with rounding pixel coordinates and such.
 
lpsy.setGammaTab(1.0, 'gray');
Screen('BlendFunction', wnd, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
%---
sci = Screen('Resolution',p.Screen);
sci.wnd = wnd;
lpsy.setScreenDistCM(p.ScreenDist);
HideCursor();
 
%% --- Initialize pretty text rendering.
lpsy.prepDrawPrettyText(wnd);
 
%% --- Unit conversion and rounding of parameter values.
% pp-units are usually pixels, refresh frames, color values ex [0,255].
% Values all possibly rounded. Most pp names are the same as p names.
 
pp.VernierDuration = max(1, round(lpsy.msec2frame(p.VernierDuration)));
pp.TestValueRange = p.TestValueRange;
pp.LineColor = p.LineColor;
pp.fc_lc = p.fc_LineColor;
pp.fc_size = lpsy.arcsec2pix(p.fc_size);
pp.fc_lw = lpsy.arcsec2pix(p.fc_lw);
% luma = 0.299 * p.col(1) + 0.587 * p.col(2) + 0.114 * p.col(3);
 
pp.x_pac_offset = round(lpsy.arcsec2pix(p.x_pac_offset));
% pp.y_pac_offset is computed!
pp.pac_radius = lpsy.arcsec2pix(p.pac_radius);
%pp.pac_mouth = lpsy.arcsec2pix(p.pac_mouth);
pp.pac_mouth = 2*round(lpsy.arcsec2pix(p.pac_mouth)/2);
pp.fl_offset = lpsy.arcsec2pix(p.fl_offset);
pp.n_flankers = p.n_flankers;
pp.flanker_len = lpsy.arcsec2pix(p.flanker_len);
 
%--- For the sake of simplicity, we never use anti-aliasing along the vertical
% axis, and we always center the stimulus between screen pixels vertically.
% This requires the vernier gap to have an even number of pixels.
if p.VernierGap==0
pp.VernierGap = 0;
else
pp.VernierGap = max(2, 2*round(lpsy.arcsec2pix(p.VernierGap)/2));
end
pp.VernierLength = max(pp.VernierGap+2, 2*round(lpsy.arcsec2pix(p.VernierLength)/2));
pp.VernierOffset = lpsy.arcsec2pix(p.VernierOffset);
pp.LineWidth = max(1, lpsy.arcsec2pix(p.LineWidth));
 
pp.fc_lc = p.fc_LineColor;
pp.fc_size = lpsy.arcsec2pix(p.fc_size);
pp.fc_lw = max(1, lpsy.arcsec2pix(p.fc_lw));
pp.periph_offs = lpsy.arcdeg2pix(p.periph_offs);
 
 
%--- Possibly expand a single value to a 2-value vector and round up the value(s).
pp.PosJitter = ceil(p.PosJitter(:))' .* [1,1];
dontRoundJitterX = false;
%--- To which multiple number of pixels do we have to round the vernier offset?
% Note that the center of a 0-offset vernier defines the center of the stimulus,
% which might either lie on a pixel or between two pixels (along X), depending
% on the line width.
 
if rem(p.PosJitter(1),1)>0
%--- Anti-aliasing and a not whole-numbered jitter value for X: Use the value "as is".
pp.PosJitter(1) = p.PosJitter(1);
dontRoundJitterX = true;
end
 
%% --- More preparations.
if p.ResponseTimeout==0
pp.ResponseTimeout = 9999;
else
pp.ResponseTimeout = p.ResponseTimeout;
end
 
pChance = 0.5;

%% --- Prepare saving user parameters (more will follow).
% Doing this here helps detecting errors early on (useful during program development).
writeppp2DV(p, pp, myPEST, myQUEST, pest, practiceTestLevels, dontRoundJitterX)
%% --- Practice/regular loop.
if p.TestMode==2
InteractiveRenderTest();
end
%--- In case the window was open already, make sure the screen is cleared.
lpsy.flip(wnd);

dataCnt = 0;
minusPlusOne = [-1,1]; % -1: left: 1:right
vOffsetDir = minusPlusOne(randi(2));
sameDirCnt = 0;
phaseName = {'regular','practice'};
maxTrialCnt = [p.Trials, length(practiceTestLevels)];
%--- The repetitions of the start trial are controlled by 'startTrialCnt'.
% We only count valid trials, and the first element counts correct responses
% in a row, the second element counts all valid repetitions.
startTrialCnt = [0,0];
for isPracticePhase=[1,0]
if maxTrialCnt(isPracticePhase+1)==0
continue
end
isEsc = p.InfoScreenFct(wnd, phaseName{isPracticePhase+1}, any(p.TestMode==[1,9]));
if isEsc
error('Aborted with ESC...');
end
endOfTrialSecs = GetSecs();
trialCntValid = 0;
%--- Trial loop
while trialCntValid < maxTrialCnt(isPracticePhase+1)
%% --- Get next test level. We shouldn't need to check here for the valid

if isPracticePhase
dataPoolBase = dataPoolBase+4;
end
%--- Get left/right stimulus alternative while obeying the no-more-than-4-of-the-same rule.
if minusPlusOne(randi(2))~=vOffsetDir || sameDirCnt==4 || (isPracticePhase && sameDirCnt==2)
vOffsetDir = -vOffsetDir;
sameDirCnt = 1;
else
sameDirCnt = sameDirCnt + 1;
end
 
%% --- Stimulus presentation and response gathering.
% Break out of the while loop on escape key press.
while true
%--- Jitter the entire stimulus, usually restricted to an integer multiple of pixels.
% For getting, after rounding, the same probability for all integer values, especially for
% the min and max integer value, we extend the range for the uniform distribution by �(0.5-epsilon).
% Note that using randi() would be another option but makes it less elegant to handle both
% coordinates, X and Y, at once.
center = [sci.width, sci.height]/2 + round((2*rand()-1)*(pp.PosJitter+0.4999));
if dontRoundJitterX
%--- Exception...
center(1) = sci.width/2 + (2*rand()-1)*pp.PosJitter(1);
end
 
%--- Await ITI.
iti = p.ITI(:)' .* [1,1];
thisITI = iti(1) + diff(iti)*rand();
isEsc = lpsy.getKey({}, max(0, endOfTrialSecs + thisITI - GetSecs()));
if isEsc
break
end
 
% decide position of the stimulus
posx = center(1);
posy = center(2);
 
if p.periphDir == 0
side = rand() < 0.5;
elseif p.periphDir == 1
side = 0;
else
side=1;
end
 
if p.periph
if side
posx = center(1) + pp.periph_offs;
else
posx = center(1) - pp.periph_offs;
end
end
 
mask = randi(255, 300,300);
texture = Screen('MakeTexture', sci.wnd, mask);
rect = CenterRectOnPoint([0 0 100 100], posx , posy);
 
if p.VernierOri==0
keyList = KbName({'left','right'});
else
keyList = KbName({'up', 'down' });
end
%% BEGIN PRESENTATION
lpsy.releaseKeyWait();
 
%% 1. Present Fixation Cross
DrawFixationCross_basic(sci, pp, center(1), center(2))
%lpsy.petWatchdog()
lpsy.flip(sci.wnd);
% rand gives a number between 0 and 1. -0.5 gives a number
% between -0.5 and 0.5. fixation duration is random between
% those bounds,
fd = p.fc_dur+(rand()-0.5);
WaitSecs(fd);
 
 
%% 2. Present central vernier.
% prep Fixation Cross
if p.periph
DrawFixationCross_basic(sci, pp,center(1), center(2))
else
Screen('FillRect', sci.wnd, pp.BgLum, [], []);
end
% prep Vernier and Flankers
if p.VernierOri==0
vert=1;
else
vert=0;
end
DrawFlankers(sci, pp, posx, posy, p.cond, vert);
DrawVernier(sci, pp, posx, posy, vOffsetDir*pp.VernierOffset, p.VernierOri);
vernierOnsetSecs = lpsy.flip(sci.wnd);
flipLastSecs = vernierOnsetSecs;
flipDelay = lpsy.frame2sec(pp.VernierDuration);
 
%% 3. Mask or No mask post-phase
% prepare the mask
if p.mask
Screen('DrawTexture', sci.wnd, texture, [], rect);
DrawFixationCross_basic(sci, pp,center(1), center(2))
else
disp("here")
Screen('FillRect', sci.wnd, pp.BgLum, [], []);
DrawFixationCross_basic(sci, pp,center(1), center(2))
 
end
% flip the mask
[VBLTimestamp,isEsc,isTimingErr] = lpsy.flip(sci.wnd ,flipLastSecs, flipDelay);
%[~, isEsc, keyIdx, responseSecs] = lpsy.flipOrGetKey(sci.wnd, flipLastSecs, flipDelay, keyList,1);
 
%% 4. Wait time
%prep
DrawFixationCross_basic(sci, pp,center(1), center(2))
maskdur = 0.1;
maskOffsetSecs = lpsy.flip(sci.wnd, VBLTimestamp, maskdur);
%% Time to answer
Screen('FillRect', sci.wnd, pp.BgLum, [], []);
%DrawFixationCross_phat(sci, pp,center(1), center(2))
waitdur = 0.5;
startclick = lpsy.flip(sci.wnd, maskOffsetSecs, waitdur);
%--- First check for response.
if isEsc
break
end
 
%--- Await response.
if any(p.TestMode==[1,9])
%--- Fake response.
if any(dataCnt==[3,20])
%--- No response
keyIdx = 0;
responseSecs = GetSecs()+20;
elseif any(dataCnt==[5,22])
%--- Early response
keyIdx = 1;
responseSecs = GetSecs() + 100e-3;
else
switch lower(p.AdaptivePmtr)
case 'vdur'
x = lpsy.frame2msec(pp.VernierDuration);
thresh = 50;
otherwise
x = lpsy.pix2arcsec(pp.VernierOffset);
thresh = 40;
end
%--- Simulate a cumulative Gaussian response curve in log10-space.
thresh = log10(thresh);
x = log10(x);
sigma = 2.5/20; % 2.5dB (sigma in log10-space).
pCorrect = 0.5 + (1 - 0.5 - 0.02) * 0.5*erfc(-(x-thresh)/sigma/sqrt(2));
keyIdx = 1 + (vOffsetDir+1)/2;
if rand() > pCorrect
%--- Wrong response: "Flip" keyIdx.
keyIdx = 3 - keyIdx;
end
responseSecs = GetSecs() + 0.1 + 0.15*rand();
end
if p.TestMode~=9
pause(0.1);
end
isEsc = lpsy.getKey({},0);
else
 
[isEsc, keyIdx, responseSecs] = lpsy.getKey(keyList, p.ResponseTimeout);
lpsy.flip(sci.wnd)
WaitSecs(0.2)
end
 
break
end
%% --- Handle response
Screen('Close', texture);
if isEsc
error('Aborted with Escape.');
end
%---
if keyIdx
%--- We have a response.
answerTi = 1000*(responseSecs-vernierOnsetSecs);
%--- Check for early response.
isValid = answerTi > 150;
isEarly = ~isValid;
answerDir = minusPlusOne(keyIdx);
isHit = answerDir==vOffsetDir;
 
if isEarly && p.WarningsVol>0
Beeper(180, p.WarningsVol/100, 0.5);
Snd('Wait');
elseif (isEarly || ~isHit) && p.NegFeedbackVol>0
Beeper(2000, p.NegFeedbackVol/100, 0.03);
Snd('Wait');
end
else
%--- We have no response
isValid = 0;
isHit = 0;
answerTi = -1;
isEarly = 0;
answerDir = -1;
if p.WarningsVol>0
Beeper(180, p.WarningsVol/100, 0.5);
end
for i=1:3
diamondPts = bsxfun(@plus, [sci.width;sci.height]/2, pp.VernierLength*[1,0; 0,1; -1,0; 0,-1]');
FillPolyAA(wnd, [0,0,255], diamondPts, 1);
lpsy.flip(wnd);
WaitSecs(0.1);
lpsy.flip(wnd);
WaitSecs(0.1);
if any(p.TestMode==[1,9])
break
end
end
Snd('Wait');
end
lpsy.releaseKeyWait();
endOfTrialSecs = GetSecs();
 
 
%--- Proceed within the trial sequence.
trialCntValid = trialCntValid + double(isValid);
if ~isPracticePhase
%--- If we are past the start trial confirmation phase, the
% following changes to 'startTrialCnt' will not change that
% and will be undone.
startTrialCnt(2) = startTrialCnt(2) + double(isValid);
if isHit
startTrialCnt(1) = startTrialCnt(1) + 1;
else
startTrialCnt(1) = 0;
end
if any(startTrialCnt >= p.StartTrialReps)
%--- Start trial confirmation phase was or is over - keep it that way
% by setting startTrialCnt(2) to the maximum. This mechanism works
% even when p.StartTrialReps(2)==0.
startTrialCnt(2) = p.StartTrialReps(2);
else
%--- Start trial not confirmed: undo the 'trialCntValid' counting from above
% and put the trial into the practice trial pool.
trialCntValid = trialCntValid - double(isValid);
dataPoolBase = dataPoolBase+4;
end
end
%% --- Log the data.
dataCnt = dataCnt + 1;
dataPoolBase = dataPoolBase + 100;
data.POOL(dataCnt) = dataPoolBase + 1 + (vOffsetDir+1)/2;
data.LEVEL(dataCnt) = round(reportLevel);
data.VALID(dataCnt) = double(isValid);
data.HITS(dataCnt) = double(isHit);
data.REACT_TI(dataCnt) = round(answerTi);
data.isEarly(dataCnt) = double(isEarly);
data.isPractice(dataCnt) = double(isPracticePhase | startTrialCnt(2)<p.StartTrialReps(2));
data.vOfsDir(dataCnt) = vOffsetDir;
data.answerDir(dataCnt) = answerDir;
data.vOfs(dataCnt) = lpsy.pix2arcsec(pp.VernierOffset);
data.vDur(dataCnt) = lpsy.frame2msec(pp.VernierDuration);
data.iti(dataCnt) = round(1000*thisITI);
data.vOnsetTi(dataCnt) = vernierOnsetSecs;
data.stimSide(dataCnt) = side;
data.fc_dur(dataCnt) = fd;

end % trial loop
end % practice/regular
p.InfoScreenFct(wnd, 'done', p.TestMode~=0);
lpsy.cleanup();
%% --- Save the data.
if ~strcmp(p.DvFilename, '')
fid = lpsy.writeDvHeader(p.DvFilename);
if fid>=0
lpsy.writeDvPmtr(fid);
data.vOnsetTi = data.vOnsetTi - data.vOnsetTi(1);
lpsy.writeDvTable(fid, 0, data, 'vOfs',1, 'LEVEL', 0, 'vDur',1, 'vOnsetTi',3 );
fclose(fid);
else
save('results-emergency-save.mat');
end
end
catch err
lpsy.cleanup(err);
end
 
end % end of main function