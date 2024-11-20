function Exp = run_saccade_detection(Exp, varargin)
% RUN SACCADE DETECTION
% Detects saccades for each trial of an Exp struct using the algorithm
% described in Engbert and Mergenthaler (2006)
%
%-------------------------------------------------------------------
% Input:
%
% Exp [struct]: the marmoview struct
%               must contain fields vpx and D
%
% Output:
%
% Exp [struct]: modified marmoview struct contains slist
%   slist [array]: K x 7, where K is number of saccades
%       [startsac endsac peaksac kstartsac kendsac kpeaksac flag]]
%           column 1:  startsac - in time (secs), start of saccade
%           column 2:  endsac - in time (secs)
%           column 3:  peaksac - time of peak velocity
%           column 4:  kstartsac - time as integer of smoEye array
%           column 5:  kstartsac - time as integer of smoEye array
%           column 6:  kstartsac - time as integer of smoEye array
%           column 7:  flag - 0 by default, 4 - if a curved saccades
%                      and custom analysis can make other distinctions
%                      like flag the saccade that hits a target
%
% Optional Arguments:
%
% 'VFactor'     (5)     Relative Velocity Threshold (from Engbert and Trukenbrod, 2014)
% 'MinDuration' (20)    Minumum duration (in ms)
% 'MinGap'      (10)    Minimum gap between saccades before checking if
%                       they are a single saccade (in ms)
% 'FlagCurve'   (1.2)   Flag curved saccades (1 = straight)
% 'SampRate'    (1000)  Sampling rate of the eye trace

% see also: microsaccMerge, +saccadeflag.flag_saccades, +saccadeflag.eyeposition_smooth


ip = inputParser();
ip.KeepUnmatched = true;
ip.addParameter('ShowTrials', false)
ip.parse(varargin{:});

% % convert unmatched arguments for input to the saccade_flag function
% args = [fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]';

ShowTrials = ip.Results.ShowTrials;

if ShowTrials
    H = figure(1); clf
    set(H,'position',[100 200 800 800]);
end

%% convert eye position to degrees
nTrials = numel(Exp.D);
validTrials = 1:nTrials;

% gain and offsets from online calibration
cx = cellfun(@(x) x.c(1), Exp.D(validTrials));
cy = cellfun(@(x) x.c(2), Exp.D(validTrials));
dx = cellfun(@(x) x.dx, Exp.D(validTrials));
dy = cellfun(@(x) x.dy, Exp.D(validTrials));

% use the most common value across trials (we should've only calibrated
% once in these sessions)
cx = mode(cx);
cy = mode(cy);
dx = mode(dx);
dy = mode(dy);

Fs = 1./median(diff(Exp.vpx.raw(:,1)));
% x and y position
vxx = Exp.vpx.raw(:,2);
vyy = Exp.vpx.raw(:,3);

% convert to d.v.a.
vxxd = (vxx - cx)/(dx * Exp.S.pixPerDeg);
vyy = 1 - vyy;
vyyd = (vyy - cy)/(dy * Exp.S.pixPerDeg);

vtt = Exp.vpx.raw(:,1);
% Fs = median(diff(vtt));
vpp = Exp.vpx.raw(:,4);

% bad = ((hypot(diff(vxxd), diff(vyyd))./Fs) > 2e4);
% bad = find(filtfilt(ones(5,1), 1, double(bad))>0);
% vxxd(bad) = nan;
% vyyd(bad) = nan;

vxxd = repnan(vxxd, 'linear');
vyyd = repnan(vyyd, 'linear');


vxx = medfilt1(vxxd, 3);
vyy = medfilt1(vyyd, 3);
% 
vxx = imgaussfilt(vxx, 3);
vyy = imgaussfilt(vyy, 3);

vx = [0; diff(vxx)];
vy = [0; diff(vyy)];

% vx = sgolayfilt(vx, 1, 3);
% vy = sgolayfilt(vx, 1, 3);

% convert to d.v.a / sec

vx = vx * Fs;
vy = vy * Fs;

spd = hypot(vx, vy);
acc = [0; diff(spd)];

% slist = pdsa.detectSaccades(vtt', [vxx vyy]', 'detectThresh', 20, 'filterLength', 25, 'minISI', 10, 'verbose', false);

slist = +saccadeflag.flag_saccades([vtt vxx vyy vpp vx vy spd], ...
    'VelThresh', 10, ...
    'VFactor', 5, ...
    'MinDuration', 4, ...
    'MinGap', 2, ...
    'FlagCurve', 1.2, ...
    'SampRate', ceil(1/Fs));
size(slist,1)

vpx2ephys = synchtime.sync_vpx_to_ephys_clock(Exp);

% get proper times for the saccades
slist(:,1) = vtt(slist(:,4));
slist(:,2) = vtt(slist(:,5));
slist(:,3) = vtt(slist(:,6));

% remove saccades with unreasonable durations
sacdur = slist(:,2) - slist(:,1);
bad = sacdur < 0.005 | sacdur > 0.12;
slist(bad,:) = [];


% store values
Exp.slist = slist;
Exp.vpx2ephys = vpx2ephys;
Exp.vpx.smo = [vtt vxxd vyyd vpp vx vy spd];

% --- Saccade / Eye Trace QA
nTimePoints = numel(vxx);
Exp.vpx.Labels = zeros(nTimePoints, 1);
Exp.vpx.Labels(:) = 4; % initialize everything to lost track
Exp.vpx.LabelIds = {'Fixation', 'Saccade', 'Blink', 'Lost'};

% fixation is when the eye velocity is less than 8 d.v.a/s
Exp.vpx.Labels( Exp.vpx.smo(:,7) < 8 ) = 1;

% insert saccades
nSaccades = size(Exp.slist,1);
fprintf(1, 'Found %d saccades\n', nSaccades)


% expand blinks with some padding
blink = Exp.vpx.smo(:,4)==0; % when pupil size is zero
padding = 5;
bc = ones(padding,1);
blink = filtfilt(bc, 1, double(blink)) > 0;
Exp.vpx.Labels(blink) = 3;

% lost track is all other times!
% losttrack = isnan(Exp.vpx.smo(:,2));
% Exp.vpx.Labels(losttrack) = 4;

% flag velocities that make no sense
badVel = Exp.vpx.smo(:,7) > 2e3; % unphysiologically plausible  velocities 
Exp.vpx.Labels(badVel) = 4;

% tracker is unreliable more than 20 d.v.a from center (flag and label)
offScreen = hypot(Exp.vpx.smo(:,2),Exp.vpx.smo(:,2)) > 20;
Exp.vpx.Labels(offScreen) = 4;

% if any saccades intersect with blinks, remove them
bad = blink(Exp.slist(:,4)) | blink(Exp.slist(:,5));
Exp.slist(bad,:) = [];

outlier = isnan(Exp.vpx.raw(:,2));
Exp.vpx.Labels(outlier) = 4;
bad = outlier(Exp.slist(:,4)) | outlier(Exp.slist(:,5));
Exp.slist(bad,:) = [];
% 
% % if any saccades intersect with fixation, remove
% fixations = (Exp.vpx.Labels==1);
% bad = fixations(Exp.slist(:,4)) | fixations(Exp.slist(:,5));
% Exp.slist(bad,:) = [];
% 
% % if any fixations have absurd velocities, remove those points
% fixations = find(Exp.vpx.Labels==1);
% bad = Exp.vpx.smo(fixations,7) > 100; % if velocity greater than 100 deg/s
% Exp.vpx.Labels(fixations(bad)) = 4; % label as lost track

% remove saccades that have peak velocities that are unreasonable

% losttrack = Exp.vpx.Labels==4;

% if any peak velocities are at the time of a lost track
% badVel = find(losttrack(Exp.slist(:,6)));
% nBad = numel(badVel);

% n = zeros(nBad,1);
% for i = 1:nBad
%     s = Exp.slist(badVel(i),:);
%     spd = Exp.vpx.smo(s(4):s(5),7);
%     spd = diff(spd);
%     % Find downward and upward going zero-crossings
%     previous = spd(1:end-1);
%     current = spd(2:end);
%     down = find(previous > 0 & current < 0);
%     up = find(previous < 0 & current > 0);
%     nzc = numel([down; up]);
%     n(i)= nzc;
% end

% Exp.slist(badVel,:) = [];


% % remove ridiculous amplitudes
% % or unreasonable amplitudes
% x = Exp.vpx.smo(:,2);
% y = Exp.vpx.smo(:,3);
% 
% amp = hypot(x(Exp.slist(:,5)) - x(Exp.slist(:,4)), y(Exp.slist(:,5)) - y(Exp.slist(:,4)));
% v = Exp.vpx.smo(Exp.slist(:,6),7);
% bad = (amp > 20 | amp == 0) | (v > 200*amp);
% Exp.slist(bad,:) = [];
% 


% insert saccades
nSaccades = size(Exp.slist,1);
fprintf(1, 'Ended with %d saccades after QA\n', nSaccades)
for iSac = 1:nSaccades
    ix = Exp.slist(iSac,4):Exp.slist(iSac,5);
    Exp.vpx.Labels(ix) = 2;
end

%% plot it
if ShowTrials
    
   Labels = Exp.vpx.Labels;
    slist = Exp.slist;
    for iSac = 1:nSaccades
        ix = slist(iSac,4):slist(iSac,5);
        Labels(ix) = 2;
    end
    
    %%
    figure(1); clf
    for  i = 1:2
        subplot(2,1,1)
        ix = find(Labels==i);
        plot(ix, vxx(ix), '.'); hold on
        plot(ix, vyy(ix), '.');
        
        subplot(2,1,2)
        plot(ix, spd(ix), '.'); hold on
        
    end
    drawnow
    
    
    
    
    figure(2); clf
    xs = Exp.vpx.smo(Exp.slist(:,4),2);
    ys = Exp.vpx.smo(Exp.slist(:,4),3);
    xe = Exp.vpx.smo(Exp.slist(:,5),2);
    ye = Exp.vpx.smo(Exp.slist(:,5),3);

    dx = xe-xs;
    dy = ye-ys;
    dd =  hypot(dx, dy);

    dt = Exp.vpx.smo(Exp.slist(:,5),1) - Exp.vpx.smo(Exp.slist(:,4),1);
    v = Exp.vpx.smo(Exp.slist(:,6),7);
    
    subplot(1,2,1)
    plot(dd, v, '.')
    
    subplot(1,2,2)
    plot(dd, dt, '.')
    
    win = [-100 500];
    st = Exp.slist(:,4);
    en = Exp.slist(:,5);
    st = st(2:end);
    en = en(1:end-1);
    dur = st - en;
    [~, ind] = sort(dur);
    [an, ~, ~, wfs] = eventTriggeredAverage(spd, en(ind), win);
    figure, imagesc(wfs)
    
    figure, imagesc(wfs(1:100,:))
end

