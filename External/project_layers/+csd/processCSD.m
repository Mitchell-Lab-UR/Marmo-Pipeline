function CSD_structs = processCSD(Sessions, processed_data_dir, varargin )

% processCSD interacts with the +csd toolbox to obtain and plot CSDs
% Inputs:
%   Sessions                list of Exp and LFP session tags
%   processed_data_dir      Location of exp and lfp processed files. Exp
%                           files should be directly in processed_data_dir
%                           as {FileTag}.mat while LFP processed files should
%                           be in processed_data_dir/lfp as
%                           {FileTag}_lfp.mat
% Arguments: 
%   server_data_dir         string - path to server files for processing
%                           LFP
%   window                  [1X2] - time window to get CSD over 
%   noisetype               '3', '6', 'saccade3', 'saccade6' - Which CSD mode
%                           to use? 
%   spatial_smoothing       double - How much spatial smoothing (units:
%                           channels)
%   temporal_smoothing      double - How much temporal smoothing (units:
%                           time points)
%   subtract_saccade_csd    TRUE/FALSE - subtract out saccade CSD? 
%   subtract_pre_csd        TRUE/FALSE - subtract out baseline CSD across
%                           channels time-average from
%                           subtract_pre_CSD_time0 to
%                           subtract_pre_CSD_time1
%   get_reversal_from_average_CSD TRUE/FALSE - If there are multiple
%                           shanks, then average across the two shanks to get same reversal for
%                           both
%   save_plots              TRUE/FALSE - save plots of CSDs
%   plots_root_folder       Folder for saving plots if save_plots is true -
%                           default is 'plots' folder in processed_data_dir
%   process_lfp             if processed LFP doesn't exist, should we
%                           process it? - depracted! USE csd.processLFPs instead!
%   remove_dead_channels_kilosort TRUE/FALSE - remove dead channels identified from
%                           kilosort (missing channels from Exp.ycoords)
%   compute_nt3_reversal    TRUE/FALSE - if true, will also add noisetype 3
%                           reversal to CSD struct as stats.reversalPointDepth_nt3
%   replace_reversal_with_nt3 TRUE/FALSE - if true, will replace the CSDs
%                           reversal with the one from noisetype 3 (for
%                           example, if you want to use noisetype 3
%                           reversal to align for noisetype 6)
%   do_align_and_average    TRUE/FALSE - Get aligned and average CSDs across session? 
%   align_to                'top', 'bottom', 'nt3' - what reversal to align to if do_align_and_average=true
%                           top = align to minimum depth reversal; 
%                           bottom = align to maximum depth reversal; nt3 = align to 
%                           nt3 = align to the noisetype 3 reversal point
%                           (note: compute_nt3_reversal mus tbe true if
%                           using nt3)
%   noisetype6_rvsl_window  window to look to get noisetype 6 reversal
%   noisetype3_rvsl_window  window to look to get noisetype 3 reversal
%   saccade_rvsl_window     window to look to get noisetype saccade reversal
%
%   rvsl_manual_FileTags    FileTags to manually adjust reversal depth 
%   rvsl_manual_location    'top' or 'bottom' - specify top of bottom
%                           reversal to manually adjust
%                           where each entry corresponds to each entry in rvsl_manual_FileTags
%                           (top = min reversal if there are multiple; bottom = max reversal)
%   rvsl_manual_shank_index shank index to manually adjust
%                           reversal depth where each entry corresponds to each rvsl_manual_FileTags
%   rvsl_manual_shank_depth reversal depth value to manually adjust
%                           reversal depth where each entry corresponds to each rvsl_manual_FileTags
% Returns:
%   CSD_structs: all CSD structs
%

ip = inputParser();
ip.addParameter('window', [-100 200])

ip.addParameter('noisetype', '3')
ip.addParameter('spatial_smoothing', 1.0) % in units of channels
ip.addParameter('temporal_smoothing', 0.0) % in units of time points
ip.addParameter('subtract_saccade_csd', false)
ip.addParameter('subtract_pre_csd', false)
ip.addParameter('subtract_pre_CSD_time0', -100)
ip.addParameter('subtract_pre_CSD_time1', 0)


ip.addParameter('get_reversal_from_average_CSD', false)
ip.addParameter('plots_root_folder', NaN)
ip.addParameter('process_lfp', true)
ip.addParameter('remove_dead_channels_kilosort', true)
ip.addParameter('save_plots', false)
% ip.addParameter('exclude', true)

ip.addParameter('compute_nt3_reversal', false)
ip.addParameter('replace_reversal_with_nt3', false)
ip.addParameter('server_data_dir', '')
% ip.addParameter('processed_data_dir', '')

ip.addParameter('do_align_and_average', false)
ip.addParameter('align_to', 'nt3')

% Noisetype 6 parameters
ip.addParameter('noisetype6_rvsl_window', [50 70])

% noisetype3
ip.addParameter('noisetype3_rvsl_window', [40 50])
ip.addParameter('get_nt3_rvsl_from_two_smooth', true)

% saccade-locked
ip.addParameter('saccade_rvsl_window', [50 60])

% manually specify reversal
ip.addParameter('rvsl_manual_FileTags', {})
ip.addParameter('rvsl_manual_location', {}) % 'top' or 'bottom'
ip.addParameter('rvsl_manual_shank_index', [])
ip.addParameter('rvsl_manual_shank_depth', [])

ip.parse(varargin{:});


%%%%%% PARAMS %%%%%%%
% pick_session = ip.Results.pick_session;
do_align_and_average = ip.Results.do_align_and_average;
align_to = ip.Results.align_to; % top (nt6), bottom (nt6), nt3
save_plots = ip.Results.save_plots;
window = ip.Results.window; % CSD window
noisetype = string(ip.Results.noisetype);
spatial_smoothing=ip.Results.spatial_smoothing; 
temporal_smoothing = ip.Results.temporal_smoothing;
subtract_saccade_csd = ip.Results.subtract_saccade_csd;
subtract_pre_csd = ip.Results.subtract_pre_csd;
get_reversal_from_average_CSD = ip.Results.get_reversal_from_average_CSD; % get CSD reversal by averaging the CSD on the shanks and finding reversal instead of individually
noisetype6_rvsl_window = ip.Results.noisetype6_rvsl_window;
plots_root_folder = ip.Results.plots_root_folder;
% process_lfp = ip.Results.process_lfp;
remove_dead_channels_kilosort = ip.Results.remove_dead_channels_kilosort;
replace_reversal_with_nt3 = ip.Results.replace_reversal_with_nt3;
compute_nt3_reversal = ip.Results.compute_nt3_reversal;
% SERVER_DATA_DIR = ip.Results.server_data_dir;
% PROCESSED_DATA_DIR = ip.Results.processed_data_dir;
noisetype3_rvsl_window = ip.Results.noisetype3_rvsl_window;
get_nt3_rvsl_from_two_smooth = ip.Results.get_nt3_rvsl_from_two_smooth;
rvsl_manual_FileTags = ip.Results.rvsl_manual_FileTags;
rvsl_manual_location = ip.Results.rvsl_manual_location;
rvsl_manual_shank_index = ip.Results.rvsl_manual_shank_index;
rvsl_manual_shank_depth = ip.Results.rvsl_manual_shank_depth;
subtract_pre_CSD_time0 = ip.Results.subtract_pre_CSD_time0;
subtract_pre_CSD_time1 = ip.Results.subtract_pre_CSD_time1;
%%%%%%%%%%%%%%%%%%%%%%

if strcmp(align_to, 'nt3')
    compute_nt3_reversal = true;
end

% make sure manual reversal inputs are correct 
manual_rvsl_check2 = false;
manual_rvsl_check1 = sum(isempty(rvsl_manual_FileTags)) + sum(isempty(rvsl_manual_shank_depth)) + sum(isempty(rvsl_manual_shank_index));
if manual_rvsl_check1 == 0
    manual_rvsl_check2 = length(rvsl_manual_FileTags)==length(rvsl_manual_shank_depth) && length(rvsl_manual_FileTags)==length(rvsl_manual_shank_index) && length(rvsl_manual_shank_depth)==length(rvsl_manual_shank_index);
    manual_rvsl_check2 = ~manual_rvsl_check2;
end
manual_rvsl_check1 = manual_rvsl_check1 == 0 || manual_rvsl_check1 == 3;
manual_rvsl_check1 = ~manual_rvsl_check1;
if manual_rvsl_check1 || manual_rvsl_check2
    error('wrong input for rvsl_manual_FileTags, rvsl_manual_location & rvsl_manual_shank_index')
end

% if save_plots && any(isnan(plots_root_folder))
%     error("Please supply root directory (plots_root_dir) to save plots")
% end

% default
if any(isnan(plots_root_folder)) && save_plots
    plots_root_folder = fullfile(processed_data_dir, 'plots');
end
if ~exist(plots_root_folder, 'dir')
    mkdir(plots_root_folder)
end
% user = 'bluethunder';
% 
% switch user
%     case 'gabehome'
%          = '/Users/gabrielsarch/Documents/MitchellLabRepo/Data';
%          = '/Users/gabrielsarch/Documents/MitchellLabRepo/Data'; %C:\Processed';
%     case 'judehome'
%         SERVER_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Raw';
%         PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Processed';
%     case 'bluethunder'
% %         SERVER_DATA_DIR = 'C:\PSA_Gravedigger\Raw';
%         SERVER_DATA_DIR = 'C:\Raw';
% %         PROCESSED_DATA_DIR = 'C:\PSA_Gravedigger\Processed';
%         PROCESSED_DATA_DIR = 'Z:\Data\Processed_Laminar'; %C:\Processed';
%     case 'jakegravedigger'
%         SERVER_DATA_DIR = 'C:\Raw';
%         PROCESSED_DATA_DIR = 'Z:\PSA_EDF\EDF_Processed_Decoding\X-Y Recordings';
% end
% 
% 
% 
% if pick_session
%     DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
% else 
%     %Sessions = {'Milo_2021-02-15_12-41-31_MT32_1b', 'Milo_2021-02-22_11-26-14_MT32_4b', 'Milo_2021-03-17_11-57-28_MT64_5b', 'Milo_2021-03-18_11-05-37_MT64_6b', 'Milo_2021-04-07_12-40-34_MT64_7b', 'Milo_2021-04-12_11-54-37_MT64_9c', 'Milo_2021-04-16_12-15-05_MT64_11b', 'Milo_2021-04-21_11-51-59_MT64_14b', 'Milo_2021-05-10_11-49-28_MT64_20b', 'Milo_2021-05-14_11-02-57_MT64_22b', 'Milo_2021-07-08_11-06-38_MT64_30b', 'Milo_2021-07-13_11-45-41_MT64_31b', 'Milo_2021-07-15_11-24-53_MT64_32b'};
%     Sessions = {'Milo_2021-02-22_11-26-14_MT32_4b', 'Milo_2021-03-17_11-57-28_MT64_5b', 'Milo_2021-03-18_11-05-37_MT64_6b', 'Milo_2021-04-07_12-40-34_MT64_7b', 'Milo_2021-04-12_11-54-37_MT64_9c', 'Milo_2021-04-16_12-15-05_MT64_11b', 'Milo_2021-04-21_11-51-59_MT64_14b', 'Milo_2021-05-10_11-49-28_MT64_20b', 'Milo_2021-05-14_11-02-57_MT64_22b', 'Milo_2021-07-08_11-06-38_MT64_30b', 'Milo_2021-07-13_11-45-41_MT64_31b', 'Milo_2021-07-15_11-24-53_MT64_32b'};
%     %     Sessions = {'Milo_2021-03-17_11-57-28_MT64_5b', 'Milo_2021-03-18_11-05-37_MT64_6b', 'Milo_2021-04-07_12-40-34_MT64_7b', 'Milo_2021-04-16_12-15-05_MT64_11b', 'Milo_2021-04-21_11-51-59_MT64_14b', 'Milo_2021-05-10_11-49-28_MT64_20b', 'Milo_2021-05-26_13-01-25_MT64_26b', 'Milo_2021-06-08_11-19-32_MT64_29b'};
%     Sessions = [Sessions {'Milo_280721','Milo_300721','Milo_170821','Milo_240821','Milo_260821', 'Milo_160921'}];
% %     sess_num = 1;
% %     sess = Sessions{sess_num};
% %     disp(['Session: ' sess])
%     DataFolder = sess;
% end
% 
% lower_rvsl_3 = [];
% upper_rvsl_6 = [];
% lower_rvsl_6 = [];

CSD_structs = [];


CSDs = [];
avgCSD = [];
for k = 1:length(Sessions) %[3 4 5] %1:length(Sessions) % %
    sess_num = k;
    sess = Sessions{sess_num};
    
%     disp(['Session: ' sess])
%     if length(sess)==11
%         DataFolder = sess;
%         [~, FileTag] = fileparts(DataFolder);
%         sess = FileTag;
%     else
%         
%     end
    FileTag = sess;
    % DataFolder = fullfile(SERVER_DATA_DIR,FileTag);


%     if ~exist(char(strcat(PROCESSED_DATA_DIR,filesep,FileTag, '.mat')), 'file')
%         tag_split = split(FileTag, '_');
%         subj_name = tag_split(1);
%         date = split(tag_split(2), '-');
%         yr = split(date(1), '0');
%         yr = yr(2);
%         day = date(3);
%         month = date(2);
%         FileTag = strcat(subj_name, '_', day, month, yr);
%         FileTag = FileTag{1};
% 
%     end

%     if FileTag=='Milo_150221'
%         FileTag = [FileTag 'b'];
%     end
    disp(['FileTag: ', FileTag])
    
%     lfpFile = fullfile(PROCESSED_DATA_DIR,'lfp','Gabe_to_download',[FileTag '_lfp.mat']);
%     if ~exist(lfpFile, 'file')
%         disp('NOT IN DOWNLOADS')
%     end
%     continue

    if ~exist(char(strcat(processed_data_dir,filesep,FileTag, '.mat')), 'file')
        error('experiment file not found')
    end

    disp('done')

    % 3 Load EXP file
    disp('loading exp file')
    ExpFile = [processed_data_dir,filesep,FileTag,'.mat'];
    load(ExpFile);
    disp('done loading exp file!')
    
    disp('loading lfp...')
    lfpFile = fullfile(processed_data_dir,'lfp',[FileTag '_lfp.mat']);

%     if ~exist(lfpFile, 'file') && process_lfp
%         print('LFP processed file not found at given file location.. generating it (this will take a while)..')
%         [data, timestamps, info] = io.getLFP(ops, true, false);
%         unique_x = unique(Exp.osp.xcoords);
%         num_shanks = size(unique_x, 1);
%         shank_len = size(data, 2)/num_shanks;
%         deadChan = []; % Add dead channels here if there are any (1-32 first shank, 33-64 second shank, etc.)
%         if num_shanks > 0
% 
%             xcoords = ones(shank_len,num_shanks);
%             for i = 1:num_shanks
%                 xcoords(:,i) = xcoords(:,i)*unique_x(i);
%             end
%         end
%         diff_ycoords = abs(mode(diff(Exp.osp.ycoords)));
%         ycoords = repmat(flip(linspace(diff_ycoords, diff_ycoords*shank_len, shank_len)), num_shanks, 1)';
% 
%         %lfpFile = fullfile(PROCESSED_DATA_DIR,'lfp',[FileTag '_lfp.mat']);
%         % lfpFile = [PROCESSED_DATA_DIR,filesep,'lfp',filesep,FileTag,'_lfp.mat'];
%         disp('Saving LFP struct');
%         save(lfpFile, '-v7.3', 'timestamps', 'info', 'data', 'xcoords', 'ycoords', 'deadChan')
%         fprintf('LFP struct saved to %s\n',lfpFile);
%         fprintf('Mat File %s\n',[FileTag,'_lfp']);
%         disp('Saving LFP done');
%     end
    
    if ~exist(lfpFile, 'file')
        error('LFP file not found... make sure paths are correct or generate LFPs with csd.processLFPs.m')
    end

    % 5 Import LFP (once saved)
    lfp = load(lfpFile);
    disp('done loading LFP!')
    
    if remove_dead_channels_kilosort
        disp('Removing dead channels identified by kilosort...')
        lfp = csd.getDeadChannels(lfp, Exp);
    else
        lfp.deadChan = [];
    end
    
    
    
    if replace_reversal_with_nt3 || compute_nt3_reversal || noisetype=='3'
        if get_nt3_rvsl_from_two_smooth
            stats_sm1 = csd.getCSD(lfp, Exp, ...
                'plotIt', false, 'noisetype', 3, ...
                'spatsmooth', 1.0, 'tempsmooth', 1.0, ...
                'rvsl_window', [40 50], 'do_spectrogram', false, ... 
                'subtract_pre_CSD', subtract_pre_csd, 'subtract_pre_CSD_time0', subtract_pre_CSD_time0, 'subtract_pre_CSD_time1', subtract_pre_CSD_time1); 

            stats_sm3 = csd.getCSD(lfp, Exp, ...
                'plotIt', false, 'noisetype', 3, ...
                'spatsmooth', 3.0, 'tempsmooth', 7.0, ...
                'rvsl_window', [50 60], 'sink_first', false, ...
                'subtract_pre_CSD', subtract_pre_csd, 'subtract_pre_CSD_time0', subtract_pre_CSD_time0, 'subtract_pre_CSD_time1', subtract_pre_CSD_time1); 

            stats_rvsl = csd.identifyMTSink(stats_sm1, stats_sm3);
        else
            stats_rvsl = csd.getCSD(lfp, Exp, ...
            'plotIt', false, 'noisetype', 3, ...
            'spatsmooth', spatial_smoothing, 'tempsmooth', temporal_smoothing, ...
            'rvsl_window', noisetype3_rvsl_window, 'do_spectrogram', false, 'window', window, ... 
            'subtract_pre_CSD', subtract_pre_csd, 'subtract_pre_CSD_time0', subtract_pre_CSD_time0, 'subtract_pre_CSD_time1', subtract_pre_CSD_time1);
%             stats2.reversalPointDepth = stats_rvsl.reversalPointDepth;
        end
            
    end

    % 6 PLOT CSD - noisetype 3 (full field flash)
    % stats_6 = csd.getCSD(lfp, Exp, 'plotIt', false, 'noisetype', 6, 'exclude', true);
    % stats_sm3 = csd.getCSD(lfp, Exp, 'plotIt', false, 'noisetype', 3,'spatsmooth', 3.0, 'rvsl_window', [40 50]); 
    %input 'check'
    if strcmp(noisetype, 'saccade3') || strcmp(noisetype,'saccade6')
        rvsl_window = [50 60];
        stats2 = csd.getCSD(lfp, Exp, ...
        'plotIt', false, 'noisetype', noisetype, ...
        'spatsmooth', spatial_smoothing, 'tempsmooth', temporal_smoothing, ...
        'rvsl_window', rvsl_window, 'sink_first', false, 'window', window, ...
        'subtract_pre_CSD', subtract_pre_csd, 'subtract_pre_CSD_time0', subtract_pre_CSD_time0, 'subtract_pre_CSD_time1', subtract_pre_CSD_time1);
        %stats2.reversalPointDepth = stats_rvsl.reversalPointDepth;
        %%removed by AB on 10/17/24 due to error
        stats = stats2;
    elseif strcmp(noisetype,'3') || strcmp(noisetype,'6')
        if strcmp(noisetype,'3')
            stats = stats_rvsl;
            stats.reversalPointDepth_nt3 = stats_rvsl.reversalPointDepth;
        elseif strcmp(noisetype,'6')
            stats2 = csd.getCSD(lfp, Exp, ...
                'plotIt', false, 'noisetype', 6, ...
                'spatsmooth', spatial_smoothing, 'tempsmooth', temporal_smoothing, ...
                'rvsl_window', noisetype6_rvsl_window, 'do_spectrogram', false, 'window', window, ...
                'subtract_pre_CSD', subtract_pre_csd, 'subtract_pre_CSD_time0', subtract_pre_CSD_time0, 'subtract_pre_CSD_time1', subtract_pre_CSD_time1);
    %         stats_lower = csd.getCSD(lfp, Exp, ...
    %             'plotIt', false, 'noisetype', 6, ...
    %             'spatsmooth', spatial_smoothing, 'tempsmooth', 1.0, 'sink_first', false, ...
    %             'rvsl_window', [50 150], 'do_spectrogram', false, 'window', window);
            %figure
            %csd.plotCSD(stats2, 'overlayLFP', false)
    %         disp('plotting')
            if compute_nt3_reversal
                stats2.reversalPointDepth_nt3 = stats_rvsl.reversalPointDepth;
            end
            if replace_reversal_with_nt3
                stats2.reversalPointDepth = stats_rvsl.reversalPointDepth;
            end
            stats = stats2;
            
            
%             if FileTag=='Milo_080721'
%                 disp('NOTE: altering top reversal for Milo_080721')
%                 [~,I] = min(stats.reversalPointDepth{1});
%                 stats.reversalPointDepth{1}(I) = 350;
%             elseif FileTag=='Milo_120421'
%                 disp('NOTE: altering top reversal for Milo_120421')
%                 [~,I] = min(stats.reversalPointDepth{2});
%                 stats.reversalPointDepth{2}(I) = 525;
%             end
            
        end
    else
        error('dont know this noisetype')
    end
    
    if get_reversal_from_average_CSD
        CSD = stats.CSD;
        if size(stats.CSD, 3) > 1
            CSD = mean(CSD, 3);
            time = stats.time;
            depth = stats.depth;
            %if noisetype=='6'
            ix = time > noisetype6_rvsl_window(1) & time < noisetype6_rvsl_window(2);
            % here we want source
            [~,id] = max(reshape(-CSD(:,ix), [], 1));
    %                 elseif ip.Results.sink_first
    %                     time1 = rvsl_window(1);
    %                     time2 = rvsl_window(2);
    %                     ix = time > time1 & time < time2;
    %                     % sink should be the minimum value
    %                     [~,id] = max(reshape(-CSD(:,ix), [], 1));
    %                 else
    %                     time1 = rvsl_window(1);
    %                     time2 = rvsl_window(2);
    %                     ix = time > time1 & time < time2;
    %                     % sink should be the minimum value
    %                     [~,id] = max(reshape(CSD(:,ix), [], 1));
    %                 end

            % convert to indices
            [depthIndex,timeIndex] = ind2sub(size(CSD(:,ix)), id);

            time_ = time(ix);
            sinkTime = time_(timeIndex);
            % find reversal point
            CSD_ = CSD(:,ix);
            reversalPoints = findZeroCrossings(CSD_(:,timeIndex), 0);

            % Take first reversal after max
            %if noisetype=='6'
            reversalPoints1 = reversalPoints(reversalPoints<depthIndex);
            reversalPoints1 = max(reversalPoints1);
            reversalPoints2 = reversalPoints(reversalPoints>depthIndex);
            reversalPoints2 = min(reversalPoints2);
            reversalPoints = [reversalPoints1 reversalPoints2];
    %                     elseif ip.Results.sink_first
    %                         reversalPoints = reversalPoints(reversalPoints<depthIndex);
    %                         reversalPoints = max(reversalPoints);
    %                     else
    %                         reversalPoints = reversalPoints(reversalPoints<depthIndex);
    %                         reversalPoints = max(reversalPoints);
    %                     end

            depthIndex = depthIndex+1;
    %         figure()
    %         imagesc(CSD_)
            for shankInd=1:size(stats.CSD, 3)
                if isempty(reversalPoints)
                    stats.reversalPointDepth{shankInd} = NaN;
                    stats.sinkDepth{shankInd} = NaN;
                    stats.sinkTime{shankInd} = NaN;
                else
                    stats.reversalPointDepth{shankInd} = depth(reversalPoints); %reversalPoints(1);
                    stats.sinkDepth{shankInd} = depth(depthIndex);
                    stats.sinkTime{shankInd} = sinkTime;
                end
            end
        end
    end
    
    
    % replace with manually specified reversals
    if ~isempty(rvsl_manual_FileTags)
        for i=1:length(rvsl_manual_FileTags)
            filetag_manual = rvsl_manual_FileTags{i};
            if strcmp(FileTag,filetag_manual)
                fprintf('NOTE: Manually changing reversal point for file: %s\n',filetag_manual)
                shankind = rvsl_manual_shank_index(i);
                new_depth = rvsl_manual_shank_depth(i);
                if strcmp(rvsl_manual_location{i},'top') || any(isnan(rvsl_manual_location))
                    [~,I] = min(stats.reversalPointDepth{shankind});
                elseif strcmp(rvsl_manual_location{i},'bottom')
                    [~,I] = max(stats.reversalPointDepth{shankind});
                else
                    error('wrong input for rvsl_manual_location')
                end
                stats.reversalPointDepth{shankind}(I) = new_depth;
            end
        end
    end


    if subtract_saccade_csd
        stats_saccade = csd.getCSD(lfp, Exp, ...
                'plotIt', false, 'noisetype', saccade_noisetype, ...
                'spatsmooth', spatial_smoothing, 'tempsmooth', temporal_smoothing, ...
                'rvsl_window', [40 50], 'do_spectrogram', false, 'window', window, ...
                'subtract_pre_CSD', subtract_pre_csd, 'subtract_pre_CSD_time0', subtract_pre_CSD_time0, 'subtract_pre_CSD_time1', subtract_pre_CSD_time1);
        stats_saccade.reversalPointDepth = stats_rvsl.reversalPointDepth;
        if (1) % saccade avg
            saccade_csd = stats_saccade.CSD;
            time_avg_saccade_csd = mean(saccade_csd, 2);
        else % stimulus-locked avg -100 to 0 
            where_0 = find(stats.time==0);
            where_minus100 = find(stats.time==-100);
            if isempty(where_0) || isempty(where_minus100)
                error('Could not find time range')
            end
            csdonehund_csd = stats.CSD(:,where_0:where_minus100, :);
            time_avg_saccade_csd = mean(csdonehund_csd, 2);
        end

        avg_adjusted_csd = stats.CSD - time_avg_saccade_csd;
        stats_ = stats;
        stats.CSD = avg_adjusted_csd;
        %disp('here')
        if (0)
            % This plots (1) saccade-locked, (2) stimulus locked, (3) saccade
            % adjusted, (4) saccade time vector
            clims = [-2 2];
            figure(1); clf()
            subplot(1,4,1)
            imagesc(stats.time, stats.depth,stats_saccade.CSD(:,:,1),clims); axis ij
            title(sprintf('Saccade-locked (Noisetype %s) CSD', saccade_noisetype))
            subplot(1,4,2)
            imagesc(stats.time, stats.depth,stats_.CSD(:,:,1),clims); axis ij
            title(sprintf('Stimulus-locked (noisetype %s) CSD', string(noisetype)))
            subplot(1,4,3)
            imagesc(stats.time, stats.depth,stats.CSD(:,:,1),clims); axis ij
            title(sprintf('Saccade-adjusted stimulus-locked (noisetype %s) CSD', string(noisetype)))
            subplot(1,4,4)
            imagesc(0, stats.depth, time_avg_saccade_csd(:,:,1),clims); axis ij
            title('Time-averaged saccade adjustment vector')
            colorbar
            name = strcat(FileTag,'_noistype=',string(noisetype),'_spatial_smoothing=',string(spatial_smoothing));
            if subtract_saccade_csd
                name = strcat(name,'_saccade_time_avg_subtract'); 
            end
            name = strcat(name,'_FOUR_PLOTS'); 
            %title(name,'Interpreter','none')
            %xlabel('Time (ms)')
            %ylabel('Depth')
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gca, fullfile('C:\Users\Jake\Dropbox\MarmoLabWebsite\PSA\Code\Import', 'plots', 'saccade_test', char(name)), 'png')
        end
%     elseif (0) %subtract_pre_csd
%     %     stats_saccade = csd.getCSD(lfp, Exp, ...
%     %             'plotIt', false, 'noisetype', saccade_noisetype, ...
%     %             'spatsmooth', spatial_smoothing, 'tempsmooth', 1.0, ...
%     %             'rvsl_window', [40 50], 'do_spectrogram', false, 'window', window);
%     %     stats_saccade.reversalPointDepth = stats_rvsl.reversalPointDepth;
%         where_0 = find(stats.time==0);
%         where_minus100 = find(stats.time==-100);
%         if isempty(where_0) || isempty(where_minus100)
%             error('Could not find time range')
%         end
%         csdonehund_csd = stats.CSD(:,where_0:where_minus100, :);
%         time_avg_saccade_csd = mean(csdonehund_csd, 2);
% 
%         avg_adjusted_csd = stats.CSD - time_avg_saccade_csd;
%         stats_ = stats;
%         stats.CSD = avg_adjusted_csd;
%         %disp('here')
%         if (1)
%             % This plots (1) saccade-locked, (2) stimulus locked, (3) saccade
%             % adjusted, (4) saccade time vector
%             clims = [-2 2];
%             figure(1); clf()
%     %         subplot(1,4,1)
%     %         imagesc(stats.time, stats.depth,stats_saccade.CSD(:,:,1),clims); axis ij
%     %         title(sprintf('Saccade-locked (Noisetype %s) CSD', saccade_noisetype))
%             subplot(1,3,1)
%             imagesc(stats.time, stats.depth,stats_.CSD(:,:,1),clims); axis ij
%             title(sprintf('Stimulus-locked (noisetype %s) CSD', string(noisetype)))
%             subplot(1,3,2)
%             imagesc(stats.time, stats.depth,stats.CSD(:,:,1),clims); axis ij
%             title(sprintf('Saccade-adjusted stimulus-locked (noisetype %s) CSD', string(noisetype)))
%             subplot(1,3,3)
%             imagesc(0, stats.depth, time_avg_saccade_csd(:,:,1),clims); axis ij
%             title('Time-averaged saccade adjustment vector')
%             colorbar
%             name = strcat(FileTag,'_noistype=',string(noisetype),'_spatial_smoothing=',string(spatial_smoothing));
%             if subtract_saccade_csd
%                 name = strcat(name,'_preCSD_time_avg_subtract'); 
%             end
%             name = strcat(name,'_THREE_PLOTS'); 
%             %title(name,'Interpreter','none')
%             %xlabel('Time (ms)')
%             %ylabel('Depth')
%             set(gcf, 'Position', get(0, 'Screensize'));
%             saveas(gca, fullfile('C:\Users\Jake\Dropbox\MarmoLabWebsite\PSA\Code\Import', 'plots', 'saccade_test', char(name)), 'png')
%         end
    end
    
    
    CSD_structs = [CSD_structs; stats];




    % figure(1); clf
    % if (0)
    %     gamma = csd.getGamma(lfp, 'method', 'weightedMin');
    %     csd.plotCSD(stats, 'overlayLFP', false, 'gamma', gamma)
    % else
    %     csd.plotCSD(stats, 'overlayLFP', false)
    % end


    % % save event-triggered lfp trace from MT 
    % figure(1); clf()
    % shankInd = 1;
    % data_lfp = bsxfun(@plus, stats.STA(:,:,shankInd), stats.chDepths);
    % time = stats.time;
    % time_keep = time>=0 & time<=80;
    % data_lfp = data_lfp(:,time_keep);
    % time = time(time_keep);
    % plot(time, data_lfp, 'Color', repmat(.1, 1, 3)); axis ij
    % xlabel('time (ms)', 'FontSize', 30)
    % ylabel('Depth (microns)', 'FontSize', 30)
    % ax=gca;
    % ax.FontSize = 30;
    % set(gca,'TickLength',[0 .01])
    % set(gca,'box','off');
    % xlim([0 80])
    % save_folder = 'C:\Users\Jake\Dropbox\MarmoLabWebsite\PSA\Code\Import\plots\poster_images';
    % save_name = 'LFP_CSD_MT.mat';
    % split_name = split(save_name, '.');
    % save_name = [split_name{1} '.png'];
    % set(gcf, 'Position', get(0, 'Screensize'));
    % saveas(gca, fullfile(save_folder, save_name));
    % close all


    numShanks = size(stats.CSD, 3);
    if do_align_and_average
        max_pad = 100;
        if numShanks == 1

            %ind_rvsl = find(stats.chDepths==stats.reversalPointDepth{1});
            if strcmp(align_to, 'top')
                ind_rvsl = find(stats.chDepths==min(stats.reversalPointDepth{1}));
            elseif strcmp(align_to, 'bottom')
                ind_rvsl = find(stats.chDepths==max(stats.reversalPointDepth{1}));
            elseif strcmp(align_to, 'nt3')
                ind_rvsl = find(stats.chDepths==max(stats.reversalPointDepth_nt3{1}));     
            else
                error('wrong input')
            end
            if ~isempty(ind_rvsl)
                CSD = stats.CSD;
                CSD = (CSD - min(min(CSD))) / ( max(max(CSD)) - min(min(CSD)) );
                pad_len = max_pad - ind_rvsl;
                disp(pad_len)
                CSD = padarray(CSD, [pad_len 0], NaN, 'pre'); % pad above
                pad_len = max_pad - (length(stats.chDepths) - ind_rvsl);
                disp(pad_len)
                CSD = padarray(CSD, [pad_len 0], NaN, 'post');
                CSDs = cat(2, CSDs, CSD);
                avgCSD = cat(3, avgCSD, CSD);
            else
                CSD = NaN(max_pad*2, size(CSD, 2));
                avgCSD = cat(3, avgCSD, CSD);
            end
        else
            for shankInd = 1:numShanks
                if strcmp(align_to, 'top')
                    ind_rvsl = find(stats.chDepths==min(stats.reversalPointDepth{shankInd}));
                elseif strcmp(align_to, 'bottom')
                    ind_rvsl = find(stats.chDepths==max(stats.reversalPointDepth{shankInd}));
                elseif strcmp(align_to, 'nt3')
                    ind_rvsl = find(stats.chDepths==max(stats.reversalPointDepth_nt3{shankInd}));
                    
                    if isempty(ind_rvsl)
                        disp('No reversal found for align. Checking if other shanks have reversal.')
                        for shankInd_ = 1:numShanks
                            ind_rvsl = find(stats.chDepths==max(stats.reversalPointDepth_nt3{shankInd_}));
                            if ~isempty(ind_rvsl)
                                break
                            end
                        end
                    end
                            
                else
                    error('wrong input')
                end
                if ~isempty(ind_rvsl)
                    CSD = stats.CSD(:,:,shankInd);
                    CSD = (CSD - min(min(CSD))) / ( max(max(CSD)) - min(min(CSD)) );
                    pad_len = max_pad - ind_rvsl;
                    disp(pad_len)
                    CSD = padarray(CSD, [pad_len 0], NaN, 'pre'); % pad above
                    pad_len = max_pad - (length(stats.chDepths) - ind_rvsl);
                    disp(pad_len)
                    CSD = padarray(CSD, [pad_len 0], NaN, 'post');
                    CSDs = cat(2, CSDs, CSD);
                    avgCSD = cat(3, avgCSD, CSD);
                else
                    CSD = NaN(max_pad*2, size(CSD, 2));
                    avgCSD = cat(3, avgCSD, CSD);
                end
            end

        end
    end


    % stats_3 = csd.getCSD(lfp, Exp, 'plotIt', false, 'noisetype', 3,'spatsmooth', 1.0, 'rvsl_window', [40 50]); 
    % csdFile = fullfile(PROCESSED_DATA_DIR,'csd',[FileTag '_csd.mat']);
    % disp('Saving CSD struct');
    % save(csdFile, 'stats')
    % fprintf('CSD struct saved to %s\n',lfpFile);
    % % fprintf('Mat File %s\n',[FileTag,'_lfp']);
    % disp('Saving CSD done');
    % for i = 1:num_shanks
    %     lower_rvsl_3 = [lower_rvsl_3; stats_3.reversalPointDepth{i}(1)];
    %     lower_rvsl_6 = [lower_rvsl_6; stats_6.reversalPointDepth{i}(1)];
    %     if length(stats_6.reversalPointDepth{i}) > 1
    %         upper_rvsl_6 = [upper_rvsl_6; stats_6.reversalPointDepth{i}(2)];
    %     else
    %         upper_rvsl_6 = [upper_rvsl_6; nan];
    %     end
    % end
    % 
    % figure(1); clf
    % csd.plotCSD(stats_6, 'overlayLFP', false)
    % title([FileTag '_noistype6'],'Interpreter','none')
    % xlabel('Time (ms)')
    % ylabel('Depth')
    % set(gcf, 'Position', get(0, 'Screensize'));
    % saveas(gca, fullfile('C:\Users\Jake\Dropbox\MarmoLabWebsite\PSA\Code\Import', 'plots', 'noisetype6', [FileTag '_noistype6']), 'png')

    if save_plots
        %try
        figure(1); clf
        if (0)
            gamma = csd.getGamma(lfp, 'method', 'weightedMin');
            csd.plotCSD(stats, 'overlayLFP', false, 'gamma', gamma)
        else
            csd.plotCSD(stats, 'overlayLFP', false)
        end
%         name = strcat(FileTag,'_noistype=',string(noisetype),'_spatial_smoothing=',string(spatial_smoothing));
%         if subtract_saccade_csd
%             name = strcat(name,'_saccade_time_avg_subtract'); 
%         end
        name = FileTag;
        title(name,'Interpreter','none')
        xlabel('Time (ms)')
        ylabel('Depth')
        set(gcf, 'Position', get(0, 'Screensize'));
        if ~exist(plots_root_folder, 'dir')
           mkdir(plots_root_folder)
        end
        saveas(gca, fullfile(plots_root_folder, char(name)), 'png')
    end
end

disp('DONE PROCESSING CSD')


if (0) % plot noisetype3 vs noisetype6 reversal 
    nt6_uppers = [];
    nt6_lowers = [];
    nt3 = [];
    for i = 1:length(CSD_structs)
        CSD_cur = CSD_structs(i);
        for j = 1:length(CSD_cur.reversalPointDepth)
            nt6_upper = min(CSD_cur.reversalPointDepth{j});
            nt6_lower = max(CSD_cur.reversalPointDepth{j});
            nt6_uppers = [nt6_uppers; nt6_upper];
            nt6_lowers = [nt6_lowers; nt6_lower];
            nt3 = [nt3; CSD_cur.reversalPointDepth_nt3{j}];
        end
    end

    figure(1); clf;
    scatter(nt3, nt6_uppers); hold on
    plot([0 1000], [0 1000])
    ylim([0 1000])
    xlim([0 1000])
    xlabel('Noisetype 3')
    ylabel('Noisetype 6 upper')
    title('Top reversal point, noisetype 6')

    figure(2); clf;
    scatter(nt3, nt6_lowers); hold on
    plot([0 1000], [0 1000])
    ylim([0 1000])
    xlim([0 1000])
    xlabel('Noisetype 3')
    ylabel('Noisetype 6 lower')
    title('Bottom reversal point, noisetype 6')
end


time = stats.time;
depth = stats.depth;


% %%%%%%######## SAVE CSDs for Averaging
% save_name = 'CSD_345.m';
% disp('Saving relevant variables for plotting');
% File_path = 'Z:\Data\Processed_Laminar\misc_files';
% save_file = fullfile(File_path, save_name);
% save(save_file, '-v7.3', 'avgCSD', 'CSDs', 'max_pad', 'depth', 'time', 'deadChan')
% fprintf('struct saved to %s\n',save_file);
% disp('Saving done');%
% %%%%%%%%############%%%%

if do_align_and_average
    num_sessions = size(avgCSD, 3);
    time_points = size(avgCSD,2);
    
%     inds_examples = 3:6;
%     CSDs2 = CSDs
    
    figure(5); clf; 
    rm_nan_rows = ~all(isnan(CSDs),2);
%     offset = find((rm_nan_rows==0));
%     offset = offset(offset>max_pad);
%     offset = length(offset);
    offset = find((rm_nan_rows==0));
    offset = offset(offset<max_pad);
    offset = length(offset);
    CSDs2 = CSDs(rm_nan_rows, :);
    CSDs2(isnan(CSDs2)) = -0.1;
    coords = flip(1:size(CSDs2, 1));
    med_diff = median(diff(depth));
    channel0 = max_pad-offset;
    yticks = [];
    for i=1:length(coords)
        yticks = [yticks (i-channel0)*med_diff];
    end
    csd_c = time_points/2; %(window(2) - window(1))/2;
    x_ticks = [];
    line_indices = [];
    for i=1:num_sessions
       x_ticks = [x_ticks csd_c+(i-1)*time_points];
       line_indices = [line_indices (i-1)*time_points];
    end
    imagesc(1:size(CSDs2, 2), yticks, CSDs2); axis ij
    xticks(x_ticks)
    xticklabels(1:num_sessions)
    cmap = parula;
    cmap = [ 1 1 1 ; cmap ];
    colormap(cmap);
    hold on
    xLimits = get(gca,'XLim');
%     plot([floor(xLimits(1)) floor(xLimits(2))], [1; 1]*max_pad+1-offset, 'r--', 'Linewidth', 2);
%    plot([floor(xLimits(1)) floor(xLimits(2))], [1; 1]*0, 'r--', 'Linewidth', 2);
    yLimits = get(gca,'YLim');
    for i=line_indices
        plot([1; 1]*i,[yLimits(1) yLimits(2)],'color','white', 'Linewidth', 2);
    end
    xlabel('Example Session')
    ylabel('Depth (microns)')
    disp('done')
    
    if save_plots
        title('Aligned CSDs','Interpreter','none')
        set(gcf, 'Position', get(0, 'Screensize'));
        if ~exist(plots_root_folder, 'dir')
           mkdir(plots_root_folder)
        end
        saveas(gca, fullfile(plots_root_folder, 'Aligned_CSDs'), 'png')
    end
        

    avgCSD = mean(avgCSD,3,'omitnan');
    figure(6); clf; 
    rm_nan_rows = ~all(isnan(avgCSD),2);
    offset = find((rm_nan_rows==0));
    offset = offset(offset<max_pad);
    offset = length(offset);
    avgCSD2 = avgCSD(rm_nan_rows, :);
    %avgCSD2 = avgCSD2(3:end-3,:);
    coords = flip(1:size(avgCSD2, 1));
    med_diff = median(diff(depth));
    channel0 = max_pad-offset;
    yticks = [];
    for i=1:length(coords)
        yticks = [yticks (i-channel0)*med_diff];
    end
    imagesc(time, yticks, avgCSD2); axis ij
%     cmap = parula;
%     cmap = [ 1 1 1 ; cmap ];
    colormap(parula);
    hold on
    xLimits = get(gca,'XLim');
%     plot([floor(xLimits(1)) floor(xLimits(2))], [1; 1]*max_pad+1-offset, 'r--', 'Linewidth', 2);
%    plot([floor(xLimits(1)) floor(xLimits(2))], [1; 1]*0, 'r--', 'Linewidth', 2);
%     title('Average MUA')
    xlabel('time (ms)')
    ylabel('Depth (microns)')
    disp('done')
    
    if save_plots
        title('Average CSDs','Interpreter','none')
        set(gcf, 'Position', get(0, 'Screensize'));
        if ~exist(plots_root_folder, 'dir')
           mkdir(plots_root_folder)
        end
        saveas(gca, fullfile(plots_root_folder, 'Average_CSD'), 'png')
    end

end
% disp('hi')
% 
% figure(1); clf
% plot(lower_rvsl_3, lower_rvsl_6, '*', 'Color', 'red'); hold on
% hline = refline([1 0]);
% hline.Color = 'b';
% title('reversal noisetype 3 vs 6','Interpreter','none')
% xlabel('Noisetype 3 reversal depth')
% ylabel('Noisetype 3 reversal depth')
% 
% figure(2); clf
% rows_without_nan = ~isnan(upper_rvsl_6);
% lower_rvsl_6_nonan = lower_rvsl_6(rows_without_nan);
% upper_rvsl_6_nonan = upper_rvsl_6(rows_without_nan);
% input_heights = lower_rvsl_6_nonan - upper_rvsl_6_nonan;
% histogram(input_heights, 20)
% title('Input height from noisetype 6','Interpreter','none')
% xlabel('Input height from noisetype 6')
% ylabel('counts')
% if do_align_and_average
% figure(1); clf; 
% rm_nan_rows = ~all(isnan(CSDs),2);
% offset = find((rm_nan_rows==0));
% offset = offset(offset>max_pad/2);
% offset = length(offset);
% CSDs2 = CSDs(rm_nan_rows, :);
% coords = flip(1:size(CSDs2, 1));
% imagesc(1:size(CSDs2, 2), coords, CSDs2); axis ij
% colormap(parula);
% hold on
% xLimits = get(gca,'XLim');
% plot([floor(xLimits(1)) floor(xLimits(2))], [1; 1]*max_pad+1-offset, 'r--', 'Linewidth', 2);
% disp('done')
% 
% avgCSD = mean(avgCSD,3,'omitnan');
% figure(2); clf; 
% rm_nan_rows = ~all(isnan(avgCSD),2);
% offset = find((rm_nan_rows==0));
% offset = offset(offset>max_pad/2);
% offset = length(offset);
% avgCSD2 = avgCSD(rm_nan_rows, :);
% coords = flip(1:size(avgCSD2, 1));
% imagesc(stats.time, coords, avgCSD2); axis ij
% colormap(parula);
% hold on
% xLimits = get(gca,'XLim');
% plot([floor(xLimits(1)) floor(xLimits(2))], [1; 1]*max_pad+1-offset, 'r--', 'Linewidth', 2);
% title('Average CSD')
% xlabel('time (ms)')
% disp('done')
% 
% end

% %% 6 PLOT CSD - noisetype 6 (moving flashes)
% stats = csd.getCSD(lfp, Exp, 'plotIt', false, 'noisetype', 6);
% figure(2); clf
% csd.plotCSD(stats, 'overlayLFP', false)
% title([FileTag '_noistype6'],'Interpreter','none')
% xlabel('Time (ms)')
% ylabel('Depth')

% %% at the end clear it up
% %*********** clear up environment when finished
% clear all;
% close all;







end

function i = findZeroCrossings(data, mode)
%FINDZEROCROSSINGS Find zero crossing points.
%   I = FINDZEROCROSSINGS(DATA,MODE) returns the indicies into the supplied
%   DATA vector, corresponding to the zero crossings.
%
%   MODE specifies the type of crossing required:
%     MODE < 0 - results in indicies for the -ve going zero crossings,
%     MODE = 0 - results in indicies for ALL zero crossings (default), and
%     MODE > 0 - results in indicies for the +ve going zero crossings.

% $Id: findZeroCrossings.m,v 1.1 2008-07-21 23:31:50 shaunc Exp $

if nargin < 2
    mode = 0;
end

[i,~,p] = find(data); % ignore zeros in the data vector

switch sign(mode)
    case -1
        % find -ve going crossings
        ii = find(diff(sign(p))==-2);
    case 0
        % find all zero crossings
        ii = find(abs(diff(sign(p)))==2);
    case 1
        % find +ve going crossings
        ii = find(diff(sign(p))==2);
end;

i = round((i(ii)+i(ii+1))/2);
end

