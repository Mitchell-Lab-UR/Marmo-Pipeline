%% Add paths
addMarmoPipe()

%% Specify user 

user = 'bluethunder';

switch user
    case 'gabehome'
        SERVER_DATA_DIR = ''; % have no server access at home
        PROCESSED_DATA_DIR = '/Users/gabrielsarch/Documents/MitchellLabRepo/Data'; %C:\Processed';
    case 'judehome'
        SERVER_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Raw';
        PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Processed';
    case 'bluethunder'
         %SERVER_DATA_DIR = 'Y:\Data\PSA';
         %PROCESSED_DATA_DIR = 'Y:\Data\Processed_Laminar'; %C:\Processed';
           
         %SERVER_DATA_DIR = 'C:\PSA_Gravedigger\Raw';  
         SERVER_DATA_DIR = 'Z:\Data\PSA'; %pull from the server  
         PROCESSED_DATA_DIR = 'Z:\Data\Processed_Laminar'; %'Z:\Data\PSA';
        
%         SERVER_DATA_DIR = 'C:\PSA_Gravedigger\Raw';
%         PROCESSED_DATA_DIR = 'E:\PSA_Gravedigger\Processed';
         
%          SERVER_DATA_DIR = 'Z:\Data\AllenV1';
%          PROCESSED_DATA_DIR = 'Z:\Data\AllenV1\Processed';
         
         %SERVER_DATA_DIR = 'C:\PSA_Gravedigger\Raw';
         %PROCESSED_DATA_DIR = 'E:\PSA_Gravedigger\Processed'; %C:\Processed';
        % PROCESSED_DATA_DIR = 'C:\Raw\Gabe_to_download_single';
    case 'jakegravedigger'
        SERVER_DATA_DIR = 'C:\Raw';
        PROCESSED_DATA_DIR = 'Z:\PSA_EDF\EDF_Processed_Decoding\X-Y Recordings';
        
end


%% Example 1 - Process raw LFPs (Note: need to get processed LFPs before running csd.processCSD)

raw_folder_names = {'Milo_2022-02-22_12-24-27_MT64_72b'};
FileTags = {'Milo_220222'}; % Filetags corresponding to each raw_folder_names (these should also correspond to experiment file tags)
csd.processLFPs(raw_folder_names, SERVER_DATA_DIR, PROCESSED_DATA_DIR, 'FileTags', FileTags, 'override_existing_lfp', false);

%% Example 2 - Noisetype 3 

% Specify session names by FileTag (name of experiment file tags)
%Sessions = {'Milo_220221', 'Milo_170321', 'Milo_180321', 'Milo_070421', 'Milo_120421', 'Milo_160421', 'Milo_210421', 'Milo_100521', 'Milo_140521', 'Milo_080721', 'Milo_130721', 'Milo_150721', 'Milo_280721','Milo_300721','Milo_170821','Milo_240821','Milo_260821', 'Milo_160921'};
%Sessions = {'Milo_070421'}; % FileTags of sessions

%Batch MTC CSDs (ALL FILES)
%Sessions = {'Milo_300921','Milo_051021', 'Milo_071021',
%'Milo_121021', 'Milo_141021', 'Milo_211021', 'Milo_091121', 'Milo_091221',
%'Milo_141221', 'Milo_180222', 'Milo_220222', 'Milo_220322',
%'Milo_290322', 'Milo_010422', 'Milo_080422', 'Milo_150422',
%'Milo_190422'}

%Sessions = {'Milo_300921','Milo_051021', 'Milo_071021','Milo_121021', 'Milo_141021', 'Milo_091121', 'Milo_091221', 'Milo_141221', 'Milo_180222','Milo_220222', 'Milo_290322', 'Milo_220322', 'Milo_010422', 'Milo_080422', 'Milo_150422'};

%Sessions = {'Milo_301121', 'Milo_060122', 'Milo_280122', 'Sprout_161222', 'Sprout_270123', 'Sprout_030323'};


%Sessions = {'Sprout_030323'}; No lfp??

plots_folder = './plots_3';

%MT Examples
Sessions = {'Milo_301121', 'Milo_060122', 'Milo_080222', 'Sprout_161222', 'Sprout_010223', 'Sprout_010323'};
 
rvsl_manual_FileTags = {'Milo_301121', 'Milo_301121', 'Milo_060122', 'Milo_060122', 'Milo_080222', 'Milo_080222', 'Sprout_161222', 'Sprout_161222', 'Sprout_010223', 'Sprout_010223' ,'Sprout_010323', 'Sprout_010323'}; % filetags to alter
 rvsl_manual_location = {'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top'}; % top or bottom reversal point for each filetag
 rvsl_manual_shank_index = [1 2 1 2 1 2 1 2 1 2 1 2]; % which shank index to alter
 rvsl_manual_shank_depth = [525 525 735 780 550 620 560 560 700 700 525 805]; % which shank index to alter
%  
%MTC Examples
%  Sessions = {'Milo_180222', 'Milo_190422', 'Sprout_201022', 'Sprout_071222', 'Sprout_031122'};
%  rvsl_manual_FileTags = {'Milo_180222', 'Milo_180222', 'Milo_190422', 'Milo_190422', 'Sprout_201022', 'Sprout_201022', 'Sprout_071222', 'Sprout_071222', 'Sprout_031122', 'Sprout_031122'}; % filetags to alter
%  rvsl_manual_location = {'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top'}; % top or bottom reversal point for each filetag
%  rvsl_manual_shank_index = [1 2 1 2 1 2 1 2 1 2 ]; % which shank index to alter
%  rvsl_manual_shank_depth = [175 290 245 350 525 560 490 525 560 595 ]; % which shank index to alter

% manual alter reversal depths on some files 
% rvsl_manual_FileTags = {'Milo_010422', 'Milo_080422', 'Milo_091221', 'Milo_121021', 'Milo_121021', 'Milo_141021','Milo_141221','Milo_180222'}; % filetags to alter
% rvsl_manual_location = {'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top'}; % top or bottom reversal point for each filetag
% rvsl_manual_shank_index = [2 2 1 1 2 1 1 2]; % which shank index to alter
% rvsl_manual_shank_depth = [275 220 875 600 520 540 850 100]; % which shank index to alter

CSDs = csd.processCSD(Sessions, PROCESSED_DATA_DIR, ...
                        'noisetype', 'saccade3', ... % specify noisetype (3, 6, saccade3, saccade6)
                        'spatial_smoothing', 1.0, ... % how much spatial smoothing (unit = 1 channel)
                        'temporal_smoothing', 1.0, ... % how much temporal smoothing (unit = 1 time point)
                        'subtract_pre_csd', true, ... % subtract off time-averaged baseline CSD (average from subtract_pre_CSD_time0 to subtract_pre_CSD_time1 relative to onset)
                        'subtract_pre_CSD_time0', -100, ... 
                        'subtract_pre_CSD_time0', 0, ... 
                        'do_align_and_average', true, ... % Plot aligned and averaged CSD across sessions
                        'get_reversal_from_average_CSD', false, ... % If there are multiple shanks, then average across the two shanks to get same reversal for both (false for this example) 
                        'align_to', 'top', ... % %'nt3' What reversal point to align to? 
                        'save_plots', true, ... % Save out plots in plots_root_folder folder
                        'plots_root_folder', plots_folder,... % folder to save plots to
                        'rvsl_manual_FileTags', rvsl_manual_FileTags, ...
                        'rvsl_manual_location', rvsl_manual_location, ...
                        'rvsl_manual_shank_index', rvsl_manual_shank_index, ...
                        'rvsl_manual_shank_depth', rvsl_manual_shank_depth);
                    
%% Example 3 - Noisetype 6 aligned to top reversal  

% Specify session names by FileTag (name of experiment file tags)
% Sessions = {'Milo_220221', 'Milo_170321', 'Milo_180321', 'Milo_070421', 'Milo_120421', 'Milo_160421', 'Milo_210421', 'Milo_100521', 'Milo_140521', 'Milo_080721', 'Milo_130721', 'Milo_150721', 'Milo_280721','Milo_300721','Milo_170821','Milo_240821','Milo_260821', 'Milo_160921'};
%Sessions = {'Milo_070421'};
%Sessions = {'Milo_091121'};

%Sessions = {'Milo_300921','Milo_051021', 'Milo_071021','Milo_121021', 'Milo_141021', 'Milo_091121', 'Milo_091221', 'Milo_141221', 'Milo_180222','Milo_220222', 'Milo_290322', 'Milo_220322', 'Milo_010422', 'Milo_080422', 'Milo_150422'};
%Sessions = {'Milo_010422'};

%Sessions= {'Milo_010422', 'Milo_051021', 'Milo_071021','Milo_080422', 'Milo_091121', 'Milo_091221', 'Milo_121021', 'Milo_141021', 'Milo_141221', 'Milo_180222', 'Milo_220222', 'Milo_220322', 'Milo_290322', 'Milo_300921'};

%Dyed Recordings
%Sessions= {'Milo_030622', 'Milo_170622', 'Milo_070722'};

%Sessions= {'Milo_301121', 'Milo_180322'};
%Sessions = {'Sprout_251022'};
%091222 141222 161222 281222 040123 060123 Sprout_270123c
Sessions = {'Milo_280421'};

plots_folder = './plots_6';

% manual alter reversal depths on some files 
% rvsl_manual_FileTags =  {'Milo_010422', 'Milo_051021', 'Milo_071021','Milo_080422', 'Milo_091121', 'Milo_091221', 'Milo_121021', 'Milo_141021', 'Milo_141221', 'Milo_180222', 'Milo_220222', 'Milo_220322', 'Milo_290322', 'Milo_300921', 'Milo_010422', 'Milo_051021', 'Milo_071021','Milo_080422', 'Milo_091121', 'Milo_091221', 'Milo_121021', 'Milo_141021', 'Milo_141221', 'Milo_180222', 'Milo_220222', 'Milo_220322', 'Milo_290322', 'Milo_300921'}; %, 'Milo_051021', 'Milo_071021','Milo_080422', 'Milo_091121', 'Milo_091221', 'Milo_121021', 'Milo_141021', 'Milo_141221', 'Milo_180222', 'Milo_220222', 'Milo_220322', 'Milo_290322', 'Milo_300921'};; % filetags to alter
% rvsl_manual_location = {'top', 'top', 'top','top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top', 'top'}; % top or bottom reversal point for each filetag
% rvsl_manual_shank_index = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2]; % which shank index to alter
% rvsl_manual_shank_depth = [NaN 455 315 475 315 560 385 490 630 175 230 500 500 350 100 455 315 NaN 445 595 385 490 595 275 245 580 350 350]; % which shank index to alter

CSDs = csd.processCSD(Sessions, PROCESSED_DATA_DIR, ...
                        'noisetype', '6', ... % specify noisetype (3, 6, saccade3, saccade6)
                        'spatial_smoothing', 1.0, ... % how much spatial smoothing (unit = 1 channel)
                        'temporal_smoothing', 1.0, ... % how much temporal smoothing (unit = 1 time point)
                        'subtract_pre_csd', true, ... % subtract off time-averaged baseline CSD (average from subtract_pre_CSD_time0 to subtract_pre_CSD_time1 relative to onset)
                        'subtract_pre_CSD_time0', -100, ... 
                        'subtract_pre_CSD_time0', 0, ... 
                        'do_align_and_average', false, ... % Plot aligned and averaged CSD across sessions
                        'get_reversal_from_average_CSD', false, ... % If there are multiple shanks, then average across the two shanks to get same reversal for both (false for this example) 
                        'align_to', 'top', ... % What reversal point to align to? (can be 'top' or 'bottom' or 'nt3')
                        'save_plots', true, ... % Save out plots in plots_root_folder folder
                        'plots_root_folder', plots_folder); % folder to save plots to
%                         'rvsl_manual_FileTags', rvsl_manual_FileTags, ...
%                         'rvsl_manual_location', rvsl_manual_location, ...
%                         'rvsl_manual_shank_index', rvsl_manual_shank_index, ...
%                         'rvsl_manual_shank_depth', rvsl_manual_shank_depth);
                    
                    
                    
                    
                    