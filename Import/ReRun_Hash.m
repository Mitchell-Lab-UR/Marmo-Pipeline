%Batch A Files for New Hash Processing

INFO_DATA_DIR = 'Z:\Users\Amy\Processed_A'; 
SERVER_DATA_DIR = 'Z:\Data\PSA';
PROCESSED_DATA_DIR = 'Z:\Users\Amy\Processed_Hash_thres'; 
InfoNames = dir([INFO_DATA_DIR]); %Find all processed A files
N = length(InfoNames)-2;

%need to do i=1;
for i = 3:N
    filepath = [INFO_DATA_DIR,filesep,InfoNames(i+2).name];
    BaseTag = InfoNames(i+2).name;
    load(filepath);
    disp(sprintf('Loading Exp: %s',BaseTag));
    DataFolder = [SERVER_DATA_DIR,filesep,Exp.FileTag];
    FileTag = Exp.FileTag; 


    %%  compute the CSD onsets if available
    CSD_Spike_Onset_Script;  % determine onsets of CSD task
    minCSD = min(CSD.Onsets_Ephys)-0.5;  % first CSD flash in secs
    maxCSD = max(CSD.Offsets_Ephys)+0.5; % last CSD flash in sec
    Exp.CSD = CSD;

    %% HASH sort all the channels and store a struct
    [hash,smallfp] = compute_hash_channels_interval2(DataFolder,[],[minCSD,maxCSD]); %does not use one threshold
    %[hash,smallfp] =
    %compute_hash_channels_interval(DataFolder,[],[minCSD,maxCSD]); %AB use
    %single threshold (20Hz)
    hashFile = [PROCESSED_DATA_DIR,filesep,FileTag,'_hashint.mat'];
    save(hashFile,'-v7.3','hash');
    fprintf('Hash struct saved to %s\n',hashFile);
    smallfpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'_smallfp.mat'];
    save(smallfpFile,'-v7.3','smallfp');
    fprintf('Smallfp struct saved to %s\n',smallfpFile);
    fprintf('Mat File %s\n',[FileTag,'_smallfp']);
    disp('Saving Hash and Smallfp done');
end