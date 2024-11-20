function processLFPs(raw_folder_names,server_data_dir, processed_data_dir,varargin)
    % processCSD interacts with the +csd toolbox to obtain and plot CSDs
    % Inputs:
    %   raw_folder_names        list of folder names in server_data_dir
    %                           with raw LFP to process
    %   server_data_dir         string - path to server files where
    %                           raw_folder_names are located
    %   processed_data_dir      Location of exp and saved lfp processed files. Exp
    %                           files should be directly in processed_data_dir
    %                           as {FileTag}.mat while LFP processed files should
    %                           be in processed_data_dir/lfp as
    %                           {FileTag}_lfp.mat
    % Arguments: 
    %   override_existing_lfp   TRUE/FALSE - whether to override existing
    %                           LFP if it exists
    
    ip = inputParser();
    ip.addParameter('override_existing_lfp', false)
    ip.addParameter('FileTags', {})
    ip.parse(varargin{:});

    override_existing_lfp = ip.Results.override_existing_lfp;
    FileTags = ip.Results.FileTags;
    
    
    for sess_idx = 1:length(raw_folder_names)
        
        DataFolder= fullfile(server_data_dir, raw_folder_names{sess_idx});
        %end
        ops = io.loadOps(DataFolder);
        % [fPath, fName, fExt] = fileparts(ops.root);
        % newStr = split(ops.root,'/');
        newStr=regexp(ops.root,'\','split');
        % if load_from_server
        %     new_root = fullfile(ROOT_RAW_DIR,newStr{end-1},newStr{end});
        % else
        %     new_root = fullfile(SERVER_DATA_DIR,newStr{end-1},newStr{end});
        % end

        new_root = fullfile(server_data_dir,newStr{end-1},newStr{end});
        ops.root = new_root;
        disp(new_root)

        newStr=regexp(ops.fbinary,'\','split');
        new_fbinary = fullfile(new_root,newStr{end});
        ops.fbinary = new_fbinary;
        disp(new_fbinary)
        
        if isempty(FileTags)
            disp('NOTE: FileTags not given.. getting FileTag from raw_folder_names')
            disp('NOTE: getting FileTags this way may cause issues if experiment file is named different')
            tag_split = split(raw_folder_names{sess_idx}, '_');
            subj_name = tag_split(1);
            date = split(tag_split(2), '-');
            yr = split(date(1), '0');
            yr = yr(2);
            day = date(3);
            month = date(2);
            FileTag = strcat(subj_name, '_', day, month, yr);
            FileTag = FileTag{1};
        else
            FileTag = FileTags{sess_idx};
        end

        disp('loading exp file')
        ExpFile = [processed_data_dir,filesep,FileTag,'.mat'];
        load(ExpFile);
        disp('done')
        
        Exp
        %input('check');

        disp('loading lfp...')
        lfpFile = fullfile(processed_data_dir,'lfp',[FileTag '_lfp.mat']);

        size(lfpFile)
        %input('check');
        
        if ~exist(lfpFile, 'file') || override_existing_lfp

            disp('LFP processed file not found at given file location.. generating it (this will take a while)..')
            [data, timestamps, info] = io.getLFP(ops, true, false);
            unique_x = unique(Exp.osp.xcoords);
            num_shanks = size(unique_x, 1);
            shank_len = size(data, 2)/num_shanks;
            deadChan = []; % Add dead channels here if there are any (1-32 first shank, 33-64 second shank, etc.)
            if num_shanks > 0

                xcoords = ones(shank_len,num_shanks);
                for i = 1:num_shanks
                    xcoords(:,i) = xcoords(:,i)*unique_x(i);
                end
            end
            diff_ycoords = abs(mode(diff(Exp.osp.ycoords)));
            ycoords = repmat(flip(linspace(diff_ycoords, diff_ycoords*shank_len, shank_len)), num_shanks, 1)';

            %lfpFile = fullfile(PROCESSED_DATA_DIR,'lfp',[FileTag '_lfp.mat']);
            % lfpFile = [PROCESSED_DATA_DIR,filesep,'lfp',filesep,FileTag,'_lfp.mat'];
            disp('Saving LFP struct');
            save(lfpFile, '-v7.3', 'timestamps', 'info', 'data', 'xcoords', 'ycoords', 'deadChan')
            fprintf('LFP struct saved to %s\n',lfpFile);
            fprintf('Mat File %s\n',[FileTag,'_lfp']);
            disp('Saving LFP done');
        elseif exist(lfpFile, 'file')
            error('LFP file already exists... either set new path or set override_existing_lfp=true to override it')
        end
    end
end

