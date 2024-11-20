%*********** ANALYSIS SCRIPT *******************************
%***
%*** SPECIAL IMPORT TO FIX THE EVENTS FILE WHEN 32bit missing
%***
%***********************************************************
%*** NOTE:  The goal here is just to analyze a single file
%***        This will not be very well organized
%***        Once we get an initial analysis done, we will
%***        build a GUI to do this systematically and
%***        provide a larger structure
%***       
%***        For now, all info is imported to an Exp struct
%***        which is saved back to the data folder and then
%***        can be loaded for any analysis
%***
%***        Revised so it concatenates all files from a session
%***
%***        Things still to do:
%***        At present, we are still not:
%***          1) storing information about electrode type
%***          2) storing information about LFP (notch filtered)
%***          3) using Kilosort for multi-channel data sorting
%***          4) storing session info to Meta-Table
%************************************************************

%% Spike File importing

%FileTag = 'Allen_2022-03-04_12-14-43_V1_64a'
FileTag = 'Allen_2022-08-05_13-24-54_V1_64b';

%************
TAGSTART = 63;   
TAGEND = 62; 
%*************
athome = 0;
if (athome == 1)
   SERVER_DATA_DIR = 'C:\Users\jmitchell\Box Sync\FovealTrans\V1Hart\Import\Raw';
   PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Box Sync\FovealTrans\V1Hart\Import\Processed';
else      
   SERVER_DATA_DIR = 'C:\PSA_Gravedigger\Raw';
   PROCESSED_DATA_DIR = 'E:\PSA_Gravedigger\Processed';
end
DataFolder = [SERVER_DATA_DIR,filesep,FileTag];


%% Events File importing
%******** now grab the events file with strobes
EventFiles = dir([DataFolder,filesep,'*.events']);
if isempty(EventFiles)
    disp('Error finding events file');
    return;
end
[evdata,evtime,evinfo] = read_ephys.load_open_ephys_data([DataFolder,filesep,EventFiles(1).name]);
%**** convert events into strobes with times
[tstrobes,strobes] = read_ephys.convert_data_to_strobes(evdata,evtime,evinfo);
% strobes = read_ephys.fix_missing_strobes(strobes);  % should not apply
disp('Stobes are loaded');


%% Loading up the MarmoView Data Files
%****************************************************************
ExpFiles = dir([DataFolder,filesep,'*z.mat']);
if isempty(ExpFiles)
    disp('Error finding *z.mat file');
    return;
end
%****** get order by date of files to import them
FileDates = cellfun(@datenum,{ExpFiles(:).date});
DoFileSort = [ (1:size(FileDates,2))' FileDates'];
FileSort = sortrows(DoFileSort,2); % sort them by date
%***** read and append files in order by date *********
BigN = 0;
for zk = FileSort(:,1)'
  fname = ExpFiles(zk).name;  
  load([DataFolder,filesep,fname]);
  if ~BigN
    Exp.D = D;
    Exp.S = S;
    BigN = size(D,1); % number of trials
  else
    for k = 1:size(D,1)
       Exp.D{BigN+k} = D{k};  % appending trials 
    end  
    BigN = BigN + size(D,1);
  end
  clear D;
  clear S;
  disp(sprintf('Experiment file %s loaded',fname));
end
%***** store spikes info in Exp struct, let's keep all info there
%***** once finished we will clear everything but the Exp struct
disp('Experiment files loaded');
%**************************

%% Synching up strobes from Ephys to MarmoView
%******************* returns start and end times in ephys record, or NaN if missing
%***** this code would look simple, but sometimes one bit is flipped in
%***** error and then you have to play catch up to find the missing code
%***** since there is some redundancy (start and end codes for each trial)
%***** this gives us a way to recover cases with just one errant bit
disp('Synching up ephys strobes to fix 32bit issue');
%***** in new code, must specify year, month, and day when calling find
%***** the find strobe routine
codelist = [];  % list of matched start and end indices, use order to repair ambiguous ones
%*********  step through Exp.D strobes
for k = 1:size(Exp.D,1)
  for tagger = 1:2  % start or end
     if (tagger == 1)
        TAG = TAGSTART-32;
        taglet = Exp.D{k}.STARTCLOCK;
        ktag = k;
     else
        TAG = TAGEND-32;
        taglet = Exp.D{k}.ENDCLOCK;
        ktag = -k;  % use negatives to indicate end tags in list
     end 
     %**** transform for -32 on the ctaglet
     ctag = taglet;
     for ck = 1:length(ctag)
        if (ctag(ck) >= 32)
          ctag(ck) = ctag(ck)-32;
        end
     end
     %***** find all possible matches (hopefully only one)
     zz1 = find( (strobes == TAG) );
     zz2 = find( (strobes(zz1+1) == ctag(1)) & (strobes(zz1+2) == ctag(2)) & ...
                 (strobes(zz1+3) == ctag(3)) & (strobes(zz1+4) == ctag(4)) & ...
                 (strobes(zz1+5) == ctag(5)) & (strobes(zz1+6) == ctag(6)) );
     zz = zz1(zz2);  % only those fitting all day codes should accept        
     %******** exact match?  then fix the strobe, else skip it for now
     if (size(zz) == 1)
        tip = [ktag zz(1) 0 0 0];
     else
        %******
        ln = length(zz);
        tip = [ktag zz'];
        for zk = 1:(4-ln)
           tip = [tip 0];
        end
     end
     codelist = [codelist ; tip];
  end
end
disp('Finished matching up unique strobes, now return to order non-uniques'); 
codelist
input('check');

%% then go back to fill in matches and fix ambiguous ones if possible
lastmatch = [];
for k = 1:size(codelist,1) 
   %****** look for start
   zz = codelist(k,2:5);
   if (zz(1) > 0) && (zz(2) == 0)   % unique match
         lastmatch = zz(1);
   else
       if (zz(1) > 0) && (zz(2) > 0)   % ambiguous match, then search
           if ~isempty(lastmatch)
               %*** pick label closest to the lastmatch if available
               dist = (zz-lastmatch);
               dist( find(dist < 0) ) = NaN;
               mindist = nanmin(dist);
               bz = find( dist == mindist);
               %*** swap out the startlist then
               zz = zz([bz,1:(bz-1),(bz+1):4]);
               codelist(k,2:5) = zz;  % swap in right order
               %***** since you found it, now this is lastmatch
               lastmatch = zz(1);
               %*****
           else
               % search forward for next match, select closest below that
               formatch = [];
               for fk = (k+1):size(codelist,1)
                   fzz = codelist(fk,2:5);
                   if (fzz(1) > 0) && (fzz(2) == 0)   % unique match
                       formatch = fzz(1);
                       break;
                   end
               end
               if ~isempty(formatch)
                  dist = (formatch-zz);
                  dist( find(dist < 0) ) = NaN;
                  mindist = nanmin(dist);
                  bz = find( dist == mindist);
                  %*** swap out the startlist then
                  zz = zz([bz,1:(bz-1),(bz+1):4]);
                  codelist(k,2:5) = zz;  % swap in right order
                  %***** since you found it, now this is lastmatch
                  lastmatch = zz(1);
               else
                  % unable to resolve the ambiguity, set all to zero
                  codelist(k,:) = 0;  % giving up on this code
                  lastmatch = [];  % all best are off, start again
               end
           end
       end
   end
   %********
end
codelist
input('check');

%% At the final step fill in the match in first position if non-zero
fixstrobes = strobes;  % replace strobes with repair
for k = 1:size(codelist,1)
   if (codelist(k,1) > 0)  % indicates a start code
      kk = codelist(k,1);
      taglet = Exp.D{kk}.STARTCLOCK;
      zz = codelist(k,2:5);
      if (zz(1) > 0) 
         fixstrobes(zz(1)+0) = TAGSTART;
         fixstrobes(zz(1)+(1:6)) = taglet;
      end
   else
      kk = -codelist(k,1);
      taglet = Exp.D{kk}.ENDCLOCK;
      zz = codelist(k,2:5);
      if (zz(1) > 0) 
         fixstrobes(zz(1)+0) = TAGEND;
         fixstrobes(zz(1)+(1:6)) = taglet;
      end       
   end
end
[strobes fixstrobes]
input('check');

%% finally swap in corrected strobes
strobes = fixstrobes;

%% Save Exp Struct to data folder for future use (PROC folder)
%****** should consider later how we want to store things (per protocol,
%****** but off hand it seems feasible given the generic structure of
%****** of the D struct maybe we could concatenate all sessions in one
ExpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'_strobes.mat'];
Exp.ProcDataFolder = PROCESSED_DATA_DIR;
Exp.DataFolder = DataFolder;
Exp.FileTag = FileTag;
save(ExpFile,'strobes','tstrobes');
disp(sprintf('Strobes and Tstrobes struct saved to %s',ExpFile));
disp(sprintf('Mat File %s',[FileTag,'_strobes']));

%% at the end clear it up
%*********** clear up environment when finished
clear all;
